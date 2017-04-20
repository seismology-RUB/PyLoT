#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyLoT: Main program
===================
PyLoT is a seismic data processing software capable of picking seismic
phases (symmetric and asymmetric error assignment), exporting these to
several common phase data formats and post process the data, e.g. locating
events, via external localization software.
Additionally PyLoT is meant as an interface to autoPyLoT which can
automatically pick seismic phases, if the parameters have properly been
chosen for the particular data set.

Some icons are out of a free of charge icon set, which can be found here:
https://www.iconfinder.com/iconsets/flavour

:authors:
    Sebastian Wehling-Benatelli / Ludger Küperkoch
:copyright:
    The PyLoT Development Team (https://ariadne.geophysik.rub.de/trac/PyLoT)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""

import os
import sys
import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from PySide import QtGui, QtCore
from PySide.QtCore import QCoreApplication, QSettings, Signal, QFile, \
    QFileInfo, Qt, QSize
from PySide.QtGui import QMainWindow, QInputDialog, QIcon, QFileDialog, \
    QWidget, QHBoxLayout, QVBoxLayout, QStyle, QKeySequence, QLabel, QFrame, QAction, \
    QDialog, QErrorMessage, QApplication, QPixmap, QMessageBox, QSplashScreen, \
    QActionGroup, QListWidget, QDockWidget, QLineEdit, QListView, QAbstractItemView, \
    QTreeView, QComboBox, QTabWidget, QPushButton, QGridLayout
import numpy as np
from obspy import UTCDateTime

from pylot.core.analysis.magnitude import RichterMagnitude, MomentMagnitude
from pylot.core.io.data import Data
from pylot.core.io.inputs import FilterOptions, AutoPickParameter
from autoPyLoT import autoPyLoT
from pylot.core.pick.compare import Comparison
from pylot.core.pick.utils import symmetrize_error
from pylot.core.io.phases import picksdict_from_picks
import pylot.core.loc.nll as nll
from pylot.core.util.defaults import FILTERDEFAULTS, COMPNAME_MAP
from pylot.core.util.errors import FormatError, DatastructureError, \
    OverwriteError, ProcessingError
from pylot.core.util.connection import checkurl
from pylot.core.util.dataprocessing import read_metadata, restitute_data
from pylot.core.util.utils import fnConstructor, getLogin, \
    full_range
from pylot.core.io.location import create_creation_info, create_event
from pylot.core.util.widgets import FilterOptionsDialog, NewEventDlg, \
    WaveformWidget, PropertiesDlg, HelpForm, createAction, PickDlg, \
    getDataType, ComparisonDialog
from pylot.core.util.map_projection import map_projection
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.thread import AutoPickThread, Thread
from pylot.core.util.version import get_git_version as _getVersionString
import icons_rc

locateTool = dict(nll=nll)


class MainWindow(QMainWindow):
    __version__ = _getVersionString()
    closing = Signal()

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        # check for default pylot.in-file
        infile = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
        if os.path.isfile(infile)== False:
            infile = QFileDialog().getOpenFileName(caption='Choose PyLoT-input file') 
            if not os.path.exists(infile[0]):
                QMessageBox.warning(self, "PyLoT Warning",
                           "No PyLoT-input file declared!")
                sys.exit(0)
            self.infile = infile[0]
        else:
             self.infile = infile

        self.project = Project()
        self.array_map = None
        self._metadata = None
        self._eventChanged = False

        self.poS_id = None
        self.ae_id = None

        # UI has to be set up before(!) children widgets are about to show up
        self.createAction = createAction
        # read settings
        settings = QSettings()
        self.recentfiles = settings.value("data/recentEvents", [])
        self.dispComponent = str(settings.value("plotting/dispComponent", "Z"))

        # initialize event data
        if self.recentfiles:
            lastEvent = self.getLastEvent()
            self.data = Data(self, lastEvent)
        else:
            self.data = Data(self)
        self.autodata = Data(self)

        self.dirty = False
        
        # setup UI
        self.setupUi()

        if settings.value("user/FullName", None) is None:
            fulluser = QInputDialog.getText(self, "Enter Name:", "Full name")
            settings.setValue("user/FullName", fulluser)
            settings.setValue("user/Login", getLogin())
        if settings.value("agency_id", None) is None:
            agency = QInputDialog.getText(self,
                                          "Enter authority/institution name:",
                                          "Authority")
            settings.setValue("agency_id", agency)
        self.fname = dict(manual=None, auto=None, loc=None)
        self.fnames = None
        self._stime = None
        structure_setting = settings.value("data/Structure", "PILOT")
        self.dataStructure = DATASTRUCTURE[structure_setting]()
        self.seismicPhase = str(settings.value("phase", "P"))
        if settings.value("data/dataRoot", None) is None:
            dirname = QFileDialog().getExistingDirectory(
                caption='Choose data root ...')
            settings.setValue("data/dataRoot", dirname)
        settings.sync()

        self.filteroptions = {}
        self.pickDlgs = {}
        self.picks = {}
        self.autopicks = {}
        self.loc = False

    def setupUi(self):

        try:
            self.startTime = min(
                [tr.stats.starttime for tr in self.data.wfdata])
        except:
            self.startTime = UTCDateTime()

        pylot_icon = QIcon()
        pylot_icon.addPixmap(QPixmap(':/icons/pylot.png'))

        self.setWindowTitle("PyLoT - do seismic processing the python way")
        self.setWindowIcon(pylot_icon)

        xlab = self.startTime.strftime('seconds since %Y/%m/%d %H:%M:%S (%Z)')

        _widget = QWidget()
        self._main_layout = QVBoxLayout()

        # add event combo box
        self.eventBox = QComboBox()
        self.eventBox.setMaxVisibleItems(30)
        self.eventBox.setEnabled(False)
        self._event_layout = QHBoxLayout()
        self._event_layout.addWidget(QLabel('Event: '))
        self._event_layout.addWidget(self.eventBox)
        self._event_layout.setStretch(1,1) #set stretch of item 1 to 1
        self._main_layout.addLayout(self._event_layout)
        self.eventBox.activated.connect(self.refreshEvents)
        
        # add tabs
        self.tabs = QTabWidget()
        self._main_layout.addWidget(self.tabs)
        self.tabs.currentChanged.connect(self.refreshTabs)

        # create central matplotlib figure canvas widget
        plottitle = "Overview: {0} components ".format(self.getComponent())
        self.dataPlot = WaveformWidget(parent=self, xlabel=xlab, ylabel=None,
                                       title=plottitle)
        self.dataPlot.setCursor(Qt.CrossCursor)
        wf_tab = QtGui.QWidget()
        array_tab = QtGui.QWidget()
        events_tab = QtGui.QWidget()
        self.wf_layout = QtGui.QVBoxLayout()
        self.array_layout = QtGui.QVBoxLayout()
        self.events_layout = QtGui.QVBoxLayout()
        wf_tab.setLayout(self.wf_layout)
        array_tab.setLayout(self.array_layout)
        events_tab.setLayout(self.events_layout)
        self.tabs.addTab(wf_tab, 'Waveform Plot')
        self.tabs.addTab(array_tab, 'Array Map')
        self.tabs.addTab(events_tab, 'Eventlist')
        
        self.wf_layout.addWidget(self.dataPlot)
        self.init_array_tab()
        self.init_event_table()
        self.tabs.setCurrentIndex(0)
        
        quitIcon = self.style().standardIcon(QStyle.SP_MediaStop)
        saveIcon = self.style().standardIcon(QStyle.SP_DriveHDIcon)
        openIcon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
        helpIcon = self.style().standardIcon(QStyle.SP_DialogHelpButton)
        newIcon = self.style().standardIcon(QStyle.SP_FileIcon)
        newFolderIcon = self.style().standardIcon(QStyle.SP_FileDialogNewFolder)
        
        # create resource icons
        locactionicon = QIcon()
        locactionicon.addPixmap(QPixmap(':/icons/locactionicon.png'))
        manupicksicon = QIcon()
        manupicksicon.addPixmap(QPixmap(':/icons/manupicsicon.png'))
        autopicksicon = QIcon()
        autopicksicon.addPixmap(QPixmap(':/icons/autopicsicon.png'))
        loadpiloticon = QIcon()
        loadpiloticon.addPixmap(QPixmap(':/icons/Matlab_PILOT_icon.png'))
        p_icon = QIcon()
        p_icon.addPixmap(QPixmap(':/icons/key_P.png'))
        s_icon = QIcon()
        s_icon.addPixmap(QPixmap(':/icons/key_S.png'))
        print_icon = QIcon()
        print_icon.addPixmap(QPixmap(':/icons/printer.png'))
        filter_icon = QIcon()
        filter_icon.addPixmap(QPixmap(':/icons/filter.png'))
        z_icon = QIcon()
        z_icon.addPixmap(QPixmap(':/icons/key_Z.png'))
        n_icon = QIcon()
        n_icon.addPixmap(QPixmap(':/icons/key_N.png'))
        e_icon = QIcon()
        e_icon.addPixmap(QPixmap(':/icons/key_E.png'))
        auto_icon = QIcon(':/icons/autopick_button.png')
        auto_icon.addPixmap(QPixmap(':/icons/autopick_button.png'))
        locate_icon = QIcon()
        locate_icon.addPixmap(QPixmap(':/icons/locate_button.png'))
        compare_icon = QIcon()
        compare_icon.addPixmap(QPixmap(':/icons/compare_button.png'))
        self.newProjectAction = self.createAction(self, "&New project ...",
                                           self.createNewProject,
                                           QKeySequence.New, newIcon,
                                           "Create a new Project.")
        self.openProjectAction = self.createAction(self, "Load project ...",
                                                   self.loadProject,
                                                   QKeySequence.Open,
                                                   openIcon,
                                                   "Load project file")
        self.saveProjectAction = self.createAction(self, "Save project ...",
                                                   self.saveProject,
                                                   QKeySequence.Save,
                                                   saveIcon,
                                                   "Save project file")
        # newEventAction = self.createAction(self, "&New event ...",
        #                                    self.createNewEvent,
        #                                    QKeySequence.New, newIcon,
        #                                    "Create a new event.")
        self.openmanualpicksaction = self.createAction(self, "Load &picks ...",
                                                  self.load_data,
                                                  QKeySequence.Open,
                                                  manupicksicon,
                                                  "Load manual picks for "
                                                  "the displayed event.")
        self.openmanualpicksaction.setEnabled(False)
        self.openmanualpicksaction.setData(None)

        self.openautopicksaction = self.createAction(self, "Load &automatic picks "
                                                      "...",
                                                self.load_autopicks,
                                                "Ctrl+A",
                                                autopicksicon,
                                                "Load automatic picks "
                                                "for the displayed event.")
        self.openautopicksaction.setEnabled(False)
        self.openautopicksaction.setData(None)

        loadlocationaction = self.createAction(self, "Load &location ...",
                                               self.load_loc, "Ctrl+L",
                                               locactionicon,
                                               "Load location information on "
                                               "the displayed event.")
        self.loadpilotevent = self.createAction(self, "Load PILOT &event ...",
                                                self.load_pilotevent, "Ctrl+E",
                                                loadpiloticon,
                                                "Load PILOT event from information "
                                                "MatLab binary collections (created"
                                                " in former MatLab based version).")
        self.loadpilotevent.setEnabled(False)

        self.saveEventAction = self.createAction(self, "&Save event ...",
                                            self.saveData, QKeySequence.Save,
                                            saveIcon, "Save actual event data.")
        self.saveEventAction.setEnabled(False)

        self.addEventDataAction = self.createAction(self, "Add &events ...",
                                                    self.add_events,
                                                    "Ctrl+W", newFolderIcon,
                                                    "Add event data")
        prefsEventAction = self.createAction(self, "Preferences",
                                             self.PyLoTprefs,
                                             QKeySequence.Preferences,
                                             QIcon(None),
                                             "Edit PyLoT app preferences.")
        quitAction = self.createAction(self, "&Quit",
                                       QCoreApplication.instance().quit,
                                       QKeySequence.Close, quitIcon,
                                       "Close event and quit PyLoT")
        self.filterAction = self.createAction(self, "&Filter ...",
                                              self.filterWaveformData,
                                              "Ctrl+F", filter_icon,
                                              """Toggle un-/filtered waveforms
                                         to be displayed, according to the
                                         desired seismic phase.""", True)
        filterEditAction = self.createAction(self, "&Filter parameter ...",
                                             self.adjustFilterOptions,
                                             "Alt+F", QIcon(None),
                                             """Adjust filter parameters.""")
        self.selectPAction = self.createAction(self, "&P", self.alterPhase,
                                               "Alt+P",
                                               p_icon,
                                               "Toggle P phase.", True)
        self.selectSAction = self.createAction(self, "&S", self.alterPhase,
                                               "Alt+S",
                                               s_icon,
                                               "Toggle S phase", True)
        self.compare_action = self.createAction(self, "&Compare picks...",
                                                self.comparePicks, "Alt+C",
                                                compare_icon, "Comparison of "
                                                              "manual and "
                                                              "automatic pick "
                                                              "data.", False)
        self.compare_action.setEnabled(False)

        printAction = self.createAction(self, "&Print event ...",
                                        self.show_event_information, QKeySequence.Print,
                                        print_icon,
                                        "Print waveform section.")
        helpAction = self.createAction(self, "&Help ...", self.helpHelp,
                                       QKeySequence.HelpContents, helpIcon,
                                       """Show either the documentation
                                       homepage (internet connection available),
                                       or shipped documentation files.""")
        self.fileMenu = self.menuBar().addMenu('&File')
        self.fileMenuActions = (self.newProjectAction, self.addEventDataAction,
                                self.openProjectAction, self.saveProjectAction,
                                self.openmanualpicksaction, self.saveEventAction, None,
                                prefsEventAction, quitAction)
        self.fileMenu.aboutToShow.connect(self.updateFileMenu)
        self.updateFileMenu()

        self.editMenu = self.menuBar().addMenu('&Edit')
        editActions = (self.filterAction, filterEditAction, None,
                       self.selectPAction, self.selectSAction, None,
                       printAction)
        self.addActions(self.editMenu, editActions)

        self.helpMenu = self.menuBar().addMenu('&Help')
        helpActions = (helpAction,)
        self.addActions(self.helpMenu, helpActions)

        fileToolBar = self.addToolBar("FileTools")
        fileToolActions = (self.newProjectAction, self.addEventDataAction,
                           self.openProjectAction, self.saveProjectAction,
                           self.openmanualpicksaction,
                           self.openautopicksaction, loadlocationaction,
                           self.loadpilotevent, self.saveEventAction)
        fileToolBar.setObjectName("FileTools")
        self.addActions(fileToolBar, fileToolActions)

        # phaseToolBar = self.addToolBar("PhaseTools")
        # phaseToolActions = (self.selectPAction, self.selectSAction)
        # phaseToolBar.setObjectName("PhaseTools")
        # self.addActions(phaseToolBar, phaseToolActions)

        # create button group for component selection

        componentGroup = QActionGroup(self)
        componentGroup.setExclusive(True)

        self.z_action = self.createAction(parent=componentGroup, text='Z',
                                     slot=self.plotZ, shortcut='Alt+Z',
                                     icon=z_icon, tip='Display the vertical (Z)'
                                                      ' component.',
                                     checkable=True)
        self.z_action.setChecked(True)
        self.z_action.setEnabled(False)

        self.n_action = self.createAction(parent=componentGroup, text='N',
                                     slot=self.plotN, shortcut='Alt+N',
                                     icon=n_icon,
                                     tip='Display the north-south (N) '
                                         'component.', checkable=True)
        self.n_action.setEnabled(False)

        self.e_action = self.createAction(parent=componentGroup, text='E',
                                     slot=self.plotE, shortcut='Alt+E',
                                     icon=e_icon,
                                     tip='Display the east-west (E) component.',
                                     checkable=True)
        self.e_action.setEnabled(False)

        componentToolBar = self.addToolBar("ComponentSelection")
        componentActions = (self.z_action, self.n_action, self.e_action)
        componentToolBar.setObjectName("PhaseTools")
        self.addActions(componentToolBar, componentActions)

        self.auto_pick = self.createAction(parent=self, text='autoPick',
                                      slot=self.autoPick, shortcut='Alt+Ctrl+A',
                                      icon=auto_icon, tip='Automatically pick'
                                                          ' the displayed waveforms.')
        self.auto_pick.setEnabled(False)

        autoPickToolBar = self.addToolBar("autoPyLoT")
        autoPickActions = (self.auto_pick, self.compare_action)
        self.addActions(autoPickToolBar, autoPickActions)

        # pickToolBar = self.addToolBar("PickTools")
        # pickToolActions = (selectStation, )
        # pickToolBar.setObjectName("PickTools")
        # self.addActions(pickToolBar, pickToolActions)
        self.locateEvent = self.createAction(parent=self, text='locate the event',
                                        slot=self.locate_event,
                                        shortcut='Alt+Ctrl+L',
                                        icon=locate_icon,
                                        tip='Locate the event using '
                                            'the displayed manual arrivals.')
        self.locateEvent.setEnabled(False)

        locationToolBar = self.addToolBar("LocationTools")
        locationToolActions = (self.locateEvent,)
        locationToolBar.setObjectName("LocationTools")
        self.addActions(locationToolBar, locationToolActions)

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 500)

        _widget.setLayout(self._main_layout)
        _widget.showFullScreen()

        self.setCentralWidget(_widget)


    @property
    def metadata(self):
        return self._metadata


    @metadata.setter
    def metadata(self, value):
        self._metadata = value


    def updateFileMenu(self):

        self.fileMenu.clear()
        for action in self.fileMenuActions[:-1]:
            if action is None:
                self.fileMenu.addSeparator()
            else:
                self.fileMenu.addAction(action)
        try:
            current = self.data.getID()
        except AttributeError:
            current = None
        recentEvents = []
        for eventID in self.recentfiles:
            fname = fnConstructor(eventID)
            if eventID != current and QFile.exists(fname):
                recentEvents.append(eventID)
                recentEvents.reverse()
                self.recentfiles = recentEvents[0:5]
                settings = QSettings()
                settings.setValue()
        if recentEvents:
            for i, eventID in enumerate(recentEvents):
                fname = fnConstructor(eventID)
                action = QAction(self.windowIcon(),
                                 "&{0} {1}".format(i + 1,
                                                   QFileInfo(fname).fileName()),
                                 self)
                action.setData(fname)
                self.connect(action, Signal("triggered()"),
                             self.load_data)
                self.fileMenu.addAction(action)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.fileMenuActions[-1])

    @property
    def inputs(self):
        return self._inputs

    def getRoot(self):
        settings = QSettings()
        return settings.value("data/dataRoot")

    def load_autopicks(self, fname=None):
        self.load_data(fname, type='auto')

    def load_loc(self, fname=None):
        type = getDataType(self)
        self.load_data(fname, type=type, loc=True)

    def load_pilotevent(self):
        filt = "PILOT location files (*LOC*.mat)"
        caption = "Select PILOT location file"
        fn_loc = QFileDialog().getOpenFileName(self, caption=caption,
                                              filter=filt, dir=self.getRoot())
        fn_loc = fn_loc[0]
        loc_dir = os.path.split(fn_loc)[0]
        filt = "PILOT phases files (*PHASES*.mat)"
        caption = "Select PILOT phases file"
        fn_phases = QFileDialog().getOpenFileName(self, caption=caption,
                                              filter=filt, dir=loc_dir)
        fn_phases = fn_phases[0]

        type = getDataType(self)

        fname_dict = dict(phasfn=fn_phases, locfn=fn_loc)
        self.load_data(fname_dict, type=type)

    def load_data(self, fname=None, type='manual', loc=False):
        if not self.okToContinue():
            return
        if fname is None:
            action = self.sender()
            if isinstance(action, QAction):
                fname = self.filename_from_action(action)
        self.set_fname(fname, type)
        data = dict(auto=self.autodata, manual=self.data)
        data[type] += Data(self, evtdata=fname)
        if not loc:
            self.updatePicks(type=type)
        self.drawPicks(picktype=type)
        self.draw()

    def getCurrentEvent(self):
        for event in self.project.eventlist:
            if event.path == self.eventBox.currentText():
                return event
        
    def getLastEvent(self):
        return self.recentfiles[0]

    def add_recentfile(self, event):
        self.recentfiles.insert(0, event)

    def getWFFnames(self):
        try:
            evt = self.get_data().get_evt_data()
            if evt.picks:
                for pick in evt.picks:
                    try:
                        if pick.waveform_id is not None:
                            fname = pick.waveform_id.getSEEDstring()
                            if fname not in self.fnames:
                                self.fnames.append(fname)
                    except:
                        continue
            else:
                if self.dataStructure:
                    searchPath = self.dataStructure.expandDataPath()
                    fnames = QFileDialog.getOpenFileNames(self,
                                                          "Select waveform "
                                                          "files:",
                                                          dir=searchPath)
                    self.fnames = fnames[0]

                else:
                    raise DatastructureError('not specified')
            if not self.fnames:
                return None
            return self.fnames
        except DatastructureError as e:
            print(e)
            props = PropertiesDlg(self)
            if props.exec_() == QDialog.Accepted:
                return self.getWFFnames()
            else:
                return

    def getWFFnames_from_eventlist(self):
        if self.dataStructure:
            searchPath = self.dataStructure.expandDataPath()
            directory = self.eventBox.currentText()
            self.fnames = [os.path.join(directory, f) for f in os.listdir(directory)]
        else:
            raise DatastructureError('not specified')
        if not self.fnames:
            return None
        return self.fnames

    def add_events(self):
        if not self.project:
            self.project = Project()
        ed = getExistingDirectories(self, 'Select event directories...')
        if ed.exec_():
            eventlist = ed.selectedFiles()
            # select only folders that start with 'e', containin two dots and have length 12
            eventlist = [item for item in eventlist if item.split('/')[-1].startswith('e')
                         and len(item.split('/')[-1].split('.')) == 3
                         and len(item.split('/')[-1]) == 12]
        else:
            return
        if not self.project:
            print('No project found.')
            return
        self.project.add_eventlist(eventlist)
        self.init_events()

    def init_events(self, new=False):
        nitems = self.eventBox.count()
        self.eventBox.clear()
        if len(self.project.eventlist) == 0:
            print('No events to init.')
            self.clearWaveformDataPlot()
            return
        self.eventBox.setEnabled(True)
        self.fill_eventbox(self.project.eventlist)
        if new:
            self.eventBox.setCurrentIndex(0)
        else:
            self.eventBox.setCurrentIndex(nitems)
        self.refreshEvents()
        tabindex = self.tabs.currentIndex()

    def fill_eventbox(self, eventlist):
        model = self.eventBox.model()
        for event in self.project.eventlist:
            event = event.path
            item = QtGui.QStandardItem(str(event))
            # if ref: set different color e.g.
            #item.setBackground(QtGui.QColor('teal'))
            item.setForeground(QtGui.QColor('black'))
            font = item.font()
            font.setPointSize(10)
            item.setFont(font)
            model.appendRow(item)
        
    def filename_from_action(self, action):
        if action.data() is None:
            filt = "Supported file formats" \
                   " (*.mat *.qml *.xml *.kor *.evt)"
            caption = "Open an event file"
            fname = QFileDialog().getOpenFileName(self, caption=caption,
                                                  filter=filt,
                                                  dir=self.getRoot())
            fname = fname[0]
        else:
            fname = str(action.data().toString())
        return fname

    def get_fnames(self, type='manual'):
        return self.fname[type]

    def set_fname(self, fname, type):
        if self.get_fnames(type) is not None:
            self.add_recentfile(self.get_fnames(type))
        self.fname[type] = fname

    def getEventFileName(self, type='manual'):
        if self.get_fnames(type) is None:
            self.set_fname(self.get_data().getEventFileName(), type)
        return self.get_fnames(type)

    def saveData(self):

        def getSavePath(e):
            print('warning: {0}'.format(e))
            directory = os.path.realpath(self.getRoot())
            file_filter = "QuakeML file (*.xml);;VELEST observation file " \
                          "format (*.cnv);;NonLinLoc observation file (*.obs)"
            title = 'Save event data ...'
            fname, selected_filter = QFileDialog.getSaveFileName(self,
                                                                 title,
                                                                 directory,
                                                                 file_filter)

            fbasename, exform = os.path.splitext(fname)

            if not exform and selected_filter:
                exform = selected_filter.split('*')[1][:-1]

            return fbasename, exform

        settings = QSettings()
        fbasename = self.getEventFileName()
        exform = settings.value('data/exportFormat', 'QUAKEML')
        try:
            self.get_data().applyEVTData(self.getPicks())
        except OverwriteError:
            msgBox = QMessageBox()
            msgBox.setText("Picks have been modified!")
            msgBox.setInformativeText(
                "Do you want to save the changes and overwrite the picks?")
            msgBox.setDetailedText(self.get_data().getPicksStr())
            msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
            msgBox.setDefaultButton(QMessageBox.Save)
            ret = msgBox.exec_()
            if ret == QMessageBox.Save:
                self.get_data().resetPicks()
                return self.saveData()
            elif ret == QMessageBox.Cancel:
                return False
        try:
            self.get_data().exportEvent(fbasename, exform)
        except FormatError as e:
            fbasename, exform = getSavePath(e)
        except AttributeError as e:
            fbasename, exform = getSavePath(e)

        # catch all possible cases before going on
        if not fbasename:
            return False
        # warn overwriting
        elif os.path.exists(fbasename + exform):
            ans = QMessageBox.question(self, self.tr("Overwrite file..."),
                               self.tr("File already exists: {0}\n".format(fbasename + exform) + \
                                  "Overwrite file anyway?"),
                               QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            # only negative answers have to be caught
            if ans == QMessageBox.No:
                self.saveData()
            elif ans == QMessageBox.Cancel:
                return False

        # export to given path
        self.get_data().exportEvent(fbasename, exform)
        # all files save (ui clean)
        self.setDirty(False)
        self.update_status('Event saved as %s' % (fbasename + exform))
        return True

    def getinfile(self):
        return self.infile

    def getComponent(self):
        return self.dispComponent

    def setComponent(self, component):
        self.dispComponent = component

    def get_data(self, type='manual'):
        if type == 'auto':
            return self.autodata
        return self.data

    def getPicks(self, type='manual'):
        rdict = dict(auto=self.autopicks, manual=self.picks)
        return rdict[type]

    def getPicksOnStation(self, station, type='manual'):
        try:
            return self.getPicks(type)[station]
        except KeyError:
            return None

    def comparePicks(self):
        if self.check4Comparison():
            co = Comparison(auto=self.getPicks('auto'), manu=self.getPicks())
            compare_dlg = ComparisonDialog(co, self)
            compare_dlg.exec_()

    def getPlotWidget(self):
        return self.dataPlot

    @staticmethod
    def getWFID(gui_event):

        ycoord = gui_event.ydata

        try:
            statID = int(round(ycoord))
        except TypeError as e:
            if 'a float is required' in e.message:
                return None
            else:
                raise e

        return statID

    def getStationID(self, station):
        for wfID in self.getPlotWidget().getPlotDict().keys():
            actual_station = self.getPlotWidget().getPlotDict()[wfID][0]
            if station == actual_station:
                return wfID
        return None

    def getStime(self):
        return self._stime

    def addActions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def okToContinue(self):
        if self.dirty:
            return self.saveData()
        return True

    def refreshEvents(self):
        self._eventChanged = True
        self.refreshTabs()
        
    def refreshTabs(self):
        if self.tabs.currentIndex() == 0:
            if hasattr(self.project, 'eventlist'):
                if len(self.project.eventlist) > 0:
                    if self._eventChanged:
                        self.newWFplot()
        if self.tabs.currentIndex() == 1:
            self.refresh_array_map()
        if self.tabs.currentIndex() == 2:
            self.init_event_table()

    def newWFplot(self):
        self.loadWaveformDataThread()
        self._eventChanged = False
    
    def loadWaveformDataThread(self):
        wfd_thread = Thread(self, self.loadWaveformData, progressText='Reading data input...')
        wfd_thread.finished.connect(self.plotWaveformDataThread)
        wfd_thread.start()
        
    def loadWaveformData(self):
        # if self.fnames and self.okToContinue():
        #     self.setDirty(True)
        #     ans = self.data.setWFData(self.fnames)
        # elif self.fnames is None and self.okToContinue():
        #     ans = self.data.setWFData(self.getWFFnames())
        # else:
        #     ans = False
        self.data.setWFData(self.getWFFnames_from_eventlist())
        self._stime = full_range(self.get_data().getWFData())[0]

    def connectWFplotEvents(self):
        if not self.poS_id:
            self.poS_id = self.dataPlot.mpl_connect('button_press_event',
                                                    self.pickOnStation)
        if not self.ae_id:
            self.ae_id = self.dataPlot.mpl_connect('axes_enter_event',
                                                   lambda event: self.tutor_user())

    def disconnectWFplotEvents(self):
        if self.poS_id:
            self.dataPlot.mpl_disconnect(self.poS_id)
        if self.ae_id:
            self.dataPlot.mpl_disconnect(self.ae_id)
        self.poS_id = None
        self.ae_id = None

    def finishWaveformDataPlot(self):
        self.connectWFplotEvents()
        self.auto_pick.setEnabled(True)
        self.z_action.setEnabled(True)
        self.e_action.setEnabled(True)
        self.n_action.setEnabled(True)
        self.openmanualpicksaction.setEnabled(True)
        self.openautopicksaction.setEnabled(True)
        self.loadpilotevent.setEnabled(True)
        self.saveEventAction.setEnabled(True)
        event = self.getCurrentEvent()
        if event.picks:
            self.picks = event.picks
            self.drawPicks(picktype='manual')
        if event.autopicks:
            self.autopicks = event.autopicks
            self.drawPicks(picktype='auto')
        self.draw()

    def clearWaveformDataPlot(self):
        self.disconnectWFplotEvents()
        self.dataPlot.getAxes().cla()
        self.auto_pick.setEnabled(False)
        self.z_action.setEnabled(False)
        self.e_action.setEnabled(False)
        self.n_action.setEnabled(False)
        self.openmanualpicksaction.setEnabled(False)
        self.openautopicksaction.setEnabled(False)
        self.loadpilotevent.setEnabled(False)
        self.saveEventAction.setEnabled(False)
        self.draw()
        
    def plotWaveformDataThread(self):
        wfp_thread = Thread(self, self.plotWaveformData, progressText='Plotting waveform data...')
        wfp_thread.finished.connect(self.finishWaveformDataPlot)
        wfp_thread.start()
        
    def plotWaveformData(self):
        zne_text = {'Z': 'vertical', 'N': 'north-south', 'E': 'east-west'}
        comp = self.getComponent()
        title = 'section: {0} components'.format(zne_text[comp])
        alter_comp = COMPNAME_MAP[comp]
        wfst = self.get_data().getWFData().select(component=comp)
        wfst += self.get_data().getWFData().select(component=alter_comp)
        self.getPlotWidget().plotWFData(wfdata=wfst, title=title, mapping=False)
        plotDict = self.getPlotWidget().getPlotDict()
        pos = plotDict.keys()
        labels = [plotDict[n][0] for n in pos]
        self.getPlotWidget().setYTickLabels(pos, labels)

    def plotZ(self):
        self.setComponent('Z')
        self.plotWaveformDataThread()
        self.drawPicks()
        self.draw()

    def plotN(self):
        self.setComponent('N')
        self.plotWaveformDataThread()
        self.drawPicks()
        self.draw()

    def plotE(self):
        self.setComponent('E')
        self.plotWaveformDataThread()
        self.drawPicks()
        self.draw()

    def pushFilterWF(self, param_args):
        self.get_data().filterWFData(param_args)

    def filterWaveformData(self):
        if self.get_data():
            if self.getFilterOptions() and self.filterAction.isChecked():
                kwargs = self.getFilterOptions().parseFilterOptions()
                self.pushFilterWF(kwargs)
            elif self.filterAction.isChecked():
                self.adjustFilterOptions()
            else:
                self.get_data().resetWFData()
        self.plotWaveformData()
        self.drawPicks()
        self.draw()

    def adjustFilterOptions(self):
        fstring = "Filter Options ({0})".format(self.getSeismicPhase())
        filterDlg = FilterOptionsDialog(titleString=fstring,
                                        parent=self)
        if filterDlg.exec_():
            filteroptions = filterDlg.getFilterOptions()
            self.setFilterOptions(filteroptions)
            if self.filterAction.isChecked():
                kwargs = self.getFilterOptions().parseFilterOptions()
                self.pushFilterWF(kwargs)
                self.plotWaveformData()

    def getFilterOptions(self):
        try:
            return self.filteroptions[self.getSeismicPhase()]
        except AttributeError as e:
            print(e)
            return FilterOptions(None, None, None)

    def getFilters(self):
        return self.filteroptions

    def setFilterOptions(self, filterOptions, seismicPhase=None):
        if seismicPhase is None:
            self.getFilters()[self.getSeismicPhase()] = filterOptions
        else:
            self.getFilters()[seismicPhase] = filterOptions

    def updateFilterOptions(self):
        try:
            settings = QSettings()
            if settings.value("filterdefaults",
                              None) is None and not self.getFilters():
                for key, value in FILTERDEFAULTS.items():
                    self.setFilterOptions(FilterOptions(**value), key)
            elif settings.value("filterdefaults", None) is not None:
                for key, value in settings.value("filterdefaults"):
                    self.setFilterOptions(FilterOptions(**value), key)
        except Exception as e:
            self.update_status('Error ...')
            emsg = QErrorMessage(self)
            emsg.showMessage('Error: {0}'.format(e))
        else:
            self.update_status('Filter loaded ... '
                              '[{0}: {1} Hz]'.format(
                self.getFilterOptions().getFilterType(),
                self.getFilterOptions().getFreq()))
        if self.filterAction.isChecked():
            self.filterWaveformData()

    def getSeismicPhase(self):
        return self.seismicPhase

    def getStationName(self, wfID):
        return self.getPlotWidget().getPlotDict()[wfID][0]

    def alterPhase(self):
        pass

    def setSeismicPhase(self, phase):
        self.seismicPhase = self.seismicPhaseButtonGroup.getValue()
        self.update_status('Seismic phase changed to '
                          '{0}'.format(self.getSeismicPhase()))

    def pickOnStation(self, gui_event):

        wfID = self.getWFID(gui_event)

        if wfID is None: return

        station = self.getStationName(wfID)
        self.update_status('picking on station {0}'.format(station))
        data = self.get_data().getWFData()
        pickDlg = PickDlg(self, infile=self.getinfile(), 
                          data=data.select(station=station),
                          station=station,
                          picks=self.getPicksOnStation(station, 'manual'),
                          autopicks=self.getPicksOnStation(station, 'auto'))
        if pickDlg.exec_():
            self.setDirty(True)
            self.update_status('picks accepted ({0})'.format(station))
            replot = self.addPicks(station, pickDlg.getPicks())
            if replot:
                self.plotWaveformData()
                self.drawPicks()
                self.draw()
            else:
                self.drawPicks(station)
                self.draw()
        else:
            self.update_status('picks discarded ({0})'.format(station))
        if not self.get_loc_flag() and self.check4Loc():
            self.locateEvent.setEnabled(True)
            self.set_loc_flag(True)
        elif self.get_loc_flag() and not self.check4Loc():
            self.set_loc_flag(False)

    def addListItem(self, text):
        self.listWidget.addItem(text)
        self.listWidget.scrollToBottom()

    def autoPick(self):
        self.autosave = QFileDialog().getExistingDirectory(caption='Select autoPyLoT output') 
        if not os.path.exists(self.autosave):
                QMessageBox.warning(self, "PyLoT Warning",
                           "No autoPyLoT output declared!")
                return
        self.listWidget = QListWidget()
        self.setDirty(True)
        self.logDockWidget = QDockWidget("AutoPickLog", self)
        self.logDockWidget.setObjectName("LogDockWidget")
        self.logDockWidget.setAllowedAreas(
            Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.logDockWidget.setWidget(self.listWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.logDockWidget)
        self.addListItem('Loading default values from PyLoT-input file %s' 
                                                             % self.infile)
        autopick_parameter = self._inputs
        self.addListItem(str(autopick_parameter))

        self.thread = AutoPickThread(parent=self,
                                     func=autoPyLoT,
                                     infile = self.infile, 
                                     fnames=self.fnames,
                                     savepath=self.autosave)

        self.thread.message.connect(self.addListItem)
        self.thread.start()
        self.thread.finished.connect(self.finalizeAutoPick)

    def finalizeAutoPick(self):
        self.drawPicks(picktype='auto')
        self.draw()
        self.thread.quit()

    def addPicks(self, station, picks, type='manual'):
        stat_picks = self.getPicksOnStation(station, type)
        rval = False
        if not stat_picks:
            stat_picks = picks
        else:
            msgBox = QMessageBox()
            msgBox.setText("The picks for station {0} have been "
                           "changed.".format(station))
            msgBox.setDetailedText("Old picks:\n"
                                   "{old_picks}\n\n"
                                   "New picks:\n"
                                   "{new_picks}".format(old_picks=stat_picks,
                                                        new_picks=picks))
            msgBox.setInformativeText("Do you want to save your changes?")
            msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
            msgBox.setDefaultButton(QMessageBox.Save)
            ret = msgBox.exec_()
            if ret == QMessageBox.Save:
                stat_picks = picks
                rval = True
            elif ret == QMessageBox.Cancel:
                pass
            else:
                raise Exception('FATAL: Should never occur!')
        self.getPicks(type=type)[station] = stat_picks
        return rval

    def updatePicks(self, type='manual'):
        picks = picksdict_from_picks(evt=self.get_data(type).get_evt_data())
        if type == 'manual':
            self.getCurrentEvent().addPicks(picks)
            self.picks.update(picks)
        elif type == 'auto':
            self.getCurrentEvent().addAutopicks(picks)            
            self.autopicks.update(picks)
        self.check4Comparison()

    def drawPicks(self, station=None, picktype='manual'):
        # if picks to draw not specified, draw all picks available
        if not station:
            for station in self.getPicks(type=picktype):
                self.drawPicks(station, picktype=picktype)
            return
        # plotting picks
        plotID = self.getStationID(station)
        if not plotID:
            return
        ax = self.getPlotWidget().axes
        ylims = np.array([-.5, +.5]) + plotID
        phase_col = {
            'P': ('c', 'c--', 'b-', 'bv', 'b^', 'b'),
            'S': ('m', 'm--', 'r-', 'rv', 'r^', 'r')
        }

        stat_picks = self.getPicks(type=picktype)[station]

        stime = self.getStime()

        for phase in stat_picks:
            picks = stat_picks[phase]
            if type(stat_picks[phase]) is not dict:
                return
            colors = phase_col[phase[0].upper()]

            mpp = picks['mpp'] - stime
            epp = picks['epp'] - stime
            lpp = picks['lpp'] - stime
            spe = picks['spe']
            if not spe:
                spe = symmetrize_error(mpp - epp, lpp - mpp)

            if picktype == 'manual':
                ax.fill_between([epp, lpp], ylims[0], ylims[1],
                                alpha=.25, color=colors[0])
                ax.plot([mpp - spe, mpp - spe], ylims, colors[1],
                        [mpp, mpp], ylims, colors[2],
                        [mpp + spe, mpp + spe], ylims, colors[1])
            elif picktype == 'auto':
                ax.plot(mpp, ylims[1], colors[3],
                        mpp, ylims[0], colors[4])
                ax.vlines(mpp, ylims[0], ylims[1], colors[5], linestyles='dotted')                
            else:
                raise TypeError('Unknown picktype {0}'.format(picktype))

    def locate_event(self):
        """
        locate event using the manually picked phases
        :return:
        """
        if not self.okToContinue():
            return
        settings = QSettings()
        # get location tool hook
        loctool = settings.value("loc/tool", "nll")
        lt = locateTool[loctool]
        # get working directory
        locroot = settings.value("{0}/rootPath".format(loctool), None)
        if locroot is None:
            self.PyLoTprefs()
            self.locate_event()

        infile = settings.value("{0}/inputFile".format(loctool), None)

        if not infile:
            caption = 'Select {0} input file'.format(loctool)
            filt = "Supported file formats" \
                   " (*.in *.ini *.conf *.cfg)"
            ans = QFileDialog().getOpenFileName(self, caption=caption,
                                                filter=filt, dir=locroot)
            if ans[0]:
                infile = ans[0]
            else:
                QMessageBox.information(self,
                                        self.tr('No infile selected'),
                                        self.tr('Inputfile necessary for localization.'))
                return
            settings.setValue("{0}/inputFile".format(loctool), infile)
            settings.sync()
        if loctool == 'nll':
            ttt = settings.value("{0}/travelTimeTables", None)
            ok = False
            if ttt is None:
                while not ok:
                    text, ok = QInputDialog.getText(self, 'Pattern for travel time tables',
                                                    'Base name of travel time tables',
                                                    echo=QLineEdit.Normal,
                                                    text="ttime")
                ttt = text

        outfile = settings.value("{0}/outputFile".format(loctool),
                                 os.path.split(os.tempnam())[-1])
        phasefile = os.path.split(os.tempnam())[-1]
        phasepath = os.path.join(locroot, 'obs', phasefile)
        locpath = os.path.join(locroot, 'loc', outfile)
        lt.export(self.getPicks(), phasepath)
        lt.modify_inputs(infile, locroot, outfile, phasefile, ttt)
        try:
            lt.locate(infile)
        except RuntimeError as e:
            print(e.message)
        finally:
            os.remove(phasepath)

        self.get_data().applyEVTData(lt.read_location(locpath), type='event')
        self.get_data().applyEVTData(self.calc_magnitude(), type='event')

    def init_array_tab(self):
        self.metadata_widget = QWidget(self)
        grid_layout = QGridLayout()
        grid_layout.setColumnStretch(0, 1)
        grid_layout.setColumnStretch(2, 1)
        grid_layout.setRowStretch(0, 1)
        grid_layout.setRowStretch(3, 1)

        label = QLabel('No inventory set...')
        new_inv_button = QPushButton('Set &inventory file')
        new_inv_button.clicked.connect(self.get_metadata)
        
        grid_layout.addWidget(label, 1, 1)
        grid_layout.addWidget(new_inv_button, 2, 1)

        self.metadata_widget.setLayout(grid_layout)
        self.array_layout.addWidget(self.metadata_widget)
    
    def init_array_map(self, index=1):
        if not self.array_map:
            self.get_metadata()
            if not self.metadata:
                return
        self.array_map = map_projection(self)
        self.array_layout.removeWidget(self.metadata_widget)
        self.array_layout.addWidget(self.array_map)
        self.tabs.setCurrentIndex(index)
        
    def refresh_array_map(self):
        if not self.array_map:
            return
        # refresh with new picks here!!!
        self.array_map.show()

    def init_event_table(self, index=2):
        def set_enabled(item, enabled=True, checkable=False):
            if enabled and not checkable:
                item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
            elif enabled and checkable:
                item_ref.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsSelectable)                
            else:
                item.setFlags(QtCore.Qt.ItemIsSelectable)

        if self.project:
            eventlist = self.project.eventlist
        else:
            eventlist = []

        def cell_changed(row=None, column=None):
            table = self.project._table
            event = self.project.getEventFromPath(table[row][0].text())
            if column == 1 or column == 2:
                #toggle checked states (exclusive)
                item_ref = table[row][1]
                item_test = table[row][2]
                if column == 1 and item_ref.checkState():
                    item_test.setCheckState(QtCore.Qt.Unchecked)
                    event.setRefEvent(True)
                elif column == 1 and not item_ref.checkState():
                    event.setRefEvent(False)                    
                elif column == 2 and item_test.checkState():
                    item_ref.setCheckState(QtCore.Qt.Unchecked)
                    event.setTestEvent(True)
                elif column == 2 and not item_test.checkState():
                    event.setTestEvent(False)                    
            elif column == 3:
                #update event notes
                notes = table[row][3].text()
                event.addNotes(notes)

        if hasattr(self, 'qtl'):
            self.events_layout.removeWidget(self.qtl)
        self.qtl = QtGui.QTableWidget()
        self.qtl.setColumnCount(4)
        self.qtl.setRowCount(len(eventlist))
        self.qtl.setHorizontalHeaderLabels(['Event', 'Reference', 'Test Set', 'Notes'])

        self.project._table = []
        for index, event in enumerate(eventlist):
            item_path = QtGui.QTableWidgetItem()
            item_ref = QtGui.QTableWidgetItem()
            item_test = QtGui.QTableWidgetItem()
            item_notes = QtGui.QTableWidgetItem()            

            item_path.setText(event.path)
            item_notes.setText(event.notes)            
            if event.picks:
                set_enabled(item_path, True, False)
                set_enabled(item_ref, True, True)
                set_enabled(item_test, True, True)
            else:
                set_enabled(item_path, False, False)
                set_enabled(item_ref, False, True)
                set_enabled(item_test, False, True)

            if event.isRefEvent():
                item_ref.setCheckState(QtCore.Qt.Checked)
            else:
                item_ref.setCheckState(QtCore.Qt.Unchecked)
            if event.isTestEvent():
                item_test.setCheckState(QtCore.Qt.Checked)
            else:
                item_test.setCheckState(QtCore.Qt.Unchecked)
                
            column=[item_path, item_ref, item_test, item_notes]
            self.project._table.append(column)

        for r_index, row in enumerate(self.project._table):
            for c_index, item in enumerate(row):
                self.qtl.setItem(r_index, c_index, item)

        header = self.qtl.horizontalHeader()
        header.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)
        self.qtl.cellChanged[int, int].connect(cell_changed)

        self.events_layout.addWidget(self.qtl)
        self.tabs.setCurrentIndex(index)
        
    def read_metadata_thread(self, fninv):
        self.rm_thread = Thread(self, read_metadata, arg=fninv, progressText='Reading metadata...')
        self.rm_thread.finished.connect(self.set_metadata)
        self.rm_thread.start()

    def set_metadata(self):
        self.metadata = self.rm_thread.data
        self.project.metadata = self.rm_thread.data
        self.init_array_map()
        
    def get_metadata(self):
        def set_inv(settings):
            fninv, _ = QFileDialog.getOpenFileName(self, self.tr(
                "Select inventory..."), self.tr("Select file"))
            if not fninv:
                return False
            ans = QMessageBox.question(self, self.tr("Make default..."),
                                       self.tr(
                                           "New inventory filename set.\n" + \
                                           "Do you want to make it the default value?"),
                                       QMessageBox.Yes | QMessageBox.No,
                                       QMessageBox.No)
            if ans == QMessageBox.Yes:
                settings.setValue("inventoryFile", fninv)
                settings.sync()
            self.read_metadata_thread(fninv)
            return True

        if hasattr(self.project, 'metadata'):
            self.metadata = self.project.metadata
            return True
        
        settings = QSettings()
        fninv = settings.value("inventoryFile", None)

        if fninv is None and not self.metadata:
            if not set_inv(settings):
                return None
        elif fninv is not None and not self.metadata:
            ans = QMessageBox.question(self, self.tr("Use default..."),
                                       self.tr(
                                           "Do you want to use the default value?"),
                                       QMessageBox.Yes | QMessageBox.No,
                                       QMessageBox.Yes)
            if ans == QMessageBox.No:
                if not set_inv(settings):
                    return None
            else:
                self.read_metadata_thread(fninv)
        
    def calc_magnitude(self, type='ML'):
        self.get_metadata()
        if not self.metadata:
            return None
            
        wf_copy = self.get_data().getWFData().copy()
        corr_wf = restitute_data(wf_copy, *self.metadata)
        # if not rest_flag:
        #     raise ProcessingError('Restitution of waveform data failed!')
        if type == 'ML':
            local_mag = RichterMagnitude(corr_wf, self.get_data().get_evt_data(), self.inputs.get('sstop'), verbosity = True)
            return local_mag.updated_event()
        elif type == 'Mw':
            moment_mag = MomentMagnitude(corr_wf, self.get_data().get_evt_data(), self.inputs.get('vp'), self.inputs.get('Qp'), self.inputs.get('rho'), verbosity = True)
            return moment_mag.updated_event()
        else:
            return None

    def check4Loc(self):
        return self.picksNum() >= 4

    def check4Comparison(self):
        mpicks = self.getPicks()
        apicks = self.getPicks('auto')
        for station, phases in mpicks.items():
            try:
                aphases = apicks[station]
                for phase in phases.keys():
                    if phase in aphases.keys():
                        return True
            except KeyError:
                continue
        return False

    def picksNum(self, type='manual'):
        num = 0
        for phases in self.getPicks(type).values():
            num += len(phases)
        return num

    def get_loc_flag(self):
        return self.loc

    def set_loc_flag(self, value):
        self.loc = value

    def check_loc_plt(self):
        evt = self.get_data().get_evt_data()
        if evt.origins and evt.magnitudes:
            return True
        return False

    def update_status(self, message, duration=5000):
        self.statusBar().showMessage(message, duration)
        if self.get_data() is not None:
            if not self.get_data().isNew():
                self.setWindowTitle(
                    "PyLoT - processing event %s[*]" % self.get_data().getID())
            elif self.get_data().isNew():
                self.setWindowTitle("PyLoT - New project [*]")
            else:
                self.setWindowTitle(
                    "PyLoT - seismic processing the python way[*]")
        self.setWindowModified(self.dirty)

    def tutor_user(self):
        self.update_status('select trace to pick on station ...', 10000)

    def show_event_information(self):
        pass

    def createNewEvent(self):
        if self.okToContinue():
            new = NewEventDlg()
            if new.exec_() != QDialog.Rejected:
                evtpar = new.getValues()
                cinfo = create_creation_info(agency_id=self.agency)
                event = create_event(evtpar['origintime'], cinfo)
                self.data = Data(self, evtdata=event)
                self.setDirty(True)

    def createNewProject(self, exists=False):
        if self.okToContinue():
            dlg = QFileDialog()
            fnm = dlg.getSaveFileName(self, 'Create a new project file...', filter='Pylot project (*.plp)')
            filename = fnm[0]
            if not filename.split('.')[-1] == 'plp':
                filename = fnm[0] + '.plp'
            if not exists:
                self.project = Project()
                self.init_events(new=True)                
            self.project.save(filename)
            
    def loadProject(self):
        if self.project:
            if self.project.dirty:
                qmb = QMessageBox(icon=QMessageBox.Question, text='Save changes in current project?')
                qmb.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
                qmb.setDefaultButton(QMessageBox.Yes)
                if qmb.exec_() == 16384:
                    self.saveProject()
                elif qmb.exec_() == 65536:
                    pass
                elif qmb.exec_() == 4194304:
                    return
        dlg = QFileDialog()
        fnm = dlg.getOpenFileName(self, 'Open project file...', filter='Pylot project (*.plp)')
        if fnm[0]:
            self.project = Project.load(fnm[0])
            self.init_events(new=True)
            if hasattr(self.project, 'metadata'):
                self.init_array_map(index=0)

    def saveProject(self):
        if self.project:
            if not self.project.location:
                self.createNewProject(exists=True)
            else:
                self.project.save()
            if not self.project.dirty:
                qmb = QMessageBox(icon=QMessageBox.Information, text='Saved back project to file:\n{}'.format(self.project.location))
                qmb.exec_()
                return
            else:
                # if still dirty because saving failed
                qmb = QMessageBox(icon=QMessageBox.Warning, text='Could not save back to original file.\n'
                                  'Choose new file')
                qmb.setStandardButtons(QMessageBox.Ok)
                qmb.exec_()
                self.createNewProject(exists=True)
                
    def draw(self):
        self.getPlotWidget().draw()

    def setDirty(self, value):
        self.dirty = value

    def closeEvent(self, event):
        if self.okToContinue():
            self.closing.emit()
            QMainWindow.closeEvent(self, event)

    def PyLoTprefs(self):
        props = PropertiesDlg(self, infile=self.getinfile())
        if props.exec_():
            return

    def helpHelp(self):
        if checkurl():
            form = HelpForm(
                'https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/wiki')
        else:
            form = HelpForm(':/help.html')
        form.show()


class Project(object):
    '''
    Pickable class containing information of a QtPyLoT project, like event lists and file locations.
    '''
    def __init__(self):
        self.eventlist = []
        self.location = None
        self.dirty = False
        self._table = None

    def add_eventlist(self, eventlist):
        if len(eventlist) == 0:
            return
        for item in eventlist:
            event = Event(item)
            if not event in self.eventlist:
                self.eventlist.append(event)
        self.setDirty()

    def setDirty(self):
        self.dirty = True

    def setClean(self):
        self.dirty = False

    def getEventFromPath(self, path):
        for event in self.eventlist:
            if event.path == path:
                return event

    def save(self, filename=None):
        '''
        Save PyLoT Project to a file. 
        Can be loaded by using project.load(filename).
        '''
        try:
            import cPickle
        except ImportError:
            import _pickle as cPickle

        if filename:
            self.location = filename
        else:
            filename = self.location

        try:
            outfile = open(filename, 'wb')
            cPickle.dump(self, outfile, -1)
            self.setClean()
        except Exception as e:
            print('Could not pickle PyLoT project. Reason: {}'.format(e))
            self.setDirty()

    @staticmethod
    def load(filename):
        try:
            import cPickle
        except ImportError:
            import _pickle as cPickle
        infile = open(filename, 'rb')
        project = cPickle.load(infile)
        print('Loaded %s' % filename)
        return project
    

class Event(object):
    '''
    Pickable class containing information on a single event.
    '''
    def __init__(self, path):
        self.path = path
        self.autopicks = None
        self.picks = None
        self.notes = None
        self._testEvent = False
        self._refEvent = False

    def addPicks(self, picks):
        self.picks = picks

    def addAutopicks(self, autopicks):
        self.autopicks = autopicks

    def addNotes(self, notes):
        self.notes = notes

    def clearNotes(self):
        self.notes = None

    def isRefEvent(self):
        return self._refEvent

    def isTestEvent(self):
        return self._testEvent

    def setRefEvent(self, bool):
        self._refEvent = bool
        if bool: self._testEvent = False

    def setTestEvent(self, bool):
        self._testEvent = bool
        if bool: self._refEvent = False
        
        
class getExistingDirectories(QFileDialog):
    def __init__(self, *args):
        super(getExistingDirectories, self).__init__(*args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setFileMode(self.Directory)
        self.setOption(self.ShowDirsOnly, True)
        self.findChildren(QListView)[0].setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.findChildren(QTreeView)[0].setSelectionMode(QAbstractItemView.ExtendedSelection)

        
def create_window():
    app_created = False
    app = QCoreApplication.instance()
    #check for existing app (when using ipython)
    if app is None:
        app = QApplication(sys.argv)
        app_created = True
    app.references = set()
    #app.references.add(window)
    #window.show()
    return app, app_created

        
def main():
    # create the Qt application
    pylot_app, app_created = create_window()
    #pylot_app = QApplication(sys.argv)
    pixmap = QPixmap(":/splash/splash.png")
    splash = QSplashScreen(pixmap)
    splash.show()

    app_icon = QIcon()
    app_icon.addPixmap(QPixmap(':/icons/pylot.png'))

    # create the main window
    pylot_form = MainWindow()
    icon = QIcon()
    pylot_form.setWindowIcon(icon)
    pylot_form.setIconSize(QSize(60, 60))
    
    splash.showMessage('Loading. Please wait ...')
    pylot_app.processEvents()

    # set Application Information
    pylot_app.setOrganizationName("Ruhr-University Bochum / BESTEC")
    pylot_app.setOrganizationDomain("rub.de")
    pylot_app.processEvents()
    pylot_app.setApplicationName("PyLoT")
    pylot_app.setApplicationVersion(pylot_form.__version__)
    pylot_app.setWindowIcon(app_icon)
    pylot_app.processEvents()

    # Show main window and run the app
    pylot_form.showMaximized()
    pylot_app.processEvents()

    splash.finish(pylot_form)
    if app_created:
        pylot_app.exec_()
    else:
        return pylot_form


if __name__ == "__main__":
    sys.exit(main())
