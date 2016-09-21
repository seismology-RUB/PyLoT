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

:author:
    Sebastian Wehling-Benatelli
:copyright:
    The PyLoT Development Team (https://ariadne.geophysik.rub.de/trac/PyLoT)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""

import matplotlib
import os
import sys

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from PySide.QtCore import QCoreApplication, QSettings, Signal, QFile, \
    QFileInfo, Qt
from PySide.QtGui import QMainWindow, QInputDialog, QIcon, QFileDialog, \
    QWidget, QHBoxLayout, QStyle, QKeySequence, QLabel, QFrame, QAction, \
    QDialog, QErrorMessage, QApplication, QPixmap, QMessageBox, QSplashScreen, \
    QActionGroup, QListWidget, QDockWidget, QLineEdit
import numpy as np
import subprocess
from obspy import UTCDateTime
from obspy.geodetics import degrees2kilometers
from obspy.core.event import Magnitude

from pylot.core.analysis.magnitude import calcsourcespec, calcMoMw
from pylot.core.io.data import Data
from pylot.core.io.inputs import FilterOptions, AutoPickParameter
from pylot.core.pick.autopick import autopickevent
from pylot.core.pick.compare import Comparison
from pylot.core.io.phases import picksdict_from_picks
import pylot.core.loc.nll as nll
from pylot.core.util.defaults import FILTERDEFAULTS, COMPNAME_MAP, \
    AUTOMATIC_DEFAULTS
from pylot.core.util.errors import FormatError, DatastructureError, \
    OverwriteError
from pylot.core.util.connection import checkurl
from pylot.core.util.dataprocessing import read_metadata
from pylot.core.util.utils import fnConstructor, getLogin, \
    full_range
from pylot.core.io.location import create_creation_info, create_event
from pylot.core.util.widgets import FilterOptionsDialog, NewEventDlg, \
    WaveformWidget, PropertiesDlg, HelpForm, createAction, PickDlg, \
    getDataType, ComparisonDialog
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.thread import AutoPickThread
from pylot.core.util.version import get_git_version as _getVersionString
import icons_rc

locateTool = dict(nll=nll)


class MainWindow(QMainWindow):
    __version__ = _getVersionString()
    closing = Signal()

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.createAction = createAction
        # read settings
        settings = QSettings()
        infile = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
        self._inputs = AutoPickParameter(infile)
        if settings.value("user/FullName", None) is None:
            fulluser = QInputDialog.getText(self, "Enter Name:", "Full name")
            settings.setValue("user/FullName", fulluser)
            settings.setValue("user/Login", getLogin())
        if settings.value("agency_id", None) is None:
            agency = QInputDialog.getText(self,
                                          "Enter authority name (e.g. BUG):",
                                          "Authority")
            settings.setValue("agency_id", agency)
        self.recentfiles = settings.value("data/recentEvents", [])
        self.fname = dict(manual=None, auto=None, loc=None)
        self.fnames = None
        self._stime = None
        structure_setting = settings.value("data/Structure", "PILOT")
        self.dataStructure = DATASTRUCTURE[structure_setting]()
        self.seismicPhase = str(settings.value("phase", "P"))
        self.dispComponent = str(settings.value("plotting/dispComponent", "Z"))
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

        # UI has to be set up before(!) children widgets are about to show up
        self.setupUi()

        # initialize event data
        if self.recentfiles:
            lastEvent = self.getLastEvent()
            self.data = Data(self, lastEvent)
        else:
            self.data = Data(self)
        self.autodata = Data(self)

        # load and display waveform data
        self.dirty = False
        self.load_data()
        if self.loadWaveformData():
            self.updateFilterOptions()
        else:
            sys.exit(0)

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
        _widget.setCursor(Qt.CrossCursor)
        _layout = QHBoxLayout()

        plottitle = "Overview: {0} components ".format(self.getComponent())

        # create central matplotlib figure canvas widget
        self.DataPlot = WaveformWidget(parent=self, xlabel=xlab, ylabel=None,
                                       title=plottitle)
        self.DataPlot.mpl_connect('button_press_event',
                                  self.pickOnStation)
        self.DataPlot.mpl_connect('axes_enter_event',
                                  lambda event: self.tutor_user())
        _layout.addWidget(self.DataPlot)

        manupicksicon = self.style().standardIcon(QStyle.SP_DialogYesButton)
        autopicksicon = self.style().standardIcon(QStyle.SP_DialogNoButton)
        locactionicon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
        loadpiloticon = self.style().standardIcon(QStyle.SP_ComputerIcon)
        quitIcon = self.style().standardIcon(QStyle.SP_MediaStop)
        saveIcon = self.style().standardIcon(QStyle.SP_DriveHDIcon)
        helpIcon = self.style().standardIcon(QStyle.SP_DialogHelpButton)
        newIcon = self.style().standardIcon(QStyle.SP_FileIcon)

        # create resource icons
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
        auto_icon = QIcon()
        auto_icon.addPixmap(QPixmap(':/icons/sync.png'))
        locate_icon = QIcon()
        locate_icon.addPixmap(QPixmap(':/icons/locate.png'))
        compare_icon = QIcon()
        compare_icon.addPixmap(QPixmap(':/icons/compare.png'))

        newEventAction = self.createAction(self, "&New event ...",
                                           self.createNewEvent,
                                           QKeySequence.New, newIcon,
                                           "Create a new event.")
        openmanualpicksaction = self.createAction(self, "Load &picks ...",
                                                  self.load_data,
                                                  QKeySequence.Open,
                                                  manupicksicon,
                                                  "Load pick data for "
                                                  "the actual event.")
        openmanualpicksaction.setData(None)
        openautopicksaction = self.createAction(self, "Load &automatic picks "
                                                      "...",
                                                self.load_autopicks,
                                                "Ctrl+A",
                                                autopicksicon,
                                                "Load automatic pick data "
                                                "for the actual event.")
        openautopicksaction.setData(None)

        loadlocationaction = self.createAction(self, "Load &location ...",
                                               self.load_loc, "Ctrl+L",
                                               locactionicon,
                                               "Load location information on "
                                               "the actual event.")
        loadpilotevent = self.createAction(self, "Load PILOT &event ...",
                                           self.load_pilotevent, "Ctrl+E",
                                           loadpiloticon,
                                           "Load PILOT event from information "
                                           "Matlab binary collections.")
        saveEventAction = self.createAction(self, "&Save event ...",
                                            self.saveData, QKeySequence.Save,
                                            saveIcon, "Save actual event data.")
        openWFDataAction = self.createAction(self, "Open &waveforms ...",
                                             self.loadWaveformData,
                                             "Ctrl+W", QIcon(":/wfIcon.png"),
                                             "Open waveform data (event will "
                                             "be closed)")
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
        printAction = self.createAction(self, "&Print event ...",
                                        self.show_event_information, QKeySequence.Print,
                                        print_icon,
                                        "Print waveform overview.")
        helpAction = self.createAction(self, "&Help ...", self.helpHelp,
                                       QKeySequence.HelpContents, helpIcon,
                                       """Show either the documentation
                                       homepage (internet connection available),
                                       or shipped documentation files.""")
        self.fileMenu = self.menuBar().addMenu('&File')
        self.fileMenuActions = (newEventAction, openmanualpicksaction,
                                saveEventAction, openWFDataAction, None,
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
        fileToolActions = (newEventAction, openmanualpicksaction,
                           openautopicksaction, loadlocationaction,
                           loadpilotevent, saveEventAction)
        fileToolBar.setObjectName("FileTools")
        self.addActions(fileToolBar, fileToolActions)

        # phaseToolBar = self.addToolBar("PhaseTools")
        # phaseToolActions = (self.selectPAction, self.selectSAction)
        # phaseToolBar.setObjectName("PhaseTools")
        # self.addActions(phaseToolBar, phaseToolActions)

        # create button group for component selection

        componentGroup = QActionGroup(self)
        componentGroup.setExclusive(True)

        z_action = self.createAction(parent=componentGroup, text='Z',
                                     slot=self.plotZ, shortcut='Alt+Z',
                                     icon=z_icon, tip='Display the vertical (Z)'
                                                      ' component.',
                                     checkable=True)
        z_action.setChecked(True)

        n_action = self.createAction(parent=componentGroup, text='N',
                                     slot=self.plotN, shortcut='Alt+N',
                                     icon=n_icon,
                                     tip='Display the north-south (N) '
                                         'component.', checkable=True)

        e_action = self.createAction(parent=componentGroup, text='E',
                                     slot=self.plotE, shortcut='Alt+E',
                                     icon=e_icon,
                                     tip='Display the east-west (E) component.',
                                     checkable=True)

        componentToolBar = self.addToolBar("ComponentSelection")
        componentActions = (z_action, n_action, e_action)
        componentToolBar.setObjectName("PhaseTools")
        self.addActions(componentToolBar, componentActions)

        auto_pick = self.createAction(parent=self, text='autoPick',
                                      slot=self.autoPick, shortcut='Alt+Ctrl+A',
                                      icon=auto_icon, tip='Automatically pick'
                                                          ' the entire dataset'
                                                          ' displayed!')

        autoPickToolBar = self.addToolBar("autoPyLoT")
        autoPickActions = (auto_pick, self.compare_action)
        self.addActions(autoPickToolBar, autoPickActions)

        # pickToolBar = self.addToolBar("PickTools")
        # pickToolActions = (selectStation, )
        # pickToolBar.setObjectName("PickTools")
        # self.addActions(pickToolBar, pickToolActions)

        locateEvent = self.createAction(parent=self, text='locate_event',
                                        slot=self.locate_event,
                                        shortcut='Alt+Ctrl+L',
                                        icon=locate_icon,
                                        tip='Locate the event using '
                                            'the picked arrivals.')

        locationToolBar = self.addToolBar("LocationTools")
        locationToolActions = (locateEvent,)
        locationToolBar.setObjectName("LocationTools")
        self.addActions(locationToolBar, locationToolActions)

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 500)

        _widget.setLayout(_layout)
        _widget.showFullScreen()

        self.setCentralWidget(_widget)

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

    def get_fnames(self):
        return self.fname

    def set_fname(self, fname, type):
        if self.get_fnames()[type] is not None:
            self.add_recentfile(self.get_fnames())
        self.fname[type] = fname

    def getEventFileName(self):
        if self.get_fnames() is None:
            self.set_fname(self.get_data().getEventFileName())
        return self.get_fnames()

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
        if not fbasename:
            return False
        self.get_data().exportEvent(fbasename, exform)
        self.setDirty(False)
        self.update_status('Event saved as %s' % (fbasename + exform))
        return True

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
        return self.DataPlot

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

    def loadWaveformData(self):
        if self.fnames and self.okToContinue():
            self.setDirty(True)
            ans = self.data.setWFData(self.fnames)
        elif self.fnames is None and self.okToContinue():
            ans = self.data.setWFData(self.getWFFnames())
        else:
            ans = False
        self._stime = full_range(self.get_data().getWFData())[0]
        if ans:
            self.plotWaveformData()
            return ans
        else:
            return ans

    def plotWaveformData(self):
        zne_text = {'Z': 'vertical', 'N': 'north-south', 'E': 'east-west'}
        comp = self.getComponent()
        title = 'section: {0} components'.format(zne_text[comp])
        alter_comp = COMPNAME_MAP[comp]
        wfst = self.get_data().getWFData().select(component=comp)
        wfst += self.get_data().getWFData().select(component=alter_comp)
        self.getPlotWidget().plotWFData(wfdata=wfst, title=title, mapping=False)
        self.draw()
        plotDict = self.getPlotWidget().getPlotDict()
        pos = plotDict.keys()
        labels = [plotDict[n][0] for n in pos]
        self.getPlotWidget().setYTickLabels(pos, labels)

    def plotZ(self):
        self.setComponent('Z')
        self.plotWaveformData()
        self.drawPicks()
        self.draw()

    def plotN(self):
        self.setComponent('N')
        self.plotWaveformData()
        self.drawPicks()
        self.draw()

    def plotE(self):
        self.setComponent('E')
        self.plotWaveformData()
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
        pickDlg = PickDlg(self, data=data.select(station=station),
                          station=station,
                          picks=self.getPicksOnStation(station))
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
            self.set_loc_flag(True)
        elif self.get_loc_flag() and not self.check4Loc():
            self.set_loc_flag(False)

    def addListItem(self, text):
        self.listWidget.addItem(text)
        self.listWidget.scrollToBottom()

    def autoPick(self):
        self.listWidget = QListWidget()
        self.setDirty(True)
        self.logDockWidget = QDockWidget("AutoPickLog", self)
        self.logDockWidget.setObjectName("LogDockWidget")
        self.logDockWidget.setAllowedAreas(
            Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.logDockWidget.setWidget(self.listWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.logDockWidget)
        self.addListItem('loading default values for local data ...')
        # may become obsolete if generalized input parameter a read from disc
        #  during initialization
        # TODO double check for obsolete read in of parameters
        autopick_parameter = AutoPickParameter(AUTOMATIC_DEFAULTS)
        self.addListItem(str(autopick_parameter))

        # Create the worker thread and run it
        self.thread = AutoPickThread(parent=self,
                                     func=autopickevent,
                                     data=self.get_data().getWFData(),
                                     param=autopick_parameter)
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
            self.picks.update(picks)
        elif type == 'auto':
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
            'P': ('c', 'c--', 'b-', 'bv', 'b^'),
            'S': ('m', 'm--', 'r-', 'rv', 'r^')
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

            if picktype == 'manual':
                ax.fill_between([epp, lpp], ylims[0], ylims[1],
                                alpha=.5, color=colors[0])
                ax.plot([mpp - spe, mpp - spe], ylims, colors[1],
                        [mpp, mpp], ylims, colors[2],
                        [mpp + spe, mpp + spe], ylims, colors[1])
            elif picktype == 'auto':
                ax.plot(mpp, ylims[1], colors[3],
                        mpp, ylims[0], colors[4])
            else:
                raise TypeError('Unknow picktype {0}'.format(picktype))

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
        self.get_data().get_evt_data().magnitudes.append(self.calc_magnitude())

    def calc_magnitude(self):
        e = self.get_data().get_evt_data()
        settings = QSettings()
        if e.origins:
            o = e.origins[0]
            mags = dict()
            fninv = settings.value("inventoryFile", None)
            if fninv is None:
                fninv, _ = QFileDialog.getOpenFileName(self, self.tr(
                    "Select inventory..."), self.tr("Select file"))
                ans = QMessageBox.question(self, self.tr("Make default..."),
                                           self.tr(
                                               "New inventory filename set.\n" + \
                                               "Do you want to make it the default value?"),
                                           QMessageBox.Yes | QMessageBox.No,
                                           QMessageBox.No)
                if ans == QMessageBox.Yes:
                    settings.setValue("inventoryFile", fninv)
                    settings.sync()
            metadata = read_metadata(fninv)
            for a in o.arrivals:
                if a.phase in 'sS':
                    continue
                pick = a.pick_id.get_referred_object()
                station = pick.waveform_id.station_code
                wf = self.get_data().getWFData().select(station=station)
                if not wf:
                    continue
                onset = pick.time
                dist = degrees2kilometers(a.distance)
                w0, fc = calcsourcespec(wf, onset, metadata, self.inputs.get('vp'), dist,
                                        a.azimuth, a.takeoff_angle,
                                        self.inputs.get('Qp'), 0)
                if w0 is None or fc is None:
                    continue
                station_mag = calcMoMw(wf, w0, self.inputs.get('rho'),
                                       self.inputs.get('vp'), dist)
                mags[station] = station_mag
            mag = np.median([M[1] for M in mags.values()])
            # give some information on the processing
            print('number of stations used: {0}\n'.format(len(mags.values())))
            print('stations used:\n')
            for s in mags.keys(): print('\t{0}'.format(s))

            return Magnitude(mag=mag, magnitude_type='Mw')
        else:
            return None

    def check4Loc(self):
        return self.picksNum() > 4

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
                self.setWindowTitle("PyLoT - New event [*]")
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

    def draw(self):
        self.getPlotWidget().draw()

    def setDirty(self, value):
        self.dirty = value

    def closeEvent(self, event):
        if self.okToContinue():
            self.closing.emit()
            QMainWindow.closeEvent(self, event)

    def PyLoTprefs(self):
        props = PropertiesDlg(self)
        if props.exec_():
            return

    def helpHelp(self):
        if checkurl():
            form = HelpForm(
                'https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/wiki')
        else:
            form = HelpForm(':/help.html')
        form.show()


def main():
    # create the Qt application
    pylot_app = QApplication(sys.argv)
    pixmap = QPixmap(":/splash/splash.png")
    splash = QSplashScreen(pixmap)
    splash.show()

    app_icon = QIcon()
    app_icon.addPixmap(QPixmap(':/icons/pylot.png'))

    # create the main window
    pylot_form = MainWindow()
    splash.showMessage('Loading. Please wait ...')
    pylot_app.processEvents()

    # set Application Information
    pylot_app.setOrganizationName("Ruhr-University Bochum / MAGS2")
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
    pylot_app.exec_()


if __name__ == "__main__":
    sys.exit(main())
