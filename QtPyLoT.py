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
    Sebastian Wehling-Benatelli / Ludger Küperkoch / Marcel Paffrath
:copyright:
    The PyLoT Development Team (https://ariadne.geophysik.rub.de/trac/PyLoT)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""

import os
import sys
import platform
import argparse
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
from obspy.core.event import Magnitude, Origin
from obspy.core.util import AttribDict

try:
    import pyqtgraph as pg
except Exception as e:
    print('QtPyLoT: Could not import pyqtgraph. {}'.format(e))
    pg = None

try:
    from matplotlib.backends.backend_qt4agg import FigureCanvas
except ImportError:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from pylot.core.analysis.magnitude import LocalMagnitude, MomentMagnitude
from pylot.core.io.data import Data
from pylot.core.io.inputs import FilterOptions, PylotParameter
from autoPyLoT import autoPyLoT
from pylot.core.pick.compare import Comparison
from pylot.core.pick.utils import symmetrize_error
from pylot.core.io.phases import picksdict_from_picks
import pylot.core.loc.nll as nll
from pylot.core.util.defaults import FILTERDEFAULTS, OUTPUTFORMATS, SetChannelComponents, \
    readFilterInformation
from pylot.core.util.errors import FormatError, DatastructureError, \
    OverwriteError, ProcessingError
from pylot.core.util.connection import checkurl
from pylot.core.util.dataprocessing import read_metadata, restitute_data
from pylot.core.util.utils import fnConstructor, getLogin, \
    full_range
from pylot.core.util.event import Event
from pylot.core.io.location import create_creation_info, create_event
from pylot.core.util.widgets import FilterOptionsDialog, NewEventDlg, \
    WaveformWidget, WaveformWidgetPG, PropertiesDlg, HelpForm, createAction, PickDlg, \
    getDataType, ComparisonDialog, TuneAutopicker, PylotParaBox
from pylot.core.util.map_projection import map_projection
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.thread import AutoPickThread, Thread
from pylot.core.util.version import get_git_version as _getVersionString

if sys.version_info.major == 3:
    import icons_rc_3 as icons_rc
elif sys.version_info.major == 2:
    import icons_rc_2 as icons_rc
else:
    raise ImportError('Could not determine python version.')

locateTool = dict(nll=nll)


class MainWindow(QMainWindow):
    __version__ = _getVersionString()
    closing = Signal()

    def __init__(self, parent=None, infile=None):
        super(MainWindow, self).__init__(parent)

        # check for default pylot.in-file
        if not infile:
            infile = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
            print('Using default input file {}'.format(infile))
        if os.path.isfile(infile) == False:
            infile = QFileDialog().getOpenFileName(caption='Choose PyLoT-input file')
            if not os.path.exists(infile[0]):
                QMessageBox.warning(self, "PyLoT Warning",
                                    "No PyLoT-input file declared!")
                sys.exit(0)
            self.infile = infile[0]
        else:
            self.infile = infile
        self._inputs = PylotParameter(infile)
        self._props = None

        self.dirty = False
        self.project = Project()
        self.project.parameter = self._inputs
        self.tap = None
        self.paraBox = None
        self.array_map = None
        self._metadata = None
        self._eventChanged = [False, False]

        self.poS_id = None
        self.ae_id = None
        self.scroll_id = None
        self._ctrl = False  # control key pressed
        self._shift = False  # shift key pressed

        # default factor for dataplot e.g. enabling/disabling scrollarea
        self.height_factor = 12

        # default colors for ref/test event
        self._colors = {
            'ref': QtGui.QColor(200, 210, 230, 255),
            'test': QtGui.QColor(200, 230, 200, 255)
        }

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
            self.data._new = False
        self.autodata = Data(self)

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
        if settings.value('compclass', None) is None:
            settings.setValue('compclass', SetChannelComponents())
        if settings.value('output/Format', None) is None:
            outformat = QInputDialog.getText(self,
                                             "Enter output format (*.xml, *.cnv, *.obs):",
                                             "Format")
            settings.setValue("output/Format", outformat)
        settings.sync()

        # setup UI
        self.setupUi()

        filter_info = readFilterInformation(self._inputs)
        p_filter = filter_info['P']
        s_filter = filter_info['S']
        self.filteroptions = {'P': FilterOptions(p_filter['filtertype'],
                                                 p_filter['freq'],
                                                 p_filter['order']),
                              'S': FilterOptions(s_filter['filtertype'],
                                                 s_filter['freq'],
                                                 s_filter['order'])}
        self.pylot_picks = {}
        self.pylot_autopicks = {}
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

        _widget = QWidget()
        self._main_layout = QVBoxLayout()

        # add event combo box and ref/test buttons
        self.eventBox = self.createEventBox()
        self.eventBox.setMaxVisibleItems(30)
        self.eventBox.setEnabled(False)
        self.init_ref_test_buttons()
        self._event_layout = QHBoxLayout()
        self._event_layout.addWidget(QLabel('Event: '))
        self._event_layout.addWidget(self.eventBox)
        self._event_layout.addWidget(self.ref_event_button)
        self._event_layout.addWidget(self.test_event_button)
        self._event_layout.setStretch(1, 1)  # set stretch of item 1 to 1
        self._main_layout.addLayout(self._event_layout)
        self.eventBox.activated.connect(self.refreshEvents)

        # add main tab widget
        self.tabs = QTabWidget()
        self._main_layout.addWidget(self.tabs)
        self.tabs.currentChanged.connect(self.refreshTabs)

        # add scroll area used in case number of traces gets too high
        self.wf_scroll_area = QtGui.QScrollArea()

        # create central matplotlib figure canvas widget
        self.pg = pg
        self.init_wfWidget()

        # init main widgets for main tabs
        wf_tab = QtGui.QWidget()
        array_tab = QtGui.QWidget()
        events_tab = QtGui.QWidget()

        # init main widgets layouts
        self.wf_layout = QtGui.QVBoxLayout()
        self.array_layout = QtGui.QVBoxLayout()
        self.events_layout = QtGui.QVBoxLayout()
        wf_tab.setLayout(self.wf_layout)
        array_tab.setLayout(self.array_layout)
        events_tab.setLayout(self.events_layout)

        # add tabs to main tab widget
        self.tabs.addTab(wf_tab, 'Waveform Plot')
        self.tabs.addTab(array_tab, 'Array Map')
        self.tabs.addTab(events_tab, 'Eventlist')

        self.wf_layout.addWidget(self.wf_scroll_area)
        self.wf_scroll_area.setWidgetResizable(True)
        self.init_array_tab()
        self.init_event_table()
        self.tabs.setCurrentIndex(0)

        quitIcon = self.style().standardIcon(QStyle.SP_MediaStop)
        helpIcon = self.style().standardIcon(QStyle.SP_DialogHelpButton)
        newFolderIcon = self.style().standardIcon(QStyle.SP_FileDialogNewFolder)

        # create resource icons
        newIcon = QIcon()
        newIcon.addPixmap(QPixmap(':/icons/newfile.png'))
        addIcon = QIcon()
        addIcon.addPixmap(QPixmap(':/icons/add.png'))
        saveIcon = QIcon()
        saveIcon.addPixmap(QPixmap(':/icons/save.png'))
        saveasIcon = QIcon()
        saveasIcon.addPixmap(QPixmap(':/icons/saveas.png'))
        saveProjectIcon = QIcon()
        saveProjectIcon.addPixmap(QPixmap(':/icons/saveproject.png'))
        saveProjectAsIcon = QIcon()
        saveProjectAsIcon.addPixmap(QPixmap(':/icons/saveprojectas.png'))
        openIcon = QIcon()
        openIcon.addPixmap(QPixmap(':/icons/openfile.png'))
        openProjectIcon = QIcon()
        openProjectIcon.addPixmap(QPixmap(':/icons/openproject.png'))
        openLocIcon = QIcon()
        openLocIcon.addPixmap(QPixmap(':/icons/openloc.png'))
        openEventIcon = QIcon()
        openEventIcon.addPixmap(QPixmap(':/icons/openpick.png'))
        openEventsIcon = QIcon()
        openEventsIcon.addPixmap(QPixmap(':/icons/openpicks.png'))
        saveEventsIcon = QIcon()
        saveEventsIcon.addPixmap(QPixmap(':/icons/savepicks.png'))
        prefIcon = QIcon()
        prefIcon.addPixmap(QPixmap(':/icons/preferences.png'))
        paraIcon = QIcon()
        paraIcon.addPixmap(QPixmap(':/icons/parameter.png'))
        self.autopicksicon_small = QIcon()
        self.autopicksicon_small.addPixmap(QPixmap(':/icons/autopicksicon_small.png'))
        self.manupicksicon_small = QIcon()
        self.manupicksicon_small.addPixmap(QPixmap(':/icons/manupicksicon_small.png'))
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
        autotune_icon = QIcon()
        autotune_icon.addPixmap(QPixmap(':/icons/tune.png'))
        autopylot_icon = QIcon()
        autopylot_icon.addPixmap(QPixmap(':/icons/autopylot_button'))
        locate_icon = QIcon()
        locate_icon.addPixmap(QPixmap(':/icons/locate_button.png'))
        compare_icon = QIcon()
        compare_icon.addPixmap(QPixmap(':/icons/compare_button.png'))
        self.newProjectAction = self.createAction(self, "&New project ...",
                                                  self.createNewProject,
                                                  QKeySequence.New, newIcon,
                                                  "Create a new Project.")
        self.openProjectAction = self.createAction(self, "&Open project ...",
                                                   self.loadProject,
                                                   QKeySequence.Open,
                                                   openProjectIcon,
                                                   "Load project file")
        self.saveProjectAction = self.createAction(self, "&Save project ...",
                                                   self.saveProject,
                                                   QKeySequence.Save,
                                                   saveProjectIcon,
                                                   "Save project file")
        self.saveProjectAction.setEnabled(False)
        self.saveProjectAsAction = self.createAction(self, "Save project as ...",
                                                     self.saveProjectAs,
                                                     QKeySequence.SaveAs,
                                                     saveProjectAsIcon,
                                                     "Save project file as...")
        self.saveProjectAsAction.setEnabled(False)
        # newEventAction = self.createAction(self, "&New event ...",
        #                                    self.createNewEvent,
        #                                    QKeySequence.New, newIcon,
        #                                    "Create a new event.")
        self.openmanualpicksaction = self.createAction(self, "Load event ...",
                                                       self.load_data,
                                                       "Ctrl+M",
                                                       openEventIcon,
                                                       "Load event information for "
                                                       "the displayed event.")
        self.openmanualpicksaction.setEnabled(False)
        self.openmanualpicksaction.setData(None)

        self.openautopicksaction = self.createAction(self, "Load event information &automatically ... ",
                                                     self.load_multiple_data,
                                                     "Ctrl+A",
                                                     openEventsIcon,
                                                     "Load event data automatically "
                                                     "for for all events.")
        self.openautopicksaction.setEnabled(False)
        self.openautopicksaction.setData(None)

        self.loadlocationaction = self.createAction(self, "Load &location ...",
                                                    self.load_loc, "Ctrl+L",
                                                    openLocIcon,
                                                    "Load location information on "
                                                    "the displayed event.")
        self.loadlocationaction.setEnabled(False)
        self.loadpilotevent = self.createAction(self, "Load PILOT &event ...",
                                                self.load_pilotevent, None,
                                                loadpiloticon,
                                                "Load PILOT event from information "
                                                "MatLab binary collections (created"
                                                " in former MatLab based version).")
        self.loadpilotevent.setEnabled(False)

        self.saveManualPicksAction = self.createAction(self, "Save &picks ...",
                                                       self.saveData, "Ctrl+P",
                                                       saveEventsIcon, "Save event pick data.")
        self.disableSaveManualPicksAction()

        self.addEventDataAction = self.createAction(self, "Add &events ...",
                                                    self.add_events,
                                                    "Ctrl+E", addIcon,
                                                    "Add event data")
        prefsEventAction = self.createAction(self, "Preferences",
                                             self.PyLoTprefs,
                                             QKeySequence.Preferences,
                                             prefIcon,
                                             "Edit PyLoT app preferences.")
        quitAction = self.createAction(self, "&Quit",
                                       QCoreApplication.instance().quit,
                                       QKeySequence.Close, quitIcon,
                                       "Close event and quit PyLoT")
        self.parameterAction = self.createAction(self, "Parameter",
                                                 self.setParameter,
                                                 None, paraIcon,
                                                 "Modify Parameter")
        self.filterAction = self.createAction(self, "&Filter ...",
                                              self.filterWaveformData,
                                              "Ctrl+F", filter_icon,
                                              """Toggle un-/filtered waveforms
                                              to be displayed, according to the
                                              desired seismic phase.""", True)
        filterEditAction = self.createAction(self, "&Filter parameter ...",
                                             self.adjustFilterOptions,
                                             "Alt+F", filter_icon,
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
                                self.saveProjectAsAction,
                                self.openmanualpicksaction, self.saveManualPicksAction, None,
                                prefsEventAction, self.parameterAction, quitAction)
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
                           self.saveProjectAsAction)
        fileToolBar.setObjectName("FileTools")
        self.addActions(fileToolBar, fileToolActions)
        eventToolActions = (self.openmanualpicksaction, self.openautopicksaction,
                            self.saveManualPicksAction, self.loadlocationaction,
                            self.loadpilotevent)
        eventToolBar = self.addToolBar("EventTools")
        eventToolBar.setObjectName("EventTools")
        self.addActions(eventToolBar, eventToolActions)

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

        self.auto_tune = self.createAction(parent=self, text='autoTune',
                                           slot=self.tune_autopicker, shortcut='Alt+Ctrl+A',
                                           icon=autotune_icon, tip='Tune autopicking algorithm.')

        self.auto_tune.setEnabled(False)

        self.auto_pick = self.createAction(parent=self, text='autoPick',
                                           slot=self.autoPick, shortcut='Alt+Ctrl+A',
                                           icon=autopylot_icon, tip='Automatically pick'
                                                                    ' the displayed waveforms.')
        self.auto_pick.setEnabled(False)

        autoPickToolBar = self.addToolBar("autoPyLoT")
        autoPickActions = (self.auto_tune, self.auto_pick, self.compare_action)
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

    def init_wfWidget(self):
        settings = QSettings()
        xlab = self.startTime.strftime('seconds since %Y/%m/%d %H:%M:%S (%Z)')
        plottitle = None  # "Overview: {0} components ".format(self.getComponent())
        self.disconnectWFplotEvents()
        if str(settings.value('pyqtgraphic')) == 'false' or not pg:
            self.pg = False
            self.dataPlot = WaveformWidget(parent=self, xlabel=xlab, ylabel=None,
                                           title=plottitle)
        else:
            self.pg = True
            self.dataPlot = WaveformWidgetPG(parent=self, xlabel=xlab, ylabel=None,
                                             title=plottitle)
        self.dataPlot.setCursor(Qt.CrossCursor)
        self.wf_scroll_area.setWidget(self.dataPlot)
        if self.get_current_event():
            self.plotWaveformDataThread()

    def init_ref_test_buttons(self):
        '''
        Initiate/create buttons for assigning events containing manual picks to reference or test set.
        '''
        self.ref_event_button = QtGui.QPushButton('Ref')
        self.test_event_button = QtGui.QPushButton('Test')
        self.ref_event_button.setToolTip('Set manual picks of current ' +
                                         'event as reference picks for autopicker tuning.')
        self.test_event_button.setToolTip('Set manual picks of current ' +
                                          'event as test picks for autopicker testing.')
        self.ref_event_button.setCheckable(True)
        self.test_event_button.setCheckable(True)
        self.set_button_color(self.ref_event_button, self._colors['ref'])
        self.set_button_color(self.test_event_button, self._colors['test'])
        self.ref_event_button.clicked.connect(self.toggleRef)
        self.test_event_button.clicked.connect(self.toggleTest)
        self.ref_event_button.setEnabled(False)
        self.test_event_button.setEnabled(False)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key.Key_Control:
            self._ctrl = True
        if event.key() == QtCore.Qt.Key.Key_Shift:
            self._shift = True

    def keyReleaseEvent(self, event):
        if event.key() == QtCore.Qt.Key.Key_Control:
            self._ctrl = False
        if event.key() == QtCore.Qt.Key.Key_Shift:
            self._shift = False

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

    def load_multiple_data(self, type='manual'):
        if not self.okToContinue():
            return
        refresh = False
        events = self.project.eventlist
        fext = '.xml'
        for event in events:
            path = event.path
            eventname = path.split('/')[-1]
            filename = os.path.join(path, 'PyLoT_' + eventname + fext)
            if os.path.isfile(filename):
                self.load_data(filename, draw=False, event=event, overwrite=True)
                refresh = True
        if not refresh:
            return
        if self.get_current_event().pylot_picks:
            self.refreshEvents()
        self.setDirty(True)

    def load_data(self, fname=None, type='manual', loc=False, draw=True, event=None, overwrite=False):
        if not overwrite:
            if not self.okToContinue():
                return
        if fname is None:
            action = self.sender()
            if isinstance(action, QAction):
                fname = self.filename_from_action(action)
                if not fname:
                    return
        self.set_fname(fname, type)
        # data = dict(auto=self.autodata, manual=self.data)
        if not event:
            event = self.get_current_event()
        data = Data(self, event)
        try:
            data_new = Data(self, evtdata=str(fname))
            data += data_new
        except ValueError:
            qmb = QMessageBox(self, icon=QMessageBox.Question,
                              text='Warning: Missmatch in event identifiers {} and {}. Continue?'.format(
                                  data_new.get_evt_data().resource_id,
                                  data.get_evt_data().resource_id),
                              windowTitle='PyLoT - Load data warning')
            qmb.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            qmb.setDefaultButton(QMessageBox.No)
            ret = qmb.exec_()
            if ret == qmb.Yes:
                data_new.setNew()
                data += data_new
            else:
                return

        self.data = data
        print('Loading {} picks from file {}.'.format(type, fname))
        if not loc:
            self.updatePicks(type=type, event=event)
        if draw:
            if self.get_current_event().pylot_picks:
                self.refreshEvents()
            self.setDirty(True)

    def add_recentfile(self, event):
        self.recentfiles.insert(0, event)

    def set_button_color(self, button, color=None):
        '''
        Set background color of a button.
        button: type = QtGui.QAbstractButton
        color: type = QtGui.QColor or type = str (RGBA)
        '''
        if type(color) == QtGui.QColor:
            palette = button.palette()
            role = button.backgroundRole()
            palette.setColor(role, color)
            button.setPalette(palette)
            button.setAutoFillBackground(True)
        elif type(color) == str or not color:
            button.setStyleSheet("background-color: {}".format(color))

    def getWFFnames(self):
        try:
            evt = self.get_data().get_evt_data()
            if evt.pylot_picks:
                for pick in evt.pylot_picks:
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
            if not self._props:
                self._props = PropertiesDlg(self, infile=self.infile)

            if self._props.exec_() == QDialog.Accepted:
                return self.getWFFnames()
            else:
                return

    def getWFFnames_from_eventbox(self, eventbox=None):
        '''
        Return waveform filenames from event in eventbox.
        '''
        if self.dataStructure:
            directory = self.get_current_event_path(eventbox)
            if not directory:
                return
            fnames = [os.path.join(directory, f) for f in os.listdir(directory)]
        else:
            raise DatastructureError('not specified')
        return fnames

    def get_current_event(self, eventbox=None):
        '''
        Return event (type QtPylot.Event) currently selected in eventbox.
        '''
        if not eventbox:
            eventbox = self.eventBox
        path = eventbox.currentText()
        return self.project.getEventFromPath(path)

    def get_current_event_path(self, eventbox=None):
        '''
        Return event path of event (type QtPylot.Event) currently selected in eventbox.
        '''
        event = self.get_current_event(eventbox)
        if event:
            return event.path

    def get_current_event_name(self, eventbox=None):
        '''
        Return event path of event (type QtPylot.Event) currently selected in eventbox.
        '''
        path = self.get_current_event_path(eventbox)
        if path:
            return path.split('/')[-1]

    def getLastEvent(self):
        return self.recentfiles[0]

    def add_events(self):
        '''
        Creates and adds events by user selection of event folders to GUI.
        '''
        if not self.project:
            self.createNewProject()
        ed = getExistingDirectories(self, 'Select event directories...')
        if ed.exec_():
            eventlist = ed.selectedFiles()
            # select only folders that start with 'e', containin two dots and have length 12
            eventlist = [item for item in eventlist if item.split('/')[-1].startswith('e')
                         and len(item.split('/')[-1].split('.')) == 3
                         and len(item.split('/')[-1]) == 12]
            if not eventlist:
                print('No events found! Expected structure for event folders: [evID.DOY.YR]')
                return
        else:
            return
        if not self.project:
            print('No project found.')
            return

        # get path from first event in list and split them
        path = eventlist[0]
        try:
            system_name = platform.system()
            if system_name in ["Linux", "Darwin"]:
                dirs = {
                    'database': path.split('/')[-2],
                    'datapath': path.split('/')[-3],
                    'rootpath': '/' + os.path.join(*path.split('/')[:-3])
                }
            elif system_name == "Windows":
                rootpath = path.split('/')[:-3]
                rootpath[0] += '/'
                dirs = {
                    'database': path.split('/')[-2],
                    'datapath': path.split('/')[-3],
                    'rootpath': os.path.join(*rootpath)
                }
        except Exception as e:
            dirs = {
                'database': '',
                'datapath': '',
                'rootpath': ''
            }
            print('Warning: Could not automatically init folder structure. ({})'.format(e))

        settings = QSettings()
        settings.setValue("data/dataRoot", dirs['rootpath'])
        settings.sync()

        if not self.project.eventlist:
            # init parameter object
            self.setParameter(show=False)
            # hide all parameter (show all needed parameter later)
            self.paraBox.hide_parameter()
            for directory in dirs.keys():
                # set parameter
                box = self.paraBox.boxes[directory]
                self.paraBox.setValue(box, dirs[directory])
                # show needed parameter in box
                self.paraBox.show_parameter(directory)
            dirs_box = self.paraBox.get_groupbox_dialog('Directories')
            if not dirs_box.exec_():
                return
            self.project.rootpath = dirs['rootpath']
        else:
            if hasattr(self.project, 'rootpath'):
                if not self.project.rootpath == dirs['rootpath']:
                    QMessageBox.warning(self, "PyLoT Warning",
                                        'Rootpath missmatch to current project!')
                    return
            else:
                self.project.rootpath = dirs['rootpath']

        self.project.add_eventlist(eventlist)
        self.init_events()
        self.setDirty(True)

    def createEventBox(self):
        '''
        Eventbox generator.
        '''
        qcb = QComboBox()
        palette = qcb.palette()
        # change highlight color:
        # palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Highlight,
        #                  QtGui.QBrush(QtGui.QColor(0,0,127,127)))
        # palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Highlight,
        #                  QtGui.QBrush(QtGui.QColor(0,0,127,127)))
        qcb.setPalette(palette)
        return qcb

    def init_events(self, new=False):
        '''
        Initiate GUI widgets in case of changed or newly added events.
        '''
        nitems = self.eventBox.count()
        if len(self.project.eventlist) == 0:
            print('No events to init.')
            self.clearWaveformDataPlot()
            return
        self.eventBox.setEnabled(True)
        self.fill_eventbox()
        if new:
            self.eventBox.setCurrentIndex(0)
        else:
            self.eventBox.setCurrentIndex(nitems)
        self.refreshEvents()
        tabindex = self.tabs.currentIndex()

    def fill_eventbox(self, event=None, eventBox=None, select_events='all'):
        '''
        (Re)fill the selected eventBox (type = QtGui.QComboBox).

        :param: select_events, can be 'all', 'ref'
        :type: str
        '''
        if not eventBox:
            eventBox = self.eventBox
        index = eventBox.currentIndex()
        tv = QtGui.QTableView()
        header = tv.horizontalHeader()
        header.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)
        header.hide()
        tv.verticalHeader().hide()
        tv.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)

        eventBox.setView(tv)
        eventBox.clear()
        model = eventBox.model()
        plmax = 0
        # set maximum length of path string
        for event in self.project.eventlist:
            pl = len(event.path)
            if pl > plmax:
                plmax = pl

        for id, event in enumerate(self.project.eventlist):
            event_path = event.path
            event_npicks = 0
            event_nautopicks = 0
            if event.pylot_picks:
                event_npicks = len(event.pylot_picks)
            if event.pylot_autopicks:
                event_nautopicks = len(event.pylot_autopicks)
            event_ref = event.isRefEvent()
            event_test = event.isTestEvent()

            # text = '{path:{plen}} | manual: [{p:3d}] | auto: [{a:3d}]'
            # text = text.format(path=event_path,
            #                    plen=plmax,
            #                    p=event_npicks,
            #                    a=event_nautopicks)

            item_path = QtGui.QStandardItem('{path:{plen}}'.format(path=event_path, plen=plmax))
            item_nmp = QtGui.QStandardItem(str(event_npicks))
            item_nmp.setIcon(self.manupicksicon_small)
            item_nap = QtGui.QStandardItem(str(event_nautopicks))
            item_nap.setIcon(self.autopicksicon_small)
            item_ref = QtGui.QStandardItem()  # str(event_ref))
            item_test = QtGui.QStandardItem()  # str(event_test))
            if event_ref:
                item_ref.setBackground(self._colors['ref'])
            if event_test:
                item_test.setBackground(self._colors['test'])
            item_notes = QtGui.QStandardItem(event.notes)

            openIcon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
            item_path.setIcon(openIcon)
            # if ref: set different color e.g.
            # if event_ref:
            #     item.setBackground(self._colors['ref'])
            # if event_test:
            #     itemt.setBackground(self._colors['test'])
            # item.setForeground(QtGui.QColor('black'))
            # font = item.font()
            # font.setPointSize(10)
            # item.setFont(font)
            # item2.setForeground(QtGui.QColor('black'))
            # item2.setFont(font)
            itemlist = [item_path, item_nmp, item_nap, item_ref, item_test, item_notes]
            if event_test and select_events == 'ref':
                for item in itemlist:
                    item.setEnabled(False)
            model.appendRow(itemlist)
            if not event.path == self.eventBox.itemText(id).strip():
                message = ('Path missmatch creating eventbox.\n'
                           '{} unequal {}.'
                           .format(event.path, self.eventBox.itemText(id)))
                raise ValueError(message)
                # not working with obspy events
                # eventBox.setItemData(id, event)
        eventBox.setCurrentIndex(index)
        self.refreshRefTestButtons()

    def filename_from_action(self, action):
        if action.data() is None:
            filt = "Supported file formats" \
                   " (*.mat *.qml *.xml *.kor *.evt)"
            caption = "Open an event file"
            fname = QFileDialog().getOpenFileName(self, caption=caption,
                                                  filter=filt,
                                                  dir=self.get_current_event_path())
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

    def saveData(self, directory=None, outformat=None):

        def getSavePath(e, directory, outformat):
            print('warning: {0}'.format(e))
            if not directory:
                dlgflag = 1
                directory = self.get_current_event_path()
            else:
                dlgflag = 0
            eventname = self.get_current_event_name()
            filename = 'PyLoT_' + eventname
            outpath = os.path.join(directory, filename)
            title = 'Save pick data ...'
            file_filter = "(*.xml *.obs *.cnv)"

            if dlgflag == 1:
                fname, selected_filter = QFileDialog.getSaveFileName(self,
                                                                     title,
                                                                     outpath,
                                                                     file_filter)

            try:
                fbasename = fname
            except:
                fbasename = outpath
            exform = ['.obs', '.xml', '.cnv']
            return fbasename, exform

        settings = QSettings()
        fbasename = self.getEventFileName()
        exform = settings.value('data/exportFormat', 'QUAKEML')

        uppererrorP = self._inputs['timeerrorsP']
        uppererrorS = self._inputs['timeerrorsS']

        try:
            self.get_data().applyEVTData(self.get_current_event(), typ='event')  # getPicks())
        except OverwriteError:
            #     msgBox = QMessageBox()
            #     msgBox.setText("Picks have been modified!")
            #     msgBox.setInformativeText(
            #         "Do you want to save the changes and overwrite the picks?")
            #     msgBox.setDetailedText(self.get_data().getPicksStr())
            #     msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
            #     msgBox.setDefaultButton(QMessageBox.Save)
            #     ret = msgBox.exec_()
            #     if ret == QMessageBox.Save:
            self.get_data().resetPicks()
            return self.saveData()
        # elif ret == QMessageBox.Cancel:
        #         return False
        # MP MP changed to suppress unnecessary user prompt
        try:
            self.get_data().exportEvent(fbasename, exform, upperErrors=[uppererrorP[3], uppererrorS[3]])
        except FormatError as e:
            fbasename, exform = getSavePath(e, directory, outformat)
        except AttributeError as e:
            fbasename, exform = getSavePath(e, directory, outformat)

        # catch all possible cases before going on
        if not fbasename:
            return False
        # warn overwriting
        # elif os.path.exists(fbasename + exform):
        #     ans = QMessageBox.question(self, self.tr("Overwrite file..."),
        #                        self.tr("File already exists: {0}\n".format(fbasename + exform) + \
        #                           "Overwrite file anyway?"),
        #                        QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        #     # only negative answers have to be caught
        #     if ans == QMessageBox.No:
        #         self.saveData()
        #     elif ans == QMessageBox.Cancel:
        #         return False

        # export to given path
        self.get_data().exportEvent(fbasename, exform[0], upperErrors=[uppererrorP[3], uppererrorS[3]])
        self.get_data().exportEvent(fbasename, exform[1], upperErrors=[uppererrorP[3], uppererrorS[3]])
        # self.get_data().exportEvent(fbasename, exform[2], upperErrors=[uppererrorP[3], uppererrorS[3]])
        # all files save (ui clean)
        self.update_status('Picks saved as %s, %s, and %s' % (fbasename + exform[0], fbasename + exform[1],
                                                              fbasename + exform[2]))
        self.disableSaveManualPicksAction()
        return True

    def enableSaveManualPicksAction(self):
        self.saveManualPicksAction.setEnabled(True)

    def disableSaveManualPicksAction(self):
        self.saveManualPicksAction.setEnabled(False)

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
        if type == 'manual':
            return self.get_current_event().getPicks()
        if type == 'auto':
            return self.get_current_event().getAutopicks()
            # rdict = dict(auto=self.pylot_autopicks, manual=self.pylot_picks)
            # return rdict[type]

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
    def getWFID(ycoord):
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
            qmb = QMessageBox(self, icon=QMessageBox.Question,
                              text='Do you wish to save changes in the current project?')
            qmb.setStandardButtons(QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            qmb.setDefaultButton(QMessageBox.Save)
            ret = qmb.exec_()
            if ret == qmb.Save:
                return self.saveProject()
            elif ret == qmb.Cancel:
                return False
        return True

    def enableRefTestButtons(self, bool):
        '''
        Enable/disable reference and test event selection buttons.
        '''
        self.ref_event_button.setEnabled(bool)
        self.test_event_button.setEnabled(bool)

    def refreshRefTestButtons(self):
        '''
        Refresh state of reference and test event selection buttons depending on current event.
        '''
        event = self.get_current_event()
        if event:
            self.ref_event_button.setChecked(event.isRefEvent())
            self.test_event_button.setChecked(event.isTestEvent())
            self.enableRefTestButtons(bool(self.get_current_event().pylot_picks))
            return
        self.ref_event_button.setChecked(False)
        self.test_event_button.setChecked(False)
        self.enableRefTestButtons(False)

    def toggleRef(self):
        '''
        Toggle ref/test buttons when reference button gets clicked.
        '''
        ref = self.ref_event_button.isChecked()
        self.test_event_button.setChecked(False)
        self.get_current_event().setTestEvent(False)
        self.get_current_event().setRefEvent(ref)
        self.fill_eventbox()
        self.refreshTabs()
        if self.tap:
            self.tap.fill_eventbox()

    def toggleTest(self):
        '''
        Toggle ref/test buttons when test button gets clicked.
        '''
        test = self.test_event_button.isChecked()
        self.ref_event_button.setChecked(False)
        self.get_current_event().setRefEvent(False)
        self.get_current_event().setTestEvent(test)
        self.fill_eventbox()
        self.refreshTabs()
        if self.tap:
            self.tap.fill_eventbox()

    def refreshEvents(self):
        '''
        Refresh GUI when events get changed.
        '''
        # Attribute _eventChanged refers to change in first _eventChanged[0]
        # and second _eventChanged[1] main tabs.
        # E.g. plotting is not necessary when changing event in second tab and
        # array_map refresh is not necessary when changing event in waveform plot tab,
        # but gets necessary when switching from one to another after changing an event.
        self._eventChanged = [True, True]
        self.refreshTabs()

    def refreshTabs(self):
        '''
        Refresh main tabs depending on change of the current event in first/second tab.
        '''
        # Logical problem appearing somewhere when calling this function:
        # Loading an existing project from array_tab leads to two calls of newWF
        # which will read in data input twice. Therefore current tab is changed to 0
        # in loadProject before calling this function.
        plotted = False
        # only refresh first/second tab when an event was changed.
        if self._eventChanged[0] or self._eventChanged[1]:
            event = self.get_current_event()
            if not event:
                return
            # update picks saved in GUI mainwindow (to be changed in future!!) MP MP
            if not event.pylot_picks:
                self.pylot_picks = {}
            else:
                self.pylot_picks = event.pylot_picks
            if not event.pylot_autopicks:
                self.pylot_autopicks = {}
            else:
                self.pylot_autopicks = event.pylot_autopicks
        # if current tab is waveformPlot-tab and the data in this tab was not yet refreshed
        if self.tabs.currentIndex() == 0:
            if self._eventChanged[0]:
                if hasattr(self.project, 'eventlist'):
                    if len(self.project.eventlist) > 0:
                        self.newWF()
                        # keep track whether event was already plotted
                        plotted = True
        # if current tab is array_map-tab and the data in this tab was not yet refreshed
        if self.tabs.currentIndex() == 1:
            if self._eventChanged[1]:
                self.refresh_array_map()
                if not plotted and self._eventChanged[0]:
                    # newWF(False) = load data without plotting
                    self.newWF(plot=False)

        if self.tabs.currentIndex() == 2:
            self.init_event_table()
        self.refreshRefTestButtons()

    def newWF(self, plot=True):
        '''
        Load new data and plot if necessary.
        '''
        self.loadWaveformDataThread(plot=plot)
        if plot:
            self._eventChanged[0] = False

    def loadWaveformDataThread(self, plot=True, load=True):
        '''
        Generates a modal thread to load waveform data and
        call modal plot thread method when finished.
        '''
        if load:
            self.wfd_thread = Thread(self, self.loadWaveformData,
                                     progressText='Reading data input...')
        if load and plot:
            self.wfd_thread.finished.connect(self.plotWaveformDataThread)

        if load:
            self.wfd_thread.start()

        if plot and not load:
            self.plotWaveformDataThread()

    def loadWaveformData(self):
        '''
        Load waveform data corresponding to current selected event.
        '''
        # if self.fnames and self.okToContinue():
        #     self.setDirty(True)
        #     ans = self.data.setWFData(self.fnames)
        # elif self.fnames is None and self.okToContinue():
        #     ans = self.data.setWFData(self.getWFFnames())
        # else:
        #     ans = False
        self.fnames = self.getWFFnames_from_eventbox()
        self.data.setWFData(self.fnames)
        self._stime = full_range(self.get_data().getWFData())[0]

    def connectWFplotEvents(self):
        '''
        Connect signals refering to WF-Dataplot (select station, tutor_user, scrolling)
        '''
        if self.pg:
            self.connect_pg()
        else:
            self.connect_mpl()

    def connect_pg(self):
        self.poS_id = self.dataPlot.plotWidget.scene().sigMouseClicked.connect(self.pickOnStation)

    def connect_mpl(self):
        if not self.poS_id:
            self.poS_id = self.dataPlot.mpl_connect('button_press_event',
                                                    self.pickOnStation)
        if not self.ae_id:
            self.ae_id = self.dataPlot.mpl_connect('axes_enter_event',
                                                   lambda event: self.tutor_user())

        if not self.scroll_id:
            self.scroll_id = self.dataPlot.mpl_connect('scroll_event',
                                                       self.scrollPlot)

    def disconnectWFplotEvents(self):
        '''
        Disconnect all signals refering to WF-Dataplot (select station, tutor_user, scrolling)
        '''
        if self.pg:
            self.disconnect_pg()
        else:
            self.disconnect_mpl()

    def disconnect_pg(self):
        if self.poS_id:
            try:
                self.dataPlot.plotWidget.scene().sigMouseClicked.disconnect()
            except:
                pass

    def disconnect_mpl(self):
        if self.poS_id:
            self.dataPlot.mpl_disconnect(self.poS_id)
        if self.ae_id:
            self.dataPlot.mpl_disconnect(self.ae_id)
        if self.scroll_id:
            self.dataPlot.mpl_disconnect(self.scroll_id)
        self.poS_id = None
        self.ae_id = None
        self.scroll_id = None

    def finish_pg_plot(self):
        self.getPlotWidget().updateWidget()
        plots = self.wfp_thread.data
        for times, data in plots:
            self.dataPlot.plotWidget.getPlotItem().plot(times, data, pen='k')
        self.dataPlot.reinitMoveProxy()
        self.dataPlot.plotWidget.showAxis('left')
        self.dataPlot.plotWidget.showAxis('bottom')

    def finishWaveformDataPlot(self):
        if self.pg:
            self.finish_pg_plot()
        else:
            self._max_xlims = self.dataPlot.getXLims()
        plotWidget = self.getPlotWidget()
        plotDict = plotWidget.getPlotDict()
        pos = plotDict.keys()
        labels = [plotDict[n][2] + '.' + plotDict[n][0] for n in pos]
        plotWidget.setYTickLabels(pos, labels)
        try:
            plotWidget.figure.tight_layout()
        except:
            pass
        self.connectWFplotEvents()
        self.loadlocationaction.setEnabled(True)
        self.auto_tune.setEnabled(True)
        self.auto_pick.setEnabled(True)
        self.z_action.setEnabled(True)
        self.e_action.setEnabled(True)
        self.n_action.setEnabled(True)
        self.openmanualpicksaction.setEnabled(True)
        self.openautopicksaction.setEnabled(True)
        self.loadpilotevent.setEnabled(True)
        event = self.get_current_event()
        if event.pylot_picks:
            self.pylot_picks = event.pylot_picks
            self.drawPicks(picktype='manual')
            self.enableSaveManualPicksAction()
            self.locateEvent.setEnabled(True)
        if event.pylot_autopicks:
            self.pylot_autopicks = event.pylot_autopicks
            self.drawPicks(picktype='auto')
            self.compare_action.setEnabled(True)
        self.draw()

    def clearWaveformDataPlot(self):
        self.disconnectWFplotEvents()
        if self.pg:
            self.dataPlot.plotWidget.getPlotItem().clear()
            self.dataPlot.plotWidget.hideAxis('bottom')
            self.dataPlot.plotWidget.hideAxis('left')
        else:
            self.dataPlot.getAxes().cla()
        self.loadlocationaction.setEnabled(False)
        self.auto_tune.setEnabled(False)
        self.auto_pick.setEnabled(False)
        self.z_action.setEnabled(False)
        self.e_action.setEnabled(False)
        self.n_action.setEnabled(False)
        self.openmanualpicksaction.setEnabled(False)
        self.openautopicksaction.setEnabled(False)
        self.loadpilotevent.setEnabled(False)
        self.disableSaveManualPicksAction()
        self.draw()

    def plotWaveformDataThread(self):
        '''
        Open a modal thread to plot current waveform data.
        '''
        self.clearWaveformDataPlot()
        self.wfp_thread = Thread(self, self.plotWaveformData,
                                 progressText='Plotting waveform data...')
        self.wfp_thread.finished.connect(self.finishWaveformDataPlot)
        self.wfp_thread.start()

    def plotWaveformData(self):
        '''
        Plot waveform data to current plotWidget.
        '''
        settings = QSettings()
        nth_sample = settings.value('nth_sample')
        if not nth_sample:
            nth_sample = 1
        zne_text = {'Z': 'vertical', 'N': 'north-south', 'E': 'east-west'}
        comp = self.getComponent()
        title = 'section: {0} components'.format(zne_text[comp])
        wfst = self.get_data().getWFData()
        # wfst = self.get_data().getWFData().select(component=comp)
        # wfst += self.get_data().getWFData().select(component=alter_comp)
        plotWidget = self.getPlotWidget()
        self.adjustPlotHeight()
        plots = plotWidget.plotWFData(wfdata=wfst, title=title, mapping=False, component=comp,
                                      nth_sample=int(nth_sample))
        return plots

    def adjustPlotHeight(self):
        if self.pg:
            return
        height_need = len(self.data.getWFData()) * self.height_factor
        plotWidget = self.getPlotWidget()
        if self.tabs.widget(0).frameSize().height() < height_need:
            plotWidget.setMinimumHeight(height_need)
        else:
            plotWidget.setMinimumHeight(0)

    def plotZ(self):
        self.setComponent('Z')
        self.plotWaveformDataThread()
        # self.drawPicks()
        # self.draw()

    def plotN(self):
        self.setComponent('N')
        self.plotWaveformDataThread()
        # self.drawPicks()
        # self.draw()

    def plotE(self):
        self.setComponent('E')
        self.plotWaveformDataThread()
        # self.drawPicks()
        # self.draw()

    def pushFilterWF(self, param_args):
        self.get_data().filterWFData(param_args)

    def filterWaveformData(self):
        if self.get_data():
            if self.getFilterOptions() and self.filterAction.isChecked():
                kwargs = self.getFilterOptions()[self.getSeismicPhase()].parseFilterOptions()
                self.pushFilterWF(kwargs)
            elif self.filterAction.isChecked():
                self.adjustFilterOptions()
            else:
                self.get_data().resetWFData()
        self.plotWaveformDataThread()
        self.drawPicks()
        self.draw()

    def adjustFilterOptions(self):
        fstring = "Filter Options"
        self.filterDlg = FilterOptionsDialog(titleString=fstring,
                                             parent=self)
        if self.filterDlg.exec_():
            filteroptions = self.filterDlg.getFilterOptions()
            self.setFilterOptions(filteroptions)
            if self.filterAction.isChecked():
                kwargs = self.getFilterOptions()[self.getSeismicPhase()].parseFilterOptions()
                self.pushFilterWF(kwargs)
                self.plotWaveformDataThread()

    def checkFilterOptions(self):
        fstring = "Filter Options"
        self.filterDlg = FilterOptionsDialog(titleString=fstring,
                                             parent=self)
        filteroptions = self.filterDlg.getFilterOptions()
        self.setFilterOptions(filteroptions)
        filterP = filteroptions['P']
        filterS = filteroptions['S']
        minP, maxP = filterP.getFreq()
        minS, maxS = filterS.getFreq()
        self.paraBox.params_to_gui()

    def getFilterOptions(self):
        return self.filteroptions
        # try:
        #     return self.filteroptions[self.getSeismicPhase()]
        # except AttributeError as e:
        #     print(e)
        #     return FilterOptions(None, None, None)

    def getFilters(self):
        return self.filteroptions

    def setFilterOptions(self, filterOptions):  # , seismicPhase=None):
        # if seismicPhase is None:
        #     self.getFilterOptions()[self.getSeismicPhase()] = filterOptions
        # else:
        #     self.getFilterOptions()[seismicPhase] = filterOptions
        self.filterOptions = filterOptions
        filterP = filterOptions['P']
        filterS = filterOptions['S']
        minP, maxP = filterP.getFreq()
        minS, maxS = filterS.getFreq()
        self._inputs.setParamKV('minfreq', (minP, minS))
        self._inputs.setParamKV('maxfreq', (maxP, maxS))
        self._inputs.setParamKV('filter_order', (filterP.getOrder(), filterS.getOrder()))
        self._inputs.setParamKV('filter_type', (filterP.getFilterType(), filterS.getFilterType()))

    def filterOptionsFromParameter(self):
        minP, minS = self._inputs['minfreq']
        maxP, maxS = self._inputs['maxfreq']
        orderP, orderS = self._inputs['filter_order']
        typeP, typeS = self._inputs['filter_type']

        filterP = self.getFilterOptions()['P']
        filterP.setFreq([minP, maxP])
        filterP.setOrder(orderP)
        filterP.setFilterType(typeP)

        filterS = self.getFilterOptions()['S']
        filterS.setFreq([minS, maxS])
        filterS.setOrder(orderS)
        filterS.setFilterType(typeS)

        self.checkFilterOptions()

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
        plot_dict = self.getPlotWidget().getPlotDict()
        if wfID in plot_dict.keys():
            return plot_dict[wfID][0]

    def getNetworkName(self, wfID):
        plot_dict = self.getPlotWidget().getPlotDict()
        if wfID in plot_dict.keys():
            return plot_dict[wfID][2]

    def alterPhase(self):
        pass

    def setSeismicPhase(self, phase):
        self.seismicPhase = self.seismicPhaseButtonGroup.getValue()
        self.update_status('Seismic phase changed to '
                           '{0}'.format(self.getSeismicPhase()))

    def scrollPlot(self, gui_event):
        '''
        Function connected to mouse wheel scrolling inside WFdataPlot.
        Only used if amount of traces exceeds a certain limit creating
        a scroll area.
        '''
        button = gui_event.button
        x, y = gui_event.xdata, gui_event.ydata
        if not button == 'up' and not button == 'down':
            return
        if not self._ctrl and not self._shift:
            vbar = self.wf_scroll_area.verticalScrollBar()
            up_down = {'up': -60,
                       'down': 60}
            if vbar.maximum():
                vbar.setValue(vbar.value() + up_down[button])
        if self._ctrl:
            factor = {'up': 5. / 4.,
                      'down': 4. / 5.}
            self.height_factor *= factor[button]
            self.adjustPlotHeight()
        if self._shift:
            factor = {'up': 1. / 2.,
                      'down': 2.}
            xlims = self.dataPlot.getXLims()
            xdiff = xlims[1] - xlims[0]
            xdiff *= factor[button]
            xl = x - 0.5 * xdiff
            xr = x + 0.5 * xdiff
            if xl < self._max_xlims[0]:
                xl = self._max_xlims[0]
            if xr > self._max_xlims[1]:
                xr = self._max_xlims[1]
            self.dataPlot.setXLims((xl, xr))
            self.dataPlot.draw()

    def pickOnStation(self, gui_event):
        if self.pg:
            if not gui_event.button() == 1:
                return
        else:
            if not gui_event.button == 1:
                return

        if self.pg:
            ycoord = self.dataPlot.plotWidget.getPlotItem().vb.mapSceneToView(gui_event.scenePos()).y()
        else:
            ycoord = gui_event.ydata

        wfID = self.getWFID(ycoord)

        if wfID is None: return
        self.pickDialog(wfID)

    def pickDialog(self, wfID, nextStation=False):
        station = self.getStationName(wfID)
        network = self.getNetworkName(wfID)
        if not station:
            return
        self.update_status('picking on station {0}'.format(station))
        data = self.get_data().getWFData()
        pickDlg = PickDlg(self, parameter=self._inputs,
                          data=data.select(station=station),
                          station=station, network=network,
                          picks=self.getPicksOnStation(station, 'manual'),
                          autopicks=self.getPicksOnStation(station, 'auto'))
        pickDlg.nextStation.setChecked(nextStation)
        if pickDlg.exec_():
            if pickDlg.getPicks():
                self.setDirty(True)
                self.update_status('picks accepted ({0})'.format(station))
                replot = self.addPicks(station, pickDlg.getPicks())
                self.get_current_event().setPick(station, pickDlg.getPicks())
                self.enableSaveManualPicksAction()
                if replot:
                    self.plotWaveformDataThread()
                    self.drawPicks()
                    self.draw()
                else:
                    self.drawPicks(station)
                    self.draw()
            if pickDlg.nextStation.isChecked():
                self.pickDialog(wfID - 1, nextStation=pickDlg.nextStation.isChecked())
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

    def tune_autopicker(self):
        '''
        Initiates TuneAutopicker widget use to interactively
        tune parameters for autopicking algorithm.
        '''
        # figures and canvas have to be iniated from the main GUI
        # thread to prevent handling of QPixmap objects outside of
        # the main thread
        self.fig_dict = {}
        self.canvas_dict = {}
        self.fig_keys = [
            'mainFig',
            'aicFig',
            'slength',
            'checkZ4s',
            'refPpick',
            'el_Ppick',
            'fm_picker',
            'el_S1pick',
            'el_S2pick',
            'refSpick',
            'aicARHfig',
        ]
        for key in self.fig_keys:
            fig = Figure()
            self.fig_dict[key] = fig

        if not self.tap:
            # init TuneAutopicker object
            self.tap = TuneAutopicker(self)
            # first call of update to init tabs with empty canvas
            self.update_autopicker()
            # connect update signal of TuneAutopicker with update function
            # creating and filling figure canvas
            self.tap.update.connect(self.update_autopicker)
            self.tap.figure_tabs.setCurrentIndex(0)
        else:
            self.update_autopicker()
            self.tap.fill_eventbox()
        self.tap.show()

    def update_autopicker(self):
        '''
        Create and fill TuneAutopicker tabs with figure canvas.
        '''
        for key in self.fig_keys:
            self.canvas_dict[key] = FigureCanvas(self.fig_dict[key])
        self.tap.fill_tabs(picked=True)

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
        receventid = self.get_current_event_path()
        self.thread = AutoPickThread(parent=self,
                                     func=autoPyLoT,
                                     infile=self.infile,
                                     fnames=self.fnames,
                                     eventid=receventid,
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
        if not stat_picks:
            rval = False
        else:
            # set picks (ugly syntax?)
            self.getPicks(type=type)[station] = picks
            rval = True
        return rval
        # if not stat_picks:
        #     stat_picks = picks
        # else:
        #     msgBox = QMessageBox(self)
        #     msgBox.setText("The picks for station {0} have been "
        #                    "changed.".format(station))
        #     msgBox.setDetailedText("Old picks:\n"
        #                            "{old_picks}\n\n"
        #                            "New picks:\n"
        #                            "{new_picks}".format(old_picks=stat_picks,
        #                                                 new_picks=picks))
        #     msgBox.setInformativeText("Do you want to save your changes?")
        #     msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
        #     msgBox.setDefaultButton(QMessageBox.Save)
        #     ret = msgBox.exec_()
        #     if ret == QMessageBox.Save:
        #         stat_picks = picks
        #         rval = True
        #     elif ret == QMessageBox.Cancel:
        #         pass
        #     else:
        #         raise Exception('FATAL: Should never occur!')
        # MP MP prompt redundant because new picks have to be accepted in the first place closing PickDlg

    def updatePicks(self, type='manual', event=None):
        if not event:
            event = self.get_current_event()
        picks = picksdict_from_picks(evt=self.get_data(type).get_evt_data())
        if type == 'manual':
            event.addPicks(picks)
            self.pylot_picks.update(picks)
        elif type == 'auto':
            event.addAutopicks(picks)
            self.pylot_autopicks.update(picks)
        self.check4Comparison()

    def drawPicks(self, station=None, picktype='manual'):
        # if picks to draw not specified, draw all picks available
        if not station:
            for station in self.getPicks(type=picktype):
                self.drawPicks(station, picktype=picktype)
            return
        # plotting picks
        plotID = self.getStationID(station)
        if plotID is None:
            return
        if self.pg:
            pw = self.getPlotWidget().plotWidget
        else:
            ax = self.getPlotWidget().axes
        ylims = np.array([-.5, +.5]) + plotID
        if self.pg:
            dashed = QtCore.Qt.DashLine
            dotted = QtCore.Qt.DotLine
            phase_col = {
                'P': (pg.mkPen('c'), pg.mkPen((0, 255, 255, 100), style=dashed), pg.mkPen('b', style=dashed),
                      pg.mkPen('b', style=dotted)),
                'S': (pg.mkPen('m'), pg.mkPen((255, 0, 255, 100), style=dashed), pg.mkPen('r', style=dashed),
                      pg.mkPen('r', style=dotted))
            }
        else:
            phase_col = {
                'P': ('c', 'c--', 'b-', 'bv', 'b^', 'b'),
                'S': ('m', 'm--', 'r-', 'rv', 'r^', 'r')
            }

        stat_picks = self.getPicks(type=picktype)[station]

        stime = self.getStime()

        for phase in stat_picks:
            picks = stat_picks[phase]
            if type(stat_picks[phase]) is not dict and type(stat_picks[phase]) is not AttribDict:
                return
            colors = phase_col[phase[0].upper()]

            mpp = picks['mpp'] - stime
            if picks['epp'] and picks['lpp']:
                epp = picks['epp'] - stime
                lpp = picks['lpp'] - stime
            else:
                epp = None
                lpp = None
            spe = picks['spe']

            if not spe and epp and lpp:
                spe = symmetrize_error(mpp - epp, lpp - mpp)

            if self.pg:
                if picktype == 'manual':
                    if picks['epp'] and picks['lpp']:
                        pw.plot([epp, epp], ylims,
                                alpha=.25, pen=colors[0], name='EPP')
                        pw.plot([lpp, lpp], ylims,
                                alpha=.25, pen=colors[0], name='LPP')
                    if spe:
                        spe_l = pg.PlotDataItem([mpp - spe, mpp - spe], ylims, pen=colors[1],
                                                name='{}-SPE'.format(phase))
                        spe_r = pg.PlotDataItem([mpp + spe, mpp + spe], ylims, pen=colors[1])
                        pw.addItem(spe_l)
                        pw.addItem(spe_r)
                        try:
                            fill = pg.FillBetweenItem(spe_l, spe_r, brush=colors[1].brush())
                            fb = pw.addItem(fill)
                        except:
                            print('Warning: drawPicks: Could not create fill for symmetric pick error.')
                        pw.plot([mpp, mpp], ylims, pen=colors[2], name='{}-Pick'.format(phase))
                    else:
                        pw.plot([mpp, mpp], ylims, pen=colors[0], name='{}-Pick (NO PICKERROR)'.format(phase))
                elif picktype == 'auto':
                    pw.plot([mpp, mpp], ylims, pen=colors[3])
                else:
                    raise TypeError('Unknown picktype {0}'.format(picktype))
            else:
                if picktype == 'manual':
                    if picks['epp'] and picks['lpp']:
                        ax.fill_between([epp, lpp], ylims[0], ylims[1],
                                        alpha=.25, color=colors[0], label='EPP, LPP')
                    if spe:
                        ax.plot([mpp - spe, mpp - spe], ylims, colors[1], label='{}-SPE'.format(phase))
                        ax.plot([mpp + spe, mpp + spe], ylims, colors[1])
                        ax.plot([mpp, mpp], ylims, colors[2], label='{}-Pick'.format(phase))
                    else:
                        ax.plot([mpp, mpp], ylims, colors[6], label='{}-Pick (NO PICKERROR)'.format(phase))
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
        obsdir = os.path.join(locroot, 'obs')
        self.saveData(directory=obsdir, outformat='.obs')
        eventname = self.get_current_event_name()
        filename = 'PyLoT_' + eventname
        locpath = os.path.join(locroot, 'loc', filename)
        phasefile = os.path.join(obsdir, filename + '.obs')
        phasepath = os.path.join(locroot, 'obs', phasefile)
        lt.modify_inputs(infile, locroot, filename, phasefile, ttt)
        try:
            lt.locate(infile)
        except RuntimeError as e:
            print(e.message)
        finally:
            os.remove(phasepath)

        self.get_data().applyEVTData(lt.read_location(locpath), typ='event')
        # self.get_data().applyEVTData(self.calc_magnitude(), typ='event')

    def init_array_tab(self):
        '''
        Init second main tab if neither metadata nor
        array map are given. A button will be show calling self.get_metadata.
        '''
        if hasattr(self, 'metadata_widget'):
            if self.metadata_widget:
                self.metadata_widget.setParent(None)
                self.array_layout.removeWidget(self.metadata_widget)
        if hasattr(self, 'array_map'):
            if self.array_map:
                self.array_map.setParent(None)
                self.array_layout.removeWidget(self.array_map)
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

        self.metadata = None
        self.metadata_widget.setLayout(grid_layout)
        self.array_layout.addWidget(self.metadata_widget)

    def init_array_map(self, index=1):
        '''
        Try to init array map widget. If no metadata are given,
        self.get_metadata will be called.
        '''
        self.tabs.setCurrentIndex(1)
        if hasattr(self, 'metadata_widget'):
            if self.metadata_widget:
                self.metadata_widget.setParent(None)
                self.array_layout.removeWidget(self.metadata_widget)
        if hasattr(self, 'array_map'):
            if self.array_map:
                self.array_map.setParent(None)
                self.array_layout.removeWidget(self.array_map)
        if not self.array_map:
            self.get_metadata()
            if not self.metadata:
                return
        self.am_figure = Figure()
        self.am_canvas = FigureCanvas(self.am_figure)
        self.am_toolbar = NavigationToolbar(self.am_canvas, self)
        self.array_map = map_projection(self)
        # self.array_map_thread()
        self.array_layout.addWidget(self.array_map)
        self.tabs.setCurrentIndex(index)
        self.refresh_array_map()

    def array_map_thread(self):
        '''
        Start modal thread to init the array_map object.
        '''
        # Note: basemap generation freezes GUI but cannot be threaded as it generates a Pixmap.
        self.amt = Thread(self, self.array_map.init_map, arg=None, progressText='Generating map...')
        self.amt.finished.connect(self.finish_array_map)
        self.amt.start()

    def finish_array_map(self):
        '''
        Add array_map to GUI tab when array_map_thread has finished.
        '''
        self.array_map = self.amt.data
        self.array_layout.addWidget(self.array_map)
        # self.tabs.setCurrentIndex(index)
        # self.refresh_array_map()

    def refresh_array_map(self):
        '''
        Refresh array map when current event is changed.
        '''
        if not self.array_map:
            return
        # refresh with new picks here!!!
        event = self.get_current_event()
        if hasattr(event, 'origins'):
            if event.origins:
                lat = event.origins[0].latitude
                lon = event.origins[0].longitude
                self.array_map.eventLoc = (lat, lon)
        self.array_map.refresh_drawings(self.get_current_event().getPicks())
        self._eventChanged[1] = False

    def init_event_table(self, tabindex=2):
        '''
        Build and initiate event table (3rd tab [index=2]) containing information of every event.
        '''

        def set_enabled(item, enabled=True, checkable=False):
            # modify item flags depending on case needed
            if enabled and not checkable:
                item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable)
            elif enabled and checkable:
                item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsSelectable)
            else:
                item.setFlags(QtCore.Qt.ItemIsSelectable)

        if self.project:
            eventlist = self.project.eventlist
        else:
            eventlist = []

        def cell_changed(row=None, column=None):
            # connected to cell changes in event table
            # changes attributes of the corresponding event
            table = self.project._table
            event = self.project.getEventFromPath(table[row][0].text())
            if column == 8 or column == 9:
                # toggle checked states (exclusive)
                item_ref = table[row][8]
                item_test = table[row][9]
                if column == 8 and item_ref.checkState():
                    item_test.setCheckState(QtCore.Qt.Unchecked)
                    event.setRefEvent(True)
                elif column == 8 and not item_ref.checkState():
                    event.setRefEvent(False)
                elif column == 9 and item_test.checkState():
                    item_ref.setCheckState(QtCore.Qt.Unchecked)
                    event.setTestEvent(True)
                elif column == 9 and not item_test.checkState():
                    event.setTestEvent(False)
                self.fill_eventbox()
            elif column == 5:
                # update event notes
                notes = table[row][5].text()
                event.addNotes(notes)
                self.fill_eventbox()

        # remove old table
        if hasattr(self, 'event_table'):
            self.event_table.setParent(None)
            self.events_layout.removeWidget(self.event_table)

        # init new qtable
        self.event_table = QtGui.QTableWidget()
        self.event_table.setColumnCount(11)
        self.event_table.setRowCount(len(eventlist))
        self.event_table.setHorizontalHeaderLabels(['Event',
                                                    'Time',
                                                    'Lat',
                                                    'Lon',
                                                    'Depth',
                                                    'Mag',
                                                    '[N] MP',
                                                    '[N] AP',
                                                    'Tuning Set',
                                                    'Test Set',
                                                    'Notes'])

        # iterate through eventlist and generate items for table rows
        self.project._table = []
        for index, event in enumerate(eventlist):
            event_npicks = 0
            event_nautopicks = 0
            if event.pylot_picks:
                event_npicks = len(event.pylot_picks)
            if event.pylot_autopicks:
                event_nautopicks = len(event.pylot_autopicks)
            item_path = QtGui.QTableWidgetItem()
            item_time = QtGui.QTableWidgetItem()
            item_lat = QtGui.QTableWidgetItem()
            item_lon = QtGui.QTableWidgetItem()
            item_depth = QtGui.QTableWidgetItem()
            item_mag = QtGui.QTableWidgetItem()
            item_nmp = QtGui.QTableWidgetItem(str(event_npicks))
            item_nmp.setIcon(self.manupicksicon_small)
            item_nap = QtGui.QTableWidgetItem(str(event_nautopicks))
            item_nap.setIcon(self.autopicksicon_small)
            item_ref = QtGui.QTableWidgetItem()
            item_test = QtGui.QTableWidgetItem()
            item_notes = QtGui.QTableWidgetItem()

            item_ref.setBackground(self._colors['ref'])
            item_test.setBackground(self._colors['test'])
            item_path.setText(event.path)
            if hasattr(event, 'origins'):
                if event.origins:
                    origin = event.origins[0]
                    item_time.setText(str(origin.time).split('.')[0])
                    item_lon.setText(str(origin.longitude))
                    item_lat.setText(str(origin.latitude))
                    item_depth.setText(str(origin.depth))
            if hasattr(event, 'magnitudes'):
                if event.magnitudes:
                    magnitude = event.magnitudes[0]
                    item_mag.setText(str(magnitude.mag))
            item_notes.setText(event.notes)

            set_enabled(item_path, True, False)
            set_enabled(item_nmp, True, False)
            set_enabled(item_nap, True, False)
            if event.pylot_picks:
                set_enabled(item_ref, True, True)
                set_enabled(item_test, True, True)
            else:
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

            column = [item_path, item_time, item_lat, item_lon, item_depth, item_mag,
                      item_nmp, item_nap, item_ref, item_test, item_notes]
            self.project._table.append(column)

        for r_index, row in enumerate(self.project._table):
            for c_index, item in enumerate(row):
                self.event_table.setItem(r_index, c_index, item)

        header = self.event_table.horizontalHeader()
        header.setResizeMode(QtGui.QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)
        self.event_table.cellChanged[int, int].connect(cell_changed)

        self.events_layout.addWidget(self.event_table)
        self.tabs.setCurrentIndex(tabindex)

    def read_metadata_thread(self, fninv):
        self.rm_thread = Thread(self, read_metadata, arg=fninv, progressText='Reading metadata...')
        self.rm_thread.finished.connect(self.set_metadata)
        self.rm_thread.start()

    def set_metadata(self):
        settings = QSettings()
        self.metadata = self.rm_thread.data
        if settings.value('saveMetadata'):
            self.project.metadata = self.rm_thread.data
        self.project.inv_path = settings.value("inventoryFile")
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

        settings = QSettings()

        if hasattr(self.project, 'metadata'):
            self.metadata = self.project.metadata
            return True
        if hasattr(self.project, 'invPath'):
            settings.setValue("inventoryFile", self.project.inv_path)

        fninv = settings.value("inventoryFile", None)

        if fninv is None and not self.metadata:
            if not set_inv(settings):
                return None
        elif fninv is not None and not self.metadata:
            if not hasattr(self.project, 'invPath'):
                ans = QMessageBox.question(self, self.tr("Use default metadata..."),
                                           self.tr(
                                               "Do you want to use the default value for metadata?"),
                                           QMessageBox.Yes | QMessageBox.No,
                                           QMessageBox.Yes)
                if ans == QMessageBox.No:
                    if not set_inv(settings):
                        return None
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
            local_mag = LocalMagnitude(corr_wf, self.get_data().get_evt_data(), self.inputs.get('sstop'),
                                       verbosity=True)
            return local_mag.updated_event()
        elif type == 'Mw':
            moment_mag = MomentMagnitude(corr_wf, self.get_data().get_evt_data(), self.inputs.get('vp'),
                                         self.inputs.get('Qp'), self.inputs.get('rho'), verbosity=True)
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
            if not self.get_current_event() or not self.project.location:
                self.setWindowTitle("PyLoT - New project [*]")
            elif self.get_current_event():
                self.setWindowTitle("PyLoT - {} [*]".format(self.project.location))
            else:
                self.setWindowTitle(
                    "PyLoT - seismic processing the python way[*]")
        self.setWindowModified(self.dirty)

    def tutor_user(self):
        trace_pick = ' select trace to pick on station ...'
        strg_key = ' - [CTRL + mousewheel] vertical spacing'
        shift_key = ' - [SHIFT + mousewheel] horizontal zoom'
        message = trace_pick + strg_key + shift_key
        self.update_status(message, 10000)

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

    def createNewProject(self):
        '''
        Create new project file.
        '''
        if not self.okToContinue():
            return
        self.project = Project()
        self.init_events(new=True)
        self.setDirty(False)
        self.project.parameter = self._inputs
        self.saveProjectAsAction.setEnabled(True)
        self.update_status('Created new project...', duration=1000)
        return True

    def loadProject(self, fnm=None):
        '''
        Load an existing project file.
        '''
        if not self.okToContinue():
            return
        if not fnm:
            dlg = QFileDialog()
            fnm = dlg.getOpenFileName(self, 'Open project file...', filter='Pylot project (*.plp)')
            if not fnm:
                return
            fnm = fnm[0]
        if fnm:
            self.project = Project.load(fnm)
            if hasattr(self.project, 'parameter'):
                if self.project.parameter:
                    self._inputs = self.project.parameter
            self.tabs.setCurrentIndex(0)  # implemented to prevent double-loading of waveform data
            self.init_events(new=True)
            self.setDirty(False)
            if hasattr(self.project, 'metadata'):
                if self.project.metadata:
                    self.init_array_map(index=0)
                    return
            self.init_array_tab()

    def saveProjectAs(self, exists=False):
        '''
        Save back project to new pickle file.
        '''
        dlg = QFileDialog()
        fnm = dlg.getSaveFileName(self, 'Create a new project file...', filter='Pylot project (*.plp)')
        filename = fnm[0]
        if not len(fnm[0]):
            return False
        if not filename.split('.')[-1] == 'plp':
            filename = fnm[0] + '.plp'
        self.project.parameter = self._inputs
        self.project.save(filename)
        self.setDirty(False)
        self.saveProjectAsAction.setEnabled(True)
        self.update_status('Saved new project to {}'.format(filename), duration=5000)
        return True

    def saveProject(self, new=False):
        '''
        Save back project to pickle file.
        '''
        if self.project and not new:
            if not self.project.location:
                if not self.saveProjectAs(exists=True):
                    self.setDirty(True)
                    return False
            else:
                self.project.parameter = self._inputs
                self.project.save()
            if not self.project.dirty:
                self.update_status('Saved back project to file:\n{}'.format(self.project.location), duration=5000)
                self.setDirty(False)
                return True
            else:
                # if still dirty because saving failed
                qmb = QMessageBox.warning(self, 'Could not save project',
                                          'Could not save back to original file.\nChoose new file')
                self.setDirty(True)
        return self.saveProjectAs(exists=True)

    def draw(self):
        self.fill_eventbox()
        self.getPlotWidget().draw()

    def _setDirty(self):
        self.setDirty(True)

    def setDirty(self, value):
        self.saveProjectAction.setEnabled(value)
        self.saveProjectAsAction.setEnabled(True)
        self.project.setDirty(value)
        self.dirty = value

    def closeEvent(self, event):
        if self.okToContinue():
            event.accept()
        else:
            event.ignore()
            # self.closing.emit()
            # QMainWindow.closeEvent(self, event)

    def setParameter(self, show=True):
        if not self.paraBox:
            self.paraBox = PylotParaBox(self._inputs)
            self.paraBox.accepted.connect(self._setDirty)
            self.paraBox.accepted.connect(self.filterOptionsFromParameter)
        if show:
            self.paraBox.params_to_gui()
            self.paraBox.show()

    def PyLoTprefs(self):
        if not self._props:
            self._props = PropertiesDlg(self, infile=self.infile)

        if self._props.exec_():
            self.init_wfWidget()
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
        self.rootpath = None
        self.dirty = False
        self.parameter = None
        self._table = None

    def add_eventlist(self, eventlist):
        '''
        Add events from an eventlist containing paths to event directories.
        Will skip existing paths.
        '''
        if len(eventlist) == 0:
            return
        for item in eventlist:
            event = Event(item)
            event.rootpath = self.parameter['rootpath']
            event.database = self.parameter['database']
            event.datapath = self.parameter['datapath']
            if not event.path in self.getPaths():
                self.eventlist.append(event)
                self.setDirty()
            else:
                print('Skipping event with path {}. Already part of project.'.format(event.path))
        self.search_eventfile_info()

    def read_eventfile_info(self, filename, separator=','):
        '''
        Try to read event information from file (:param:filename) comparing specific event datetimes.
        File structure (each row): event, date, time, magnitude, latitude, longitude, depth
        separated by :param:separator each.
        '''
        infile = open(filename, 'r')
        for line in infile.readlines():
            event, date, time, mag, lat, lon, depth = line.split(separator)[:7]
            # skip first line
            try:
                month, day, year = date.split('/')
            except:
                continue
            year = int(year)
            # hardcoded, if year only consists of 2 digits (e.g. 16 instead of 2016)
            if year < 100:
                year += 2000
            datetime = '{}-{}-{}T{}'.format(year, month, day, time)
            try:
                datetime = UTCDateTime(datetime)
            except Exception as e:
                print(e, datetime, filename)
                continue
            for event in self.eventlist:
                if event.origins:
                    origin = event.origins[0]  # should have only one origin
                    if origin.time == datetime:
                        origin.latitude = float(lat)
                        origin.longitude = float(lon)
                        origin.depth = float(depth)
                elif not event.origins:
                    origin = Origin(resource_id=event.resource_id,
                                    time=datetime, latitude=float(lat),
                                    longitude=float(lon), depth=float(depth))
                    event.origins.append(origin)
                event.magnitudes.append(Magnitude(resource_id=event.resource_id,
                mag=float(mag),
                mag_type='M'))

    def search_eventfile_info(self):
        '''
        Search all datapaths in rootpath for filenames with given file extension fext
        and try to read event info from it
        '''
        datapaths = []
        fext = '.csv'
        for event in self.eventlist:
            if not event.datapath in datapaths:
                datapaths.append(event.datapath)
        for datapath in datapaths:
            datapath = os.path.join(self.rootpath, datapath)
            for filename in os.listdir(datapath):
                filename = os.path.join(datapath, filename)
                if os.path.isfile(filename) and filename.endswith(fext):
                    try:
                        self.read_eventfile_info(filename)
                    except Exception as e:
                        print('Failed on reading eventfile info from file {}: {}'.format(filename, e))

    def getPaths(self):
        '''
        Returns paths (eventlist) of all events saved in the project.
        '''
        paths = []
        for event in self.eventlist:
            paths.append(event.path)
        return paths

    def setDirty(self, value=True):
        self.dirty = value

    def getEventFromPath(self, path):
        '''
        Search for an event in the project by event path.
        '''
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
            self.setDirty(False)
        except Exception as e:
            print('Could not pickle PyLoT project. Reason: {}'.format(e))
            self.setDirty()

    @staticmethod
    def load(filename):
        '''
        Load project from filename.
        '''
        try:
            import cPickle
        except ImportError:
            import _pickle as cPickle
        infile = open(filename, 'rb')
        project = cPickle.load(infile)
        project.location = filename
        print('Loaded %s' % filename)
        return project


class getExistingDirectories(QFileDialog):
    '''
    File dialog with possibility to select multiple folders.
    '''

    def __init__(self, *args):
        super(getExistingDirectories, self).__init__(*args)
        self.setOption(self.DontUseNativeDialog, True)
        self.setOption(self.ReadOnly, True)
        self.setFileMode(self.Directory)
        self.setOption(self.ShowDirsOnly, True)
        self.findChildren(QListView)[0].setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.findChildren(QTreeView)[0].setSelectionMode(QAbstractItemView.ExtendedSelection)


def create_window():
    app_created = False
    app = QCoreApplication.instance()
    # check for existing app (when using ipython)
    if app is None:
        app = QApplication(sys.argv)
        app_created = True
    app.setOrganizationName("QtPyLoT");
    app.setOrganizationDomain("rub.de");
    app.setApplicationName("RUB");
    app.references = set()
    # app.references.add(window)
    # window.show()
    return app, app_created


def main(args=None):
    project_filename = None
    pylot_infile = None
    if args:
        if args.project_filename:
            project_filename = args.project_filename
        if args.input_filename:
            pylot_infile = args.input_filename

    # create the Qt application
    pylot_app, app_created = create_window()
    # pylot_app = QApplication(sys.argv)
    pixmap = QPixmap(":/splash/splash.png")
    splash = QSplashScreen(pixmap)
    splash.show()

    app_icon = QIcon()
    app_icon.addPixmap(QPixmap(':/icons/pylot.png'))

    # create the main window
    pylot_form = MainWindow(infile=pylot_infile)
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

    if project_filename:
        pylot_form.loadProject(args.project_filename)

    if app_created:
        pylot_app.exec_()
    else:
        return pylot_form


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Welcome to PyLoT.')
    parser.add_argument('-p', dest='project_filename', help='load project file',
                        default=None)
    parser.add_argument('-in', dest='input_filename', help='set pylot input file',
                        default=None)
    args = parser.parse_args()
    sys.exit(main(args))
