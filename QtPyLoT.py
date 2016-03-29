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

import os
import sys

import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from PySide.QtCore import QCoreApplication, QSettings, Signal, QFile, \
    QFileInfo, Qt
from PySide.QtGui import QMainWindow, QInputDialog, QIcon, QFileDialog, \
    QWidget, QHBoxLayout, QStyle, QKeySequence, QLabel, QFrame, QAction, \
    QDialog, QErrorMessage, QApplication, QPixmap, QMessageBox, QSplashScreen, \
    QActionGroup, QListWidget, QDockWidget
import numpy as np
from obspy import UTCDateTime

from pylot.core.read.data import Data
from pylot.core.read.inputs import FilterOptions, AutoPickParameter
from pylot.core.pick.autopick import autopickevent
from pylot.core.loc.nll import locate as locateNll
from pylot.core.util.defaults import FILTERDEFAULTS
from pylot.core.util.errors import FormatError, DatastructureError, \
    OverwriteError
from pylot.core.util.connection import checkurl
from pylot.core.util.utils import fnConstructor, createEvent, getLogin, \
    createCreationInfo, getGlobalTimes
from pylot.core.util.widgets import FilterOptionsDialog, NewEventDlg, \
    MPLWidget, PropertiesDlg, HelpForm, createAction, PickDlg
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.thread import AutoPickThread
from pylot.core.util.version import get_git_version as _getVersionString
import icons_rc

locateTool = dict(nll=locateNll)


class MainWindow(QMainWindow):
    __version__ = _getVersionString()
    closing = Signal()

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        self.createAction = createAction
        settings = QSettings()
        if settings.value("user/FullName", None) is None:
            fulluser = QInputDialog.getText(self, "Enter Name:", "Full name")
            settings.setValue("user/FullName", fulluser)
            settings.setValue("user/Login", getLogin())
        if settings.value("agency_id", None) is None:
            agency = QInputDialog.getText(self,
                                          "Enter authority name (e.g. BUG):",
                                          "Authority")
            settings.setValue("agency_id", agency)
        self.recentEvents = settings.value("data/recentEvents", [])
        self.fname = None
        self.fnames = None
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
        if self.recentEvents:
            lastEvent = self.getLastEvent()
            self.data = Data(self, lastEvent)
        else:
            self.data = Data(self)

        # load and display waveform data
        self.dirty = False
        self.loadData()
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
        self.DataPlot = MPLWidget(parent=self, xlabel=xlab, ylabel=None,
                                  title=plottitle)
        self.DataPlot.mpl_connect('button_press_event',
                                  self.pickOnStation)
        self.DataPlot.mpl_connect('axes_enter_event',
                                  lambda event: self.tutorUser())
        _layout.addWidget(self.DataPlot)

        openIcon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
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

        newEventAction = self.createAction(self, "&New event ...",
                                           self.createNewEvent,
                                           QKeySequence.New, newIcon,
                                           "Create a new event.")
        openEventAction = self.createAction(self, "&Open event ...",
                                            self.loadData, QKeySequence.Open,
                                            openIcon, "Open an event.")
        openEventAction.setData(None)
        saveEventAction = self.createAction(self, "&Save event ...",
                                            self.saveData, QKeySequence.Save,
                                            saveIcon, "Save actual event data.")
        openWFDataAction = self.createAction(self, "Open &waveforms ...",
                                             self.loadWaveformData,
                                             "Ctrl+W", QIcon(":/wfIcon.png"),
                                             """Open waveform data (event will
                                             be closed).""")
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
        printAction = self.createAction(self, "&Print event ...",
                                        self.printEvent, QKeySequence.Print,
                                        print_icon,
                                        "Print waveform overview.")
        helpAction = self.createAction(self, "&Help ...", self.helpHelp,
                                       QKeySequence.HelpContents, helpIcon,
                                       """Show either the documentation
                                       homepage (internet connection available),
                                       or shipped documentation files.""")
        self.fileMenu = self.menuBar().addMenu('&File')
        self.fileMenuActions = (newEventAction, openEventAction,
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
        fileToolActions = (newEventAction, openEventAction, saveEventAction)
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
        autoPickActions = (auto_pick,)
        self.addActions(autoPickToolBar, autoPickActions)

        # pickToolBar = self.addToolBar("PickTools")
        # pickToolActions = (selectStation, )
        # pickToolBar.setObjectName("PickTools")
        # self.addActions(pickToolBar, pickToolActions)

        locateEvent = self.createAction(parent=self, text='locateEvent',
                                        slot=self.locateEvent, shortcut='Alt+Ctrl+L',
                                        icon=locate_icon, tip='Locate the event using '
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
        for eventID in self.recentEvents:
            fname = fnConstructor(eventID)
            if eventID != current and QFile.exists(fname):
                recentEvents.append(eventID)
                recentEvents.reverse()
                self.recentEvents = recentEvents[0:5]
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
                             self.loadData)
                self.fileMenu.addAction(action)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.fileMenuActions[-1])

    def getRoot(self):
        settings = QSettings()
        return settings.value("data/dataRoot")

    def loadAutoPicks(self):
        self.loadData(type='auto')

    def loadData(self, fname=None, type='manual'):
        if not self.okToContinue():
            return
        if fname is None:
            action = self.sender()
            if isinstance(action, QAction):
                if action.data() is None:
                    filt = "Supported event formats (*.mat *.qml *.xml *.kor *.evt)"
                    caption = "Open an event file"
                    fname = QFileDialog().getOpenFileName(self,
                                                          caption=caption,
                                                          filter=filt)
                    fname = fname[0]
                else:
                    fname = unicode(action.data().toString())
        self.setFileName(fname)
        self.data += Data(self, evtdata=self.getFileName())
        self.updatePicks(type=type)
        self.drawPicks()

    def getLastEvent(self):
        return self.recentEvents[0]

    def addRecentEvent(self, event):
        self.recentEvents.insert(0, event)

    def getWFFnames(self):
        try:
            evt = self.getData().getEvtData()
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

    def getFileName(self):
        return self.fname

    def setFileName(self, fname):
        if self.getFileName() is not None:
            self.addRecentEvent(self.getFileName())
        self.fname = fname

    def getEventFileName(self):
        if self.getFileName() is None:
            self.setFileName(self.getData().getEventFileName())
        return self.getFileName()

    def saveData(self):

        def getSavePath(e):
            print('warning: {0}'.format(e))
            directory = os.path.join(self.getRoot(), self.getEventFileName())
            file_filter = "QuakeML file (*.xml);;VELEST observation file format (*.cnv);;NonLinLoc observation file (*.obs)"
            fname, selected_filter = QFileDialog.getSaveFileName(self, 'Save event data ...',
                                                directory, file_filter)

            fbasename, exform = os.path.splitext(fname)

            if not exform and selected_filter:
                exform = selected_filter.split('*')[1][:-1]

            return fbasename, exform

        settings = QSettings()
        fbasename = self.getEventFileName()
        exform = settings.value('data/exportFormat', 'QUAKEML')
        try:
            self.getData().applyEVTData(self.getPicks())
        except OverwriteError:
            msgBox = QMessageBox()
            msgBox.setText("Picks have been modified!")
            msgBox.setInformativeText("Do you want to save the changes and overwrite the picks?")
            msgBox.setDetailedText(self.getData().getPicksStr())
            msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Cancel)
            msgBox.setDefaultButton(QMessageBox.Save)
            ret = msgBox.exec_()
            if ret == QMessageBox.Save:
                self.getData().resetPicks()
                return self.saveData()
            elif ret == QMessageBox.Cancel:
                return False
        try:
            self.getData().exportEvent(fbasename, exform)
        except FormatError as e:
            fbasename, exform = getSavePath(e)
        except AttributeError as e:
            fbasename, exform = getSavePath(e)
        if not fbasename:
            return False
        self.getData().exportEvent(fbasename, exform)
        self.setDirty(False)
        self.updateStatus('Event saved as %s' % (fbasename + exform))
        return True

    def getComponent(self):
        return self.dispComponent

    def setComponent(self, component):
        self.dispComponent = component

    def getData(self):
        return self.data

    def getPicks(self, type='manual'):
        rdict = dict(auto=self.autopicks, manual=self.picks)
        return rdict[type]

    def getPicksOnStation(self, station, type='manual'):
        try:
            return self.getPicks(type)[station]
        except KeyError:
            return None

    def getPlotWidget(self):
        return self.DataPlot

    @staticmethod
    def getWFID(gui_event):

        ycoord = gui_event.ydata

        statID = int(round(ycoord))

        return statID

    def getStationID(self, station):
        for wfID in self.getPlotWidget().getPlotDict().keys():
            actual_station = self.getPlotWidget().getPlotDict()[wfID][0]
            if station == actual_station:
                return wfID
        return None

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
        if ans:
            self.plotWaveformData()
            return ans
        else:
            return ans

    def plotWaveformData(self):
        zne_text = {'Z': 'vertical', 'N': 'north-south', 'E': 'east-west'}
        comp = self.getComponent()
        title = 'section: {0} components'.format(zne_text[comp])
        wfst = self.getData().getWFData().select(component=comp)
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

    def plotN(self):
        self.setComponent('N')
        self.plotWaveformData()
        self.drawPicks()

    def plotE(self):
        self.setComponent('E')
        self.plotWaveformData()
        self.drawPicks()

    def pushFilterWF(self, param_args):
        self.getData().filterWFData(param_args)

    def filterWaveformData(self):
        if self.getData():
            if self.getFilterOptions() and self.filterAction.isChecked():
                kwargs = self.getFilterOptions().parseFilterOptions()
                self.pushFilterWF(kwargs)
            elif self.filterAction.isChecked():
                self.adjustFilterOptions()
            else:
                self.getData().resetWFData()
        self.plotWaveformData()
        self.drawPicks()

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
            self.updateStatus('Error ...')
            emsg = QErrorMessage(self)
            emsg.showMessage('Error: {0}'.format(e))
        else:
            self.updateStatus('Filter loaded ... '
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
        self.updateStatus('Seismic phase changed to '
                          '{0}'.format(self.getSeismicPhase()))

    def pickOnStation(self, gui_event):

        wfID = self.getWFID(gui_event)

        station = self.getStationName(wfID)
        self.updateStatus('picking on station {0}'.format(station))
        data = self.getData().getWFData()
        pickDlg = PickDlg(self, data=data.select(station=station),
                          station=station,
                          picks=self.getPicksOnStation(station))
        if pickDlg.exec_():
            self.setDirty(True)
            self.updateStatus('picks accepted ({0})'.format(station))
            replot = self.addPicks(station, pickDlg.getPicks())
            if replot:
                self.plotWaveformData()
                self.drawPicks()
            else:
                self.drawPicks(station)
        else:
            self.updateStatus('picks discarded ({0})'.format(station))
        if not self.getLocflag() and self.check4Loc():
            self.setLocflag(True)
        elif self.getLocflag() and not self.check4Loc():
            self.setLocflag(False)

    def addListItem(self, text):
        self.listWidget.addItem(text)
        self.listWidget.scrollToBottom()

    def autoPick(self):
        self.listWidget = QListWidget()
        self.setDirty(True)
        self.logDockWidget = QDockWidget("AutoPickLog", self)
        self.logDockWidget.setObjectName("LogDockWidget")
        self.logDockWidget.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.logDockWidget.setWidget(self.listWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.logDockWidget)
        self.addListItem('loading default values for local data ...')
        home = os.path.expanduser("~")
        autopick_parameter = AutoPickParameter('%s/.pylot/autoPyLoT_local.in' % home)
        self.addListItem(str(autopick_parameter))

        # Create the worker thread and run it
        self.thread = AutoPickThread(parent=self,
                                     func=autopickevent,
                                     data=self.getData().getWFData(),
                                     param=autopick_parameter)
        self.thread.message.connect(self.addListItem)
        self.thread.start()
        self.thread.finished.connect(self.finalizeAutoPick)

    def finalizeAutoPick(self):
        self.drawPicks(picktype='auto')
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
        evt = self.getData().getEvtData()
        picks = {}
        for pick in evt.picks:
            phase = {}
            station = pick.waveform_id.station_code
            try:
                onsets = picks[station]
            except KeyError as e:
                print(e)
                onsets = {}
            mpp = pick.time
            lpp = mpp + pick.time_errors.upper_uncertainty
            epp = mpp - pick.time_errors.lower_uncertainty
            spe = pick.time_errors.uncertainty
            phase['mpp'] = mpp
            phase['epp'] = epp
            phase['lpp'] = lpp
            phase['spe'] = spe
            try:
                picker = str(pick.method_id)
                if picker.startswith('smi:local/'):
                    picker = picker.split('smi:local/')[1]
                phase['picker'] = picker
            except IndexError:
                pass

            onsets[pick.phase_hint] = phase.copy()
            picks[station] = onsets.copy()
        if type == 'manual':
            self.picks.update(picks)
        elif type == 'auto':
            self.autopicks.update(picks)

    def drawPicks(self, station=None, picktype='manual'):
        # if picks to draw not specified, draw all picks available
        if not station:
            for station in self.getPicks(type=picktype):
                self.drawPicks(station, picktype=picktype)
            return
        # plotting picks
        plotID = self.getStationID(station)
        ax = self.getPlotWidget().axes
        ylims = np.array([-.5, +.5]) + plotID
        phase_col = {'P': ('c', 'c--', 'b-', 'bv', 'b^'),
                     'S': ('m', 'm--', 'r-', 'rv', 'r^')}

        stat_picks = self.getPicks(type=picktype)[station]

        for phase in stat_picks:
            picks = stat_picks[phase]
            if type(stat_picks[phase]) is not dict:
                return
            colors = phase_col[phase[0].upper()]

            stime = getGlobalTimes(self.getData().getWFData())[0]

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
        self.draw()

    def locateEvent(self):
        settings = QSettings()
        loctool = settings.value("loc/tool", "nll")
        extlocpath = settings.value("%s/binPath".format(loctool), None)
        locroot = settings.value("%s/rootPath".format(loctool), None)
        if extlocpath is None or locroot is None:
            self.PyLoTprefs()

    def check4Loc(self):
        return self.picksNum() > 4

    def picksNum(self):
        num = 0
        for phases in self.getPicks().values():
            num += len(phases)
        return num

    def getLocflag(self):
        return self.loc

    def setLocflag(self, value):
        self.loc = value

    def updateStatus(self, message, duration=5000):
        self.statusBar().showMessage(message, duration)
        if self.getData() is not None:
            if not self.getData().isNew():
                self.setWindowTitle(
                        "PyLoT - processing event %s[*]" % self.getData().getID())
            elif self.getData().isNew():
                self.setWindowTitle("PyLoT - New event [*]")
            else:
                self.setWindowTitle(
                        "PyLoT - seismic processing the python way[*]")
        self.setWindowModified(self.dirty)

    def tutorUser(self):
        self.updateStatus('select trace to pick on station ...', 10000)

    def printEvent(self):
        pass

    def createNewEvent(self):
        if self.okToContinue():
            new = NewEventDlg()
            if new.exec_() != QDialog.Rejected:
                evtpar = new.getValues()
                cinfo = createCreationInfo(agency_id=self.agency)
                event = createEvent(evtpar['origintime'], cinfo)
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
