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

import os, sys
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
from obspy.core import UTCDateTime

from pylot.core.read.data import Data
from pylot.core.read.inputs import FilterOptions, AutoPickParameter
from pylot.core.pick.autopick import autopickevent
from pylot.core.util.defaults import FILTERDEFAULTS
from pylot.core.util.errors import FormatError, DatastructureError
from pylot.core.util.connection import checkurl
from pylot.core.util.utils import fnConstructor, createEvent, getLogin,\
    createCreationInfo, getGlobalTimes
from pylot.core.util.widgets import FilterOptionsDialog, NewEventDlg,\
    MPLWidget, PropertiesDlg, HelpForm, createAction, PickDlg
from pylot.core.util.structure import DATASTRUCTURE
from pylot.core.util.thread import WorkerThread
from pylot.core.util.version import get_git_version as _getVersionString
import icons_rc

# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):
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
        self.loadWaveformData()
        self.loadData()
        self.updateFilterOptions()

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

        phaseToolBar = self.addToolBar("PhaseTools")
        phaseToolActions = (self.selectPAction, self.selectSAction)
        phaseToolBar.setObjectName("PhaseTools")
        self.addActions(phaseToolBar, phaseToolActions)

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
                                                          ' displayed!',
                                      checkable=False)

        autoPickToolBar = self.addToolBar("autoPyLoT")
        autoPickActions = (auto_pick,)
        self.addActions(autoPickToolBar, autoPickActions)

        # pickToolBar = self.addToolBar("PickTools")
        # pickToolActions = (selectStation, )
        # pickToolBar.setObjectName("PickTools")
        # self.addActions(pickToolBar, pickToolActions)

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

    def loadData(self, fname=None):
        if fname is None:
            try:
                self.data = Data(self, evtdata=self.fname)
            except AttributeError:
                action = self.sender()
                if isinstance(action, QAction):
                    if action.data() is None:
                        filt = "Supported event formats (*.mat *.qml *.xml *.kor *.evt)"
                        caption = "Open an event file"
                        fname = QFileDialog().getOpenFileName(self,
                                                              caption=caption,
                                                              filter=filt)
                        self.fname = fname[0]
                    else:
                        self.fname = unicode(action.data().toString())
                if not self.okToContinue():
                    return
        else:
            self.fname = fname
            self.data = Data(self, evtdata=self.fname)

    def getLastEvent(self):
        return self.recentEvents[0]

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
            return self.fnames
        except DatastructureError, e:
            print e
            props = PropertiesDlg(self)
            if props.exec_() == QDialog.Accepted:
                return self.getWFFnames()
            else:
                return

    def getEventFileName(self):
        return self.getData().getEventFileName()

    def saveData(self):
        settings = QSettings()
        exform = settings.value('data/exportFormat', 'QUAKEML')
        self.getData().applyEVTData(self.getPicks())
        try:
            self.getData().exportEvent(self.fname, exform)
        except FormatError:
            return False
        except AttributeError, e:
            print 'warning: {0}'.format(e)
            directory = os.path.join(self.getRoot(), self.getEventFileName())
            file_filter = "Seismic observation files (*.cnv *.obs *.xml)"
            fname = QFileDialog.getSaveFileName(self, 'Save event data ...',
                                                directory, file_filter)
            fbasename, exform = os.path.splitext(fname[0])
            if not fbasename:
                return False
            self.getData().exportEvent(fbasename, exform)
        return True

    def getComponent(self):
        return self.dispComponent

    def setComponent(self, component):
        self.dispComponent = component

    def getData(self):
        return self.data

    def getPicks(self):
        return self.picks

    def getPicksOnStation(self, station):
        try:
            return self.getPicks()[station]
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
            self.data.setWFData(self.fnames)
        elif self.fnames is None and self.okToContinue():
            self.data.setWFData(self.getWFFnames())
        self.plotWaveformData()

    def plotWaveformData(self):
        zne_text = {'Z': 'vertical', 'N': 'north-south', 'E': 'east-west'}
        comp = self.getComponent()
        title = 'overview: {0} components'.format(zne_text[comp])
        wfst = self.getData().getWFData().select(component=comp)
        self.getPlotWidget().plotWFData(wfdata=wfst, title=title)
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

    def filterWaveformData(self):
        if self.getData():
            def hasfreq(kwdict):
                for key in kwdict.keys():
                    if not key.startswith('freq'):
                        return True
                return False

            if self.filterAction.isChecked():
                kwargs = {}
                freq = self.getFilterOptions().getFreq()
                if freq is not None and len(freq) > 1:
                    kwargs['freqmin'] = freq[0]
                    kwargs['freqmax'] = freq[1]
                elif freq is not None and len(freq) == 1:
                    kwargs['freq'] = freq
                if hasfreq(kwargs):
                    kwargs['type'] = self.getFilterOptions().getFilterType()
                    kwargs['corners'] = self.getFilterOptions().getOrder()
                    self.getData().filterWFData(kwargs)
            else:
                self.getData().resetWFData()
        self.plotWaveformData()

    def adjustFilterOptions(self):
        filteroptions = self.getFilterOptions()
        fstring = "Filter Options ({0})".format(self.getSeismicPhase())
        filterDlg = FilterOptionsDialog(titleString=fstring,
                                        parent=self,
                                        filterOptions=filteroptions)
        if filterDlg.exec_():
            filteroptions = filterDlg.getFilterOptions()
            self.setFilterOptions(filteroptions)

    def getFilterOptions(self):
        try:
            return self.filteroptions[self.getSeismicPhase()]
        except AttributeError, e:
            print e
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
                for key, value in FILTERDEFAULTS.iteritems():
                    self.setFilterOptions(FilterOptions(**value), key)
            elif settings.value("filterdefaults", None) is not None:
                for key, value in settings.value("filterdefaults"):
                    self.setFilterOptions(FilterOptions(**value), key)
        except Exception, e:
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
            self.addPicks(station, pickDlg.getPicks())
            self.drawPicks(station)
        else:
            self.updateStatus('picks discarded ({0})'.format(station))

    def autoPick(self):
        list = QListWidget()
        self.setDirty(True)
        logDockWidget = QDockWidget("AutoPickLog", self)
        logDockWidget.setObjectName("LogDockWidget")
        logDockWidget.setAllowedAreas(Qt.LeftDockWidgetArea)
        logDockWidget.setWidget(list)
        logDockWidget.show()
        logDockWidget.setFloating(False)
        list.addItem('loading default values for local data ...')
        autopick_parameter = AutoPickParameter('autoPyLoT_local.in')
        list.addItem(str(autopick_parameter))

        # Create the worker thread and run it
        self.thread = WorkerThread(parent=self,
                                   func=autopickevent,
                                   data=self.getData().getWFData(),
                                   param=autopick_parameter)
        self.thread.message.connect(list.addItem)
        self.thread.start()

        self.drawPicks()


    def addPicks(self, station, picks):
        stat_picks = self.getPicksOnStation(station)
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
            elif ret == QMessageBox.Cancel:
                pass
            else:
                raise Exception('FATAL: Should never occur!')
        self.getPicks()[station] = stat_picks

    def drawPicks(self, station=None):
        # if picks to draw not specified, draw all picks available
        if not station:
            for station in self.getPicks():
                self.drawPicks(station)
            return
        # plotting picks
        plotID = self.getStationID(station)
        ax = self.getPlotWidget().axes
        ylims = np.array([-.5, +.5]) + plotID
        phase_col = {'P': ('c', 'c--', 'b-'),
                     'S': ('m', 'm--', 'r-')}

        stat_picks = self.getPicks()[station]

        for phase in stat_picks:

            picks = stat_picks[phase]
            colors = phase_col[phase[0].upper()]

            stime = getGlobalTimes(self.getData().getWFData())[0]

            mpp = picks['mpp'] - stime
            epp = picks['epp'] - stime
            lpp = picks['lpp'] - stime
            spe = picks['spe']

            ax.fill_between([epp, lpp], ylims[0], ylims[1],
                            alpha=.5, color=colors[0])
            ax.plot([mpp - spe, mpp - spe], ylims, colors[1],
                    [mpp, mpp], ylims, colors[2],
                    [mpp + spe, mpp + spe], ylims, colors[1])
        self.draw()

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

    # set Application Information
    pylot_app.setOrganizationName("Ruhr-University Bochum / MAGS2")
    pylot_app.setOrganizationDomain("rub.de")
    pylot_app.setApplicationName("PyLoT")
    pylot_app.setApplicationVersion(__version__)
    pylot_app.setWindowIcon(app_icon)

    # create the main window
    pylot_form = MainWindow()
    splash.showMessage('Loading. Please wait ...')
    pylot_app.processEvents()

    # Show main window and run the app
    pylot_form.showMaximized()
    splash.finish(pylot_form)
    pylot_app.exec_()


if __name__ == "__main__":
    sys.exit(main())
