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

import sys

from PySide.QtCore import QCoreApplication, QSettings, Signal, QFile, QFileInfo
from PySide.QtGui import QMainWindow, QInputDialog, QIcon, QFileDialog, \
    QWidget, QHBoxLayout, QStyle, QKeySequence, QLabel, QFrame, QAction, \
    QDialog, QErrorMessage, QApplication
from obspy.core import UTCDateTime

from pylot.core.read import Data, FilterOptions
from pylot.core.util import _getVersionString, FILTERDEFAULTS, fnConstructor, \
    checkurl, FormatError, FilterOptionsDialog, \
    NewEventDlg, createEvent, MPLWidget, PropertiesDlg, HelpForm, \
    DatastructureError, createAction, getLogin, createCreationInfo
from pylot.core.util.structure import DATASTRUCTURE



# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):
    closing = Signal()

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)


        self.createAction = createAction
        self.dirty = False
        settings = QSettings()
        if settings.value("user/FullName", None) is None:
            fulluser = QInputDialog.getText(self, "Enter Name:", "Full name")
            settings.setValue("user/FullName", fulluser)
            settings.setValue("user/Login", getLogin())
        if settings.value("agency_id", None) is None:
            agency = QInputDialog.getText(self, "Enter authority name (e.g. BUG):", "Authority")
        self.recentEvents = settings.value("data/recentEvents", [])
        self.fnames = None
        self.dataStructure = DATASTRUCTURE[
            settings.value("data/Structure", "PILOT")]()
        self.seismicPhase = str(settings.value("phase", "P"))
        self.dispComponent = str(settings.value("plotting/dispComponent", "Z"))
        if settings.value("data/dataRoot", None) is None:
            dirname = QFileDialog().getExistingDirectory(
                caption='Choose data root ...')
            settings.setValue("data/dataRoot", dirname)
        settings.sync()

        self.filteroptions = {}

        # UI has to be set up before(!) children widgets are about to show up
        self.setupUi()

        # initialize event data
        if self.recentEvents:
            lastEvent = self.getLastEvent()
            self.data = Data(self, lastEvent)
        else:
            self.data = Data(self)

        # load and display waveform data
        self.loadWaveformData()
        self.dirty = False
        self.loadData()
        self.updateFilterOptions()


    def setupUi(self):

        try:
            self.startTime = min(
                [tr.stats.starttime for tr in self.data.wfdata])
        except:
            self.startTime = UTCDateTime()

        self.setWindowTitle("PyLoT - do seismic processing the python way")
        self.setWindowIcon(QIcon(":/icon.ico"))

        xlab = self.startTime.strftime('seconds since %d %b %Y %H:%M:%S (%Z)')

        _widget = QWidget()
        _layout = QHBoxLayout()

        plottitle = "Overview: {0} components ".format(self.getComponent())

        # create central matplotlib figure canvas widget
        self.DataPlot = MPLWidget(parent=self, xlabel=xlab, ylabel=None,
                                  title=plottitle)

        _layout.addWidget(self.DataPlot)

        openIcon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
        quitIcon = self.style().standardIcon(QStyle.SP_MediaStop)
        saveIcon = self.style().standardIcon(QStyle.SP_DriveHDIcon)
        helpIcon = self.style().standardIcon(QStyle.SP_DialogHelpButton)
        newIcon = self.style().standardIcon(QStyle.SP_FileIcon)
        newEventAction = self.createAction(self, "&New event ...",
                                           self.createNewEvent,
                                           QKeySequence.New, newIcon,
                                           "Create a new event.")
        openEventAction = self.createAction(self, "&Open event ...", self.loadData,
                                            QKeySequence.Open, openIcon,
                                            "Open an event.")
        openEventAction.setData(None)
        saveEventAction = self.createAction(self, "&Save event ...", self.saveData,
                                            QKeySequence.Save, saveIcon,
                                            "Save actual event data.")
        openWFDataAction = self.createAction(self, "Open &waveforms ...",
                                             self.loadWaveformData,
                                             "Ctrl+W", QIcon(":/wfIcon.png"),
                                             """Open waveform data (event will
                                             be closed).""")

        prefsEventAction = self.createAction(self, "Preferences", self.PyLoTprefs,
                                             QKeySequence.Preferences,
                                             QIcon(None),
                                             "Edit PyLoT app preferences.")
        quitAction = self.createAction(self, "&Quit",
                                       QCoreApplication.instance().quit,
                                       QKeySequence.Close, quitIcon,
                                       "Close event and quit PyLoT")
        self.filterAction = self.createAction(self, "&Filter ...", self.filterWaveformData,
                                         "Ctrl+F", QIcon(":/filter.png"),
                                         """Toggle un-/filtered waveforms
                                         to be displayed, according to the
                                         desired seismic phase.""", True)
        filterEditAction = self.createAction(self, "&Filter parameter ...",
                                             self.adjustFilterOptions,
                                             "Alt+F", QIcon(None),
                                             """Adjust filter parameters.""")
        self.selectPAction = self.createAction(self, "&P", self.alterPhase, "Alt+P",
                                          QIcon(":/picon.png"),
                                          "Toggle P phase.", True)
        self.selectSAction = self.createAction(self, "&S", self.alterPhase, "Alt+S",
                                          QIcon(":/sicon.png"),
                                          "Toggle S phase", True)
        printAction = self.createAction(self, "&Print event ...",
                                        self.printEvent, QKeySequence.Print,
                                        QIcon(":/printer.png"),
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
        helpActions = (helpAction, )
        self.addActions(self.helpMenu, helpActions)

        fileToolBar = self.addToolBar("FileTools")
        fileToolActions = (newEventAction, openEventAction, saveEventAction)
        fileToolBar.setObjectName("FileTools")
        self.addActions(fileToolBar, fileToolActions)

        phaseToolBar = self.addToolBar("PhaseTools")
        phaseToolActions = (self.selectPAction, self.selectSAction)
        phaseToolBar.setObjectName("PhaseTools")
        self.addActions(phaseToolBar, phaseToolActions)

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 500)

        _widget.setLayout(_layout)
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
                action = QAction(QIcon(":/icon.png"),
                                 "&{0} {1}".format(i + 1,
                                                   QFileInfo(fname).fileName()),
                                 self)
                action.setData(fname)
                self.connect(action, Signal("triggered()"),
                             self.loadData)
                self.fileMenu.addAction(action)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.fileMenuActions[-1])

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

    def saveData(self):
        settings = QSettings()
        exform = settings.value('data/exportFormat', 'None')
        try:
            self.data.exportEvent(self.fname, exform)
        except FormatError:
            return False
        except AttributeError, e:
            print 'warning: {0}'.format(e)
            fname = QFileDialog.getSaveFileName(self, 'Save event')
            fname = fname[0]
            self.data.exportEvent(fname, exform)
        return True

    def getComponent(self):
        return self.dispComponent

    def getData(self):
        return self.data

    def getPlotWidget(self):
        return self.DataPlot

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
            self.dirty = True
            self.data.setWFData(self.fnames)
        elif self.fnames is None and self.okToContinue():
            self.data.setWFData(self.getWFFnames())
        self.plotWaveformData()

    def plotWaveformData(self):
        self.getData().plotWFData(self.getPlotWidget())

    def filterWaveformData(self):
        if self.getData():
            def hasfreq(kwargs):
                for key in kwargs.keys():
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
        filteroptions = None
        fstring = "Filter Options ({0})".format(self.getSeismicPhase())
        filterDlg = FilterOptionsDialog(titleString=fstring,
                                        parent=self,
                                        filterOptions=self.getFilterOptions())
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
            if settings.value("filterdefaults", None) is None and not self.getFilters():
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
                              '[{0}: {1} Hz]'.format(self.getFilterOptions().getFilterType(), self.getFilterOptions().getFreq()))
        if self.filterAction.isChecked():
            self.filterWaveformData()

    def getSeismicPhase(self):
        return self.seismicPhase

    def alterPhase(self):
        pass

    def setSeismicPhase(self, phase):
        self.seismicPhase = self.seismicPhaseButtonGroup.getValue()
        self.updateStatus('Seismic phase changed to '
                          '{0}'.format(self.getSeismicPhase()))

    def updateStatus(self, message):
        self.statusBar().showMessage(message, 5000)
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


    def printEvent(self):
        pass

    def createNewEvent(self):
        if self.okToContinue():
            new = NewEventDlg()
            if new.exec_() != QDialog.Rejected:
                evtpar = new.getValues()
                cinfo = createCreationInfo(agency_id=self.agency)
                event = createEvent(evtpar['origintime'])
                self.data = Data(self, evtdata=createEvent(**evtpar))
                self.dirty = True

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
    pylot_app = QApplication(sys.argv[0])

    # set Application Information
    pylot_app.setOrganizationName("Ruhr-University Bochum / MAGS2")
    pylot_app.setOrganizationDomain("rub.de")
    pylot_app.setApplicationName("PyLoT")
    pylot_app.setWindowIcon(QIcon(":/icon.ico"))

    # create the main window
    pylot_form = MainWindow()

    # Show main window and run the app
    pylot_form.show()
    pylot_app.exec_()


if __name__ == "__main__":
    sys.exit(main())
