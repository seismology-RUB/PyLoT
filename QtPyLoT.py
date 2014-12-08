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

from PySide.QtCore import *
from PySide.QtGui import *
from obspy.core import (UTCDateTime)

from pylot.core.util import _getVersionString
from pylot.core.read import (Data,
                             FilterOptions)
from pylot.core.util import FILTERDEFAULTS
from pylot.core.util import fnConstructor
from pylot.core.util import checkurl
from pylot.core.util import layoutStationButtons
from pylot.core.util import (FilterOptionsDialog,
                             MPLWidget,
                             HelpForm)


# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):

    closing = Signal()

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        settings = QSettings()
        if settings.value("user/FullName", None) is None:
            fulluser = QInputDialog.getText(self, "Enter Name:", "Full name")
            settings.setValue("user/FullName", fulluser)
            settings.setValue("user/Login", os.getlogin())
            settings.sync()
        self.recentEvents = settings.value("data/recentEvents", [])
        self.setWindowTitle("PyLoT - do seismic processing the python way")
        self.setWindowIcon(QIcon(":/icon.ico"))
        self.seismicPhase = str(settings.value("phase", "P"))
        if settings.value("data/dataRoot", None) is None:
            dirname = QFileDialog().getExistingDirectory(caption = 'Choose data root ...')
            settings.setValue("data/dataRoot", dirname)
            settings.sync()

        # initialize filter parameter
        filterOptionsP = FILTERDEFAULTS['P']
        filterOptionsS = FILTERDEFAULTS['S']
        # print filterOptionsP, "\n", filterOptionsS
        self.filterOptionsP = FilterOptions(**filterOptionsP)
        self.filterOptionsS = FilterOptions(**filterOptionsS)

        # initialize data
        self.data = None
        self.dirty = False
        self.loadData()
        self.updateFilterOptions()
        # print self.filteroptions
        try:
            self.startTime = min([tr.stats.starttime for tr in self.data.wfdata])
        except:
            self.startTime = UTCDateTime()

        self.setupUi()

    def _getCurrentPlotType(self):
        return 'TestType'

    def createAction(self, text, slot=None, shortcut=None, icon=None,
                     tip=None, checkable=False):
        """
        :rtype : ~PySide.QtGui.QAction
        """
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(icon)
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action

    def updateFileMenu(self):

        self.fileMenu.clear()
        self.addActions(self.fileMenu, self.fileMenuActions[:-1])
        current = self.data.evtdata.getID()
        recentEvents = []
        for eventID in self.recentEvents:
            fname = fnConstructor(eventID)
            if eventID != current and QFile.exists(fname):
                recentEvents.append(eventID)
        if recentEvents:
            self.fileMenu.addSeparator()
            for i, eventID in enumerate(recentEvents):
                fname = fnConstructor(eventID)
                action = QAction(QIcon(":/icon.png"),
                                 "&{0} {1}".format(i + 1,
                                                   QFileInfo(fname).fileName()),
                                 self)
                action.setData(fname)
                self.connect(action, SIGNAL("triggered()"),
                             self.loadData)
                self.fileMenu.addAction(action)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.fileMenuActions[-1])

    def loadData(self, fname=None):
        if fname is None:
            action = self.sender()
            if isinstance(action, QAction):
                if action.data() is None:
                    fname = QFileDialog()
                else:
                    fname = unicode(action.data().toString())
                if not self.okToContinue():
                    return
            else:
                return
        if fname:
            self.data = Data(evtdata=fname)

    def saveData(self):
        return True

    def getComponent(self):
        return self

    def getData(self):
        return self.data

    def getDataWidget(self):
        return self.DataPlot

    def setupUi(self):
        self.setWindowIcon(QIcon(":/icon.ico"))

        xlab = self.startTime.strftime('seconds since %d %b %Y %H:%M:%S (%Z)')
        plottitle = self._getCurrentPlotType()

        _widget = QWidget()
        _layout = QHBoxLayout()

        # create central matplotlib figure widget
        self.DataPlot = MPLWidget(parent=self,
                                  xlabel=xlab,
                                  ylabel=None,
                                  title=plottitle)
        statsButtons = layoutStationButtons(self.getData(), self.getComponent())
        _layout.addLayout(statsButtons)
        _layout.addWidget(self.DataPlot)

        openIcon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
        quitIcon = self.style().standardIcon(QStyle.SP_MediaStop)
        saveIcon = self.style().standardIcon(QStyle.SP_DriveHDIcon)
        openEventAction = self.createAction("&Open event ...", self.loadData,
                                            QKeySequence.Open, openIcon,
                                            "Open an event.")
        openEventAction.setData(None)
        saveEventAction = self.createAction("&Save event ...", self.saveData,
                                            QKeySequence.Save, saveIcon,
                                            "Save actual event data.")
        quitAction = self.createAction("&Quit",
                                       QCoreApplication.instance().quit,
                                       QKeySequence.Close, quitIcon,
                                       "Close event and quit PyLoT")
        filterAction = self.createAction("&Filter ...", self.filterData,
                                         "Ctrl+F", QIcon(":/filter.png"),
                                         """Toggle un-/filtered waveforms 
                                         to be displayed, according to the 
                                         desired seismic phase.""", True)
        filterEditAction = self.createAction("&Filter parameter ...",
                                             self.adjustFilterOptions,
                                             "Alt+F", QIcon(None),
                                             """Adjust filter parameters.""")
        selectPAction = self.createAction("&P", self.alterPhase, "Alt+P",
                                          QIcon(":/picon.png"),
                                          "Toggle P phase.", True)
        selectSAction = self.createAction("&S", self.alterPhase, "Alt+S",
                                          QIcon(":/sicon.png"),
                                          "Toggle S phase", True)
        printAction = self.createAction("&Print event ...",
                                        self.printEvent, QKeySequence.Print,
                                        QIcon(":/printer.png"),
                                        "Print waveform overview.")
        self.fileMenu = self.menuBar().addMenu('&File')
        self.fileMenuActions = (openEventAction, saveEventAction, None,
                                quitAction)
        self.fileMenu.aboutToShow.connect(self.updateFileMenu)

        self.editMenu = self.menuBar().addMenu('&Edit')
        for action in (filterAction, filterEditAction, None, selectPAction,
                       selectSAction, None, printAction):
            if action is None:
                self.editMenu.addSeparator()
            else:
                self.editMenu.addAction(action)

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 500)

        statsButtons = layoutStationButtons(self.getData(), self.getComponent())
        _layout.addLayout(statsButtons)
        _layout.addWidget(self.DataPlot)
        _widget.setLayout(_layout)
        self.setCentralWidget(_widget)

    def okToContinue(self):
        if self.dirty:
            return self.saveData()
        return True

    def plotData(self):
        pass #self.data.plotData(self.DataPlot)

    def filterData(self):
        pass

    def adjustFilterOptions(self):
        filterOptions = None
        fstring = "Filter Options ({0})".format(self.getSeismicPhase())
        filterDlg = FilterOptionsDialog(titleString=fstring,
                                        parent=self,
                                        filterOptions=self.getFilterOptions())
        if filterDlg.exec_():
            filterOptions = filterDlg.getFilterOptions()

        assert isinstance(filterOptions, FilterOptions)
        self.setFilterOptions(filterOptions)

    def getFilterOptions(self):
        return self.filteroptions

    def setFilterOptions(self, filterOptions):
        cases = {'P':self.filterOptionsP,
                 'S':self.filterOptionsS}
        cases[self.getSeismicPhase()] = filterOptions
        self.updateFilterOptions()

    def updateFilterOptions(self):
        try:
            self.filteroptions = [self.filterOptionsP
                                  if not self.seismicPhase == 'S'
                                  else self.filterOptionsS][0]
        except Exception, e:
            self.updateStatus('Error ...')
            emsg = QErrorMessage(self)
            emsg.showMessage('Error: {0}'.format(e))
        else:
            self.updateStatus('Filter loaded ...')

    def getSeismicPhase(self):
        return self.seismicPhase

    def alterPhase(self):
        pass

    def setSeismicPhase(self, phase):
        self.seismicPhase = self.seismicPhaseButtonGroup.getValue()

    def updateStatus(self, message):
        self.statusBar().showMessage(message, 5000)
        if self.getData() is not None:
            if not self.getData().isNew():
                self.setWindowTitle("PyLoT - processing event %s[*]" % self.getData().getID())
            elif self.getData().isNew():
                self.setWindowTitle("PyLoT - New event [*]")
            else:
                self.setWindowTitle("PyLoT - seismic processing the python way[*]")
        self.setWindowTitle("PyLoT - seismic processing the python way[*]")
        self.setWindowModified(self.dirty)

        self.statusBar().showMessage(message, 5000)

    def printEvent(self):
        pass

    def closeEvent(self, event):
        if self.okToContinue():
            self.closing.emit()
            QMainWindow.closeEvent(self, event)

    def helpHelp(self):
        if checkurl():
            form = HelpForm('https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/wiki')
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
