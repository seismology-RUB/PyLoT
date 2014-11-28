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
import platform
import sys
from PySide.QtCore import *
from PySide.QtGui import *
from obspy.core import (read, UTCDateTime)
from pylot import *
from pylot.core.util import _getVersionString
from pylot.core.read import (Data,
                             FilterOptions)
from pylot.core.util import FILTERDEFAULTS
from pylot.core.util import checkurl
from pylot.core.util import (PickDlg,
                             FilterOptionsDialog,
                             PropertiesDlg,
                             MPLWidget,
                             HelpForm)
from pylot.core.util import layoutStationButtons
import qrc_resources

# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        settings = QSettings()
        self.setWindowTitle("PyLoT - do seismic processing the pythonic way")
        self.setWindowIcon(QIcon(":/icon.ico"))
        self.seismicPhase = str(settings.value("phase", "P"))

        # initialize filter parameter
        filterOptionsP = FILTERDEFAULTS['P']
        filterOptionsS = FILTERDEFAULTS['S']
        print filterOptionsP, "\n", filterOptionsS
        self.filterOptionsP = FilterOptions(**filterOptionsP)
        self.filterOptionsS = FilterOptions(**filterOptionsS)

        # initialize data
        self.data = None
        self.loadData()
        self.updateFilterOptions()
        print self.filteroptions
        try:
            self.startTime = min([tr.stats.starttime for tr in self.data.wfdata])
        except:
            self.startTime = UTCDateTime()

        self.setupUi()

    def _getCurrentPlotType(self):
        return 'TestType'

    def createAction(self, text, slot=None, shortcut=None, icon=None,
                     tip=None, checkable=False):
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

    def createMenus(self):

        fileMenu = self.menuBar().addMenu("&File")
        fileMenu.addAction(self.openEventAction)
        fileMenu.addAction(self.saveEventAction)
        fileMenu.addAction(self.printAction)
        fileMenu.addSeparator()
        fileMenu.addAction(self.quitAction)
        
        editMenu = self.menuBar().addMenu("&Edit")
        editMenu.addAction(self.filterAction)
        editMenu.addAction(self.filterEditAction)        
        editMenu.addSeparator()
        editMenu.addAction(self.selectPAction)
        editMenu.addAction(self.selectSAction)

    def loadData(self):
        self.data = None

    def saveData(self):
        pass

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

        # create central matplotlib figure widget
        self.DataPlot = MPLWidget(parent=self,
                                  xlabel=xlab,
                                  ylabel=None,
                                  title=plottitle)
        
        self.setCentralWidget(self.getDataWidget())
        
        openIcon = self.style().standardIcon(QStyle.SP_DirOpenIcon)
        quitIcon = self.style().standardIcon(QStyle.SP_MediaStop)
        saveIcon = self.style().standardIcon(QStyle.SP_DriveHDIcon)
        self.openEventAction = self.createAction("&Open event ...",
                                                 self.loadData,
                                                 QKeySequence.Open,
                                                 openIcon,
                                                 "Open an event.")
        self.saveEventAction = self.createAction("&Save event ...",
                                                 self.saveData,
                                                 QKeySequence.Save, saveIcon,
                                                 "Save actual event data.")
        self.quitAction = self.createAction("&Quit", self.cleanUp,
                                            QKeySequence.Close,
                                            quitIcon,
                                            "Close event and quit PyLoT")
        self.filterAction = self.createAction("&Filter ...", self.filterData,
                                              "Ctrl+F", QIcon(":/filter.png"),
                                              """Toggle un-/filtered waveforms 
                                              to be displayed, according to the 
                                              desired seismic phase.""", True)
        self.filterEditAction = self.createAction("&Filter ...",
                                                  self.adjustFilterOptions,
                                                  "Alt+F", QIcon(None),
                                                  """Adjust filter
                                                  parameters.""")
        self.selectPAction = self.createAction("&P", self.alterPhase, "Alt+P",
                                               QIcon(":/picon.png"),
                                               "Toggle P phase.", True)
        self.selectSAction = self.createAction("&S", self.alterPhase, "Alt+S",
                                               QIcon(":/sicon.png"),
                                               "Toggle S phase", True)
        self.printAction = self.createAction("&Print event ...",
                                             self.printEvent,
                                             QKeySequence.Print,
                                             QIcon(":/printer.png"),
                                             "Print waveform overview.")
        self.createMenus()

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel|QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 10000)

        #statLayout = layoutStationButtons(self.getData(), self.getComponent())
        #dataLayout = self.getDataWidget()

#         maingrid = QGridLayout()
#         maingrid.setSpacing(10)
#         maingrid.addLayout(statLayout, 0, 0)
#         maingrid.addWidget(dataLayout, 1, 0)
        #self.setLayout(maingrid)

    def plotData(self):
        pass #self.data.plotData(self.DataPlot)

    def filterData(self):
        pass

    def adjustFilterOptions(self):
        fstring = "Filter Options ({0})".format(self.getSeismicPhase())
        filterDlg = FilterOptionsDialog(titleString=fstring,
                                        parent=self,
                                        filterOptions=self.getFilterOptions())
        filterDlg.accepted.connect(filterDlg.getFilterOptions)

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

    def printEvent(self):
        pass

    def cleanUp(self):
        pass

    def helpHelp(self):
        if checkurl():
            form = HelpForm('https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/wiki')
        else:
            form = HelpForm(':/help.html')
        form.show()


def main():
    # create the Qt application
    pylot_app = QApplication(['PyLoT'])

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
    main()
