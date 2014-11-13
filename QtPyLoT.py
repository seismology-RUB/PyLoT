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

# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        # initialize filter parameter
        filterOptionsP = FILTERDEFAULTS['P']
        filterOptionsS = FILTERDEFAULTS['S']
        self.filterOptionsP = FilterOptions(**filterOptionsP)
        self.filterOptionsS = FilterOptions(**filterOptionsS)

        # initialize data
        self.data = None
        self.loadData()
        self.updateFilterOptions()
        self.startTime = min([tr.stats.starttime for tr in self.data])

        self.setupUi()

    def _getCurrentPlotType(self):
        return 'TestType'

    def loadData(self):
        self.data = Data()

    def getData(self):
        return self.data

    def getDataWidget(self):
        return self.DataPlot

    def setupUi(self):
        self.setWindowIcon(QIcon(":/pylot.ico"))

        xlab = self.startTime.strftime('seconds since %d %b %Y %H:%M:%S (%Z)')
        plottitle = self._getCurrentPlotType()

        # create central matplotlib figure widget
        self.DataPlot = MPLWidget(parent=self,
                                  xlabel=xlab,
                                  ylabel=None,
                                  title=plottitle)

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 5000)

        statLayout = layoutStationButtons(self.getData(), self.getComponent())
        dataLayout = self.getDataWidget()

        maingrid = QGridLayout()
        maingrid.setSpacing(10)
        maingrid.addLayout(statLayout, 0, 0)
        maingrid.addLayout(dataLayout, 1, 0)
        maingrid.setCentralWidget(dataLayout)

    def plotData(self):
        self.data.plotData(self.DataPlot)

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
                                  else self.filterOptionsS]
        except Exception, e:
            self.updateStatus('Error: %s' % e + ' ... no filteroptions loaded')
        else:
            self.updateStatus('Filteroptions succesfully loaded ...')

    def getSeismicPhase(self):
        return self.seismicPhase

    def setSeismicPhase(self, phase):
        self.seismicPhase = self.seismicPhaseButton.getValue()

    def updateStatus(self, message):
        self.statusBar().showMessage(message, 5000)

    def helpHelp(self):
        if checkurl():
            form = HelpForm('https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/wiki')
        else:
            form = HelpForm(':/help.html')
        form.show()


def main():
    # create the Qt application
    pylot_app = QApplication(sys.argv)

    # set Application Information
    pylot_app.setOrganizationName("Ruhr-University Bochum / MAGS2")
    pylot_app.setOrganizationDomain("rub.de")
    pylot_app.setApplicationName("PyLoT")
    pylot_app.setWindowIcon(QIcon(":/pylot.ico"))

    # create the main window
    pylot_form = MainWindow()

    # Show main window and run the app
    pylot_form.show()
    pylot_app.exec_()

if __name__ == "__main__":
    main()
