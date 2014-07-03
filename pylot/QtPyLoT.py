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
from pylot.core.util import _getVersionString
from pylot.core.read import FilterOptions
from pylot.core.util import FILTERDEFAULTS
from pylot.core.util import checkurl
from pylot.core.util import (PickDlg,
                             FilterOptionsDock,
                             PropertiesDlg,
                             MPLWidget,
                             HelpForm)

# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        filterOptionsP = FILTERDEFAULTS['P']
        filterOptionsS = FILTERDEFAULTS['S']
        self.filterOptionsP = FilterOptions(**filterOptionsP)
        self.filterOptionsS = FilterOptions(**filterOptionsS)

        self.loadData()
        self.updateFilterOptions()

        self.setupUi()

    def loadData(self):
        pass

    def setupUi(self):
        self.setWindowIcon(QIcon("PyLoT.ico"))

        # create central matplotlib figure widget
        self.DataPlot = MPLWidget(parent=self)

        filterDockWidget = FilterOptionsDock(titleString="Filter Options",
                                             parent=self,
                                             filterOptions=filteroptions)

        self.eventLabel = QLabel()
        self.eventLabel.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        status = self.statusBar()
        status.setSizeGripEnabled(False)
        status.addPermanentWidget(self.eventLabel)
        status.showMessage("Ready", 5000)

        statLayout = self.layoutStationButtons(self.numStations)

        maingrid = QGridLayout()
        maingrid.setSpacing(10)
        maingrid.addLayout(statLayout, 0, 0)
        maingrid.addLayout(dataLayout, 1, 0)
        maingrid.setCentralWidget(self.DataPlot)

    def plotData(self, data):
        if data is not None and isinstance(data, Stream):
            pass

    def updateFilterOptions(self):
        self.filteroptions = [self.filterOptionsP
                              if not self.seismicPhase == 'S'
                              else self.filterOptionsS]

    def updateStatus(self, message):
        self.statusBar().showMessage(message, 5000)

    def layoutStationButtons(self, numStations):
        layout = QVBoxLayout()
        for n in range(numStations):
            tr = data.select(component=self.dispOptions.comp)
            try:
                stationButtons[n] = QPushButton('%s'.format(
                                                tr[n].stats.station))
            except IndexError:
                error = QErrorMessage(self)
                errorString = QString()
                errorString.setText('''Number of stations does not match number
                                    of traces!''')
                error.showMessage(errorString)
                self.__del__()
        layout.addWidget(stationButtons)

        return layout

    def helpHelp(self):
        if checkurl():
            form = HelpForm('https://ariadne.geophysik.rub.de/trac/PyLoT/wiki')
        else:
            form = HelpForm(':/help.html')
        form.show()


def main():
    # create th Qt application
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
