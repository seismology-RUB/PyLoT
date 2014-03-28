#!/usr/bin/env python
#
#
# Main program: QtPyLoT.py

import os
import platform
import sys
from PySide.QtCore import *
from PySide.QtGui import *
import helpform
from pylot.core.util import _getVersionString
from pylot.core.read.inputs import FilterOptions
from pylot.core.util import FILTERDEFAULTS
from pylot.core.util import MPLWidget

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
        dataLayout = setupPlot()

        filterDockWidget = FilterOptionsDock(titleString="Filter Options",
                                             parent=self,
                                             filterOptions=filteroptions)

        statLayout = self.layoutStationButtons(self.numStations)

        maingrid = QGridLayout()
        maingrid.setSpacing(10)
        maingrid.addLayout(statLayout, 0, 0)
        maingrid.addWidget()

    def setupPlot(self):
        # create a matplotlib widget
        self.DataPlot = MPLWidget()
        # create a layout inside the blank widget and add the matplotlib widget
        layout = QtGui.QVBoxLayout(self.ui.widget_PlotArea)
        layout.addWidget(self.DataPlot, 1)

        return layout

    def plotData(self, data):
        if data is not None and isinstance(data, Stream):
            self.stats.numStations = data.

    def updateFilterOptions(self):
        self.filteroptions = [self.filterOptionsP
                              if not self.seismicPhase == 'S'
                              else self.filterOptionsS]

    def layoutStationButtons(self, numStations):
        layout = QVBoxLayout()
        for n in range(numStations):
            stationButtons[n] = QPushButton('%s'.format(self.))

    def helpHelp(self):
        if internet_on():
            pass


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


main()
