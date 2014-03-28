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

    def setupUi(self):
        self.setWindowIcon(QIcon("PyLoT.ico"))

        # create central matplotlib figure widget
        dataLayout = setupPlot()

        filterDockWidget = FilterOptionsDock(titleString="Filter Options",
                                             parent=self,
                                             filterOptions=filteroptions)

        statLayout = self.layoutStationButtons()

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
            self.DataPlot.height

    def updateFilterOptions(self):
        self.filteroptions = [self.filterOptionsP if not self.seismicPhase == 'S'
                         else self.filterOptionsS]

if __name__ == '__main__':
    # Creating a Qt application
    pylot_app = QApplication(sys.argv)

    pylot_main = MainWindow()
    pylot_main.setWindowTitle('PyLoT-The Picking and Localization Tool')

    # Show main window and run the app
    pylot_main.show()
    pylot_app.exec_()
