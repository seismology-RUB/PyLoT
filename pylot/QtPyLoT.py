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

# Version information
__version__ = _getVersionString()


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        filterOptionsP = FILTERDEFAULTS['P']
        filterOptionsS = FILTERDEFAULTS['S']
        self.filterOptionsP = FilterOptions(**filterOptionsP)
        self.filterOptionsS = FilterOptions(**filterOptionsS)

        filteroptions = [self.filterOptionsP if not self.seismicPhase == 'S'
                         else self.filterOptionsS]
        filterDockWidget = FilterOptionsDock(titleString="Filter Options",
                                             parent=self,
                                             filterOptions=filteroptions)
        self.


class PickWindow(QDialog):

    def __init__(self, station=None, parent=None):
        super(PickWindow, self).__init__(parent)

        filterDockWidget = FilterOptionsDock(titleString="Filter Options",
                                             parent=self,
                                             filterOptions=filteroptions)


class PropertiesWindow(QDialog):

    def __init__(self, parent=None):
        super(PropertiesWindow, self).__init__(parent)


class FilterOptionsDock(QDockWidget):

    def __init__(self, parent=None, titleString="Filter options",
                 filterOptions=None):
        super(FilterOptionsDock, self).__init__()

        if filterOptions and not isinstance(filterOptions, FilterOptions):
            try:
                fOptions = FilterOptions(**filterOptions)
                filterOptions = fOptions
            except e:
                raise OptionsError('%s' % e)

        


class OptionsError(Exception):
    pass


if __name__ == '__main__':
    # Creating a Qt application
    pylot_app = QApplication(sys.argv)

    pylot_main = MainWindow()
    pylot_main.setWindowTitle('PyLoT-The Picking and Localization Tool')

    # Show main window and run the app
    pylot_main.show()
    pylot_app.exec_()
