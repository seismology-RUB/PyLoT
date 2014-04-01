# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:27:35 2014

@author: sebastianw
"""

import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from PySide.QtGui import (QDialog, QDockWidget)

from pylot.core.util import OptionsError
from pylot.core.util import FilterOptions


class MPLWidget(FigureCanvasQTAgg):

    def __init__(self, parent=None, xlabel='x', ylabel='y', title='Title'):
        super(MPLWidget, self).__init__(Figure())

        self.setParent(parent)
        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.add_subplot(111)

        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel(ylabel)
        self.axes.set_title(title)


class PickDlg(QDialog):

    def __init__(self, station=None, parent=None):
        super(PickDlg, self).__init__(parent)

        filterDockWidget = FilterOptionsDock(titleString="Filter Options",
                                             parent=self,
                                             filterOptions=filteroptions)


class PropertiesDlg(QDialog):

    def __init__(self, parent=None):
        super(PropertiesDlg, self).__init__(parent)


class FilterOptionsDock(QDockWidget):

    def __init__(self, parent=None, titleString="Filter options",
                 filterOptions=None):
        super(FilterOptionsDock, self).__init__()

        if filterOptions and not isinstance(filterOptions, FilterOptions):
            try:
                fOptions = FilterOptions(**filterOptions)
                filterOptions = fOptions
            except Exception, e:
                raise OptionsError('%s' % e)


class LoadDataDlg(QDialog):

    def __init__(self, parent=None):
        super(LoadDataDlg, self).__init__(parent)
