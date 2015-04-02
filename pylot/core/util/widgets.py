# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:27:35 2014

@author: sebastianw
"""

import datetime
import numpy as np
import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvas
from matplotlib.widgets import MultiCursor
from PySide.QtGui import (QAction,
                          QApplication,
                          QComboBox,
                          QDateTimeEdit,
                          QDialog,
                          QDialogButtonBox,
                          QDoubleSpinBox,
                          QGroupBox,
                          QGridLayout,
                          QHBoxLayout,
                          QIcon,
                          QKeySequence,
                          QLabel,
                          QLineEdit,
                          QMessageBox,
                          QSpinBox,
                          QTabWidget,
                          QToolBar,
                          QVBoxLayout,
                          QWidget)
from PySide.QtCore import (QSettings,
                           Qt,
                           QUrl,
                           Signal,
                           Slot)
from PySide.QtWebKit import QWebView
from pylot.core.read import FilterOptions
from pylot.core.util.defaults import OUTPUTFORMATS
from pylot.core.util import prepTimeAxis, getGlobalTimes


def createAction(parent, text, slot=None, shortcut=None, icon=None,
                 tip=None, checkable=False):
    """
    :rtype : ~PySide.QtGui.QAction
    """
    action = QAction(text, parent)
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

class MPLWidget(FigureCanvas):

    def __init__(self, parent=None, xlabel='x', ylabel='y', title='Title'):

        self._parent = None
        self.setParent(parent)
        self.figure = Figure()
        # attribute plotdict is an dictionary connecting position and a name
        self.plotdict = dict()
        # create axes
        self.axes = self.figure.add_subplot(111)
        # clear axes each time plot is called
        self.axes.hold(True)
        # initialize super class
        FigureCanvas.__init__(self, self.figure)
        # add an cursor for station selection
        self.multiCursor = MultiCursor(self.figure.canvas, (self.axes,), horizOn=True,
                                       color='m', lw=1)
        # update labels of the entire widget
        self.updateWidget(xlabel, ylabel, title)

    def getPlotDict(self):
        return self.plotdict

    def setPlotDict(self, key, value):
        self.plotdict[key] = value

    def getParent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    def plotWFData(self, wfdata, title = None):
        self.axes.lines = []
        wfstart = getGlobalTimes(wfdata)[0]
        for n, trace in enumerate(wfdata):
            station = trace.stats.station
            print('plotting station: %s' % station)
            stime = trace.stats.starttime - wfstart
            time_ax = prepTimeAxis(stime, trace)
            trace.detrend()
            trace.detrend('demean')
            trace.normalize(trace.data.max() * 2)
            self.axes.plot(time_ax, trace.data + n, 'k')
            self.axes.hold(True)
            xlabel = 'seconds since {0}'.format(wfstart)
            ylabel = ''
            self.updateWidget(xlabel, ylabel, title)
            self.setPlotDict(n, station)
        self.axes.autoscale(tight=True)

    def updateXLabel(self, text):
        self.axes.set_xlabel(text)

    def updateYLabel(self, text):
        self.axes.set_ylabel(text)

    def updateTitle(self, text):
        self.axes.set_title(text)

    def updateWidget(self, xlabel, ylabel, title):
        self.updateXLabel(xlabel)
        self.updateYLabel(ylabel)
        self.updateTitle(title)

    def insertLabel(self, pos, text):
        pos = pos / max(self.axes.ylim)
        axann = self.axes.annotate(text, xy=(.03, pos), xycoords='axes fraction')
        axann.set_bbox(dict(facecolor='lightgrey', alpha=.6))


class multiComponentPlot(FigureCanvas):

    def __init__(self, data, parent=None, components='ZNE'):

        self.data = data
        self._parent = parent
        self.components = components

        self.figure = Figure()

        self.noc = len(components)
        FigureCanvas.__init__(self, self.figure)
        self.multiCursor = None
        self.resetPlot(components, data)

    def getData(self):
        return self.data

    def setData(self, data):
        self.data = data

    def getParent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    def getComponents(self):
        return self.components

    def setComponents(self, components):
        self.components = components

    def getNoC(self):
        return self.noc

    def setNoC(self, noc):
        self.noc = noc

    def resetPlot(self, components=None, data=None):

        # clear figure
        self.figure.clf()

        # delete multiCursor if existing
        if self.multiCursor is not None:
            self.multiCursor = None

        # set new attribute values
        if data is not None:
            self.setData(data)
        if components is not None:
            self.setComponents(components)
            noc = len(self.getComponents())
            if self.getNoC() != noc:
                self.setNoC(noc)
        self.axesdict = dict()

        # prepare variables for plotting
        stime = getGlobalTimes(self.getData())[0]

        xlabel = 'time since {0} [s]'.format(stime)

        # plot individual component traces in separate subplots
        for n, comp in enumerate(components):
            nsub = '{0}1{1}'.format(self.noc, n+1)
            if n >= 1:
                subax = self.figure.add_subplot(nsub, sharex=self.axesdict[0])
            else:
                subax = self.figure.add_subplot(nsub)
                subax.autoscale(tight=True)
            subset = data.copy().select(component=comp)[0]
            time_ax = prepTimeAxis(subset.stats.starttime - stime, subset)
            subax.plot(time_ax, subset.data)
            self.axesdict[n] = subax
            self.updateYLabel(n, comp)
            if n == self.noc:
                self.updateXLabel(self.noc, xlabel)
            else:
                self.updateXLabel(n, '')

        self.multiCursor = MultiCursor(self.figure.canvas,
                                       tuple(self.axesdict.values()),
                                       color='r', lw=1)

    def updateXLabel(self, pos, text):
        self.axesdict[pos].set_xlabel(text)

    def updateYLabel(self, pos, text):
        self.axesdict[pos].set_ylabel(text)


class PickDlg(QDialog):

    def __init__(self, parent=None, data=None, station=None, rotate=False):
        super(PickDlg, self).__init__(parent)

        # initialize attributes
        self.station = station
        self.rotate = rotate
        self.components = 'ZNE'

        # set attribute holding data
        if data is None:
            try:
                data = parent.getData().getWFData().copy()
                self.data = data.select(station=station)
            except AttributeError, e:
                errmsg = 'You either have to put in a data or an appropriate ' \
                         'parent (PyLoT MainWindow) object: {0}'.format(e)
                raise Exception(errmsg)
        else:
            self.data = data

        # initialize plotting widget
        self.multicompfig = MPLWidget(self)

        # setup ui
        self.setupUi()

        # plot data
        self.getPlotWidget().plotWFData(wfdata=self.getWFData(),
                                        title=self.getStation())

    def setupUi(self):

        # create icons
        filter_icon = QIcon()
        filter_icon.addPixmap(QPixmap(':/icons/filter.png'))

        # create actions
        self.filterAction = createAction(parent=self, text='Filter',
                                         slot=self.filterWFData,
                                         icon=filter_icon,
                                         tip='Filter waveforms',
                                         checkable=True)
        self.selectPhase = QComboBox()
        self.selectPhase.addItems(['Pn', 'Pg', 'P1', 'P2'])

        # layout the outermost appearance of the Pick Dialog
        _outerlayout = QVBoxLayout()
        _dialtoolbar = QToolBar()

        # fill toolbar with content

        _dialtoolbar.addAction(self.filterAction)
        _dialtoolbar.addWidget(self.selectPhase)

        _innerlayout = QHBoxLayout()

        _toolslayout = QVBoxLayout()
        _toolslabel = QLabel('Place for Tools')
        _toolslayout.addWidget(_toolslabel)

        _innerlayout.addLayout(_toolslayout)
        _innerlayout.addWidget(self.multicompfig)

        _outerlayout.addWidget(_dialtoolbar)
        _outerlayout.addLayout(_innerlayout)
        self.setLayout(_outerlayout)


    def getComponents(self):
        return self.components

    def getStation(self):
        return self.station

    def getPlotWidget(self):
        return self.multicompfig

    def getWFData(self):
        return self.data

    def filterWFData(self):
        data = self.getWFData().copy().filter(type='bandpass', freqmin=.5, freqmax=15.)
        title = self.getStation() + ' (filtered)'
        self.getPlotWidget().plotWFData(wfdata=data, title=title)


class PropertiesDlg(QDialog):

    def __init__(self, parent=None):
        super(PropertiesDlg, self).__init__(parent)

        appName = QApplication.applicationName()

        self.setWindowTitle("{0} Properties".format(appName))

        self.tabWidget = QTabWidget()
        self.tabWidget.addTab(InputsTab(self), "Inputs")
        self.tabWidget.addTab(OutputsTab(self), "Outputs")
        self.tabWidget.addTab(PhasesTab(self), "Phases")
        self.tabWidget.addTab(GraphicsTab(self), "Graphics")
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok |
                                          QDialogButtonBox.Apply |
                                          QDialogButtonBox.Close)

        layout = QVBoxLayout()
        layout.addWidget(self.tabWidget)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(self.apply)

    def accept(self, *args, **kwargs):
        self.apply()
        QDialog.accept(self)

    def apply(self):
        for widint in range(self.tabWidget.count()):
            curwid = self.tabWidget.widget(widint)
            values = curwid.getValues()
            if values is not None: self.setValues(values)

    def setValues(self, tabValues):
        settings = QSettings()
        for setting, value in tabValues.iteritems():
            settings.setValue(setting, value)
        settings.sync()


class PropTab(QWidget):

    def __init__(self, parent=None):
        super(PropTab, self).__init__(parent)

    def getValues(self):
        return None


class InputsTab(PropTab):

    def __init__(self, parent):
        super(InputsTab, self).__init__(parent)

        settings = QSettings()
        fulluser = settings.value("user/FullName")
        login = settings.value("user/Login")

        fullNameLabel = QLabel("Full name for user '{0}': ".format(login))

        # get the full name of the actual user
        self.fullNameEdit = QLineEdit()
        self.fullNameEdit.setText(fulluser)

        # information about data structure
        dataroot = settings.value("data/dataRoot")
        dataDirLabel = QLabel("data root directory: ")
        self.dataDirEdit = QLineEdit()
        self.dataDirEdit.setText(dataroot)
        self.dataDirEdit.selectAll()
        structureLabel = QLabel("data structure: ")
        self.structureSelect = QComboBox()

        from pylot.core.util.structure import DATASTRUCTURE

        self.structureSelect.addItems(DATASTRUCTURE.keys())

        layout = QGridLayout()
        layout.addWidget(dataDirLabel, 0, 0)
        layout.addWidget(self.dataDirEdit, 0, 1)
        layout.addWidget(fullNameLabel, 1, 0)
        layout.addWidget(self.fullNameEdit, 1, 1)
        layout.addWidget(structureLabel, 2, 0)
        layout.addWidget(self.structureSelect, 2, 1)

        self.setLayout(layout)

    def getValues(self):
        values = {}
        values["data/dataRoot"] = self.dataDirEdit.text()
        values["user/FullName"] = self.fullNameEdit.text()
        values["data/Structure"] = self.structureSelect.currentText()
        return values


class OutputsTab(PropTab):

    def __init__(self, parent=None):
        super(OutputsTab, self).__init__(parent)

        settings = QSettings()
        curval = settings.value("output/Format", None)

        eventOutputLabel = QLabel("event ouput format")
        self.eventOutputComboBox = QComboBox()
        eventoutputformats = OUTPUTFORMATS.keys()
        self.eventOutputComboBox.addItems(eventoutputformats)

        if curval is None:
            ind = 0
        else:
            ind = self.eventOutputComboBox.findText(curval)

        self.eventOutputComboBox.setCurrentIndex(ind)
        layout = QGridLayout()
        layout.addWidget(eventOutputLabel, 0, 0)
        layout.addWidget(self.eventOutputComboBox, 0, 1)

        self.setLayout(layout)

    def getValues(self):
        values = {}
        values["output/Format"] = self.eventOutputComboBox.currentText()
        return values

class PhasesTab(PropTab):

    def __init__(self, parent=None):
        super(PhasesTab, self).__init__(parent)

        pass


class GraphicsTab(PropTab):

    def __init__(self, parent=None):
        super(GraphicsTab, self).__init__(parent)

        pass


class NewEventDlg(QDialog):
    def __init__(self, parent=None, titleString="Create a new event"):
        """
        QDialog object utilized to create a new event manually.
        """
        super(NewEventDlg, self).__init__()

        self.setupUI()

        now = datetime.datetime.now()
        self.eventTimeEdit.setDateTime(now)
        # event dates in the future are forbidden
        self.eventTimeEdit.setMaximumDateTime(now)

        self.latEdit.setText("51.0000")
        self.lonEdit.setText("7.0000")
        self.depEdit.setText("10.0")

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

    def getValues(self):
        return {'origintime' : self.eventTimeEdit.dateTime().toPython(),
                'latitude' : self.latEdit.text(),
                'longitude' : self.lonEdit.text(),
                'depth' : self.depEdit.text()}

    def setupUI(self):

        # create widget objects
        timeLabel = QLabel()
        timeLabel.setText("Select time: ")
        self.eventTimeEdit = QDateTimeEdit()
        latLabel = QLabel()
        latLabel.setText("Latitude: ")
        self.latEdit = QLineEdit()
        lonLabel = QLabel()
        lonLabel.setText("Longitude: ")
        self.lonEdit = QLineEdit()
        depLabel = QLabel()
        depLabel.setText("Depth: ")
        self.depEdit = QLineEdit()

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok |
                                          QDialogButtonBox.Cancel)

        grid = QGridLayout()
        grid.addWidget(timeLabel, 0, 0)
        grid.addWidget(self.eventTimeEdit, 0, 1)
        grid.addWidget(latLabel, 1, 0)
        grid.addWidget(self.latEdit, 1, 1)
        grid.addWidget(lonLabel, 2, 0)
        grid.addWidget(self.lonEdit, 2, 1)
        grid.addWidget(depLabel, 3, 0)
        grid.addWidget(self.depEdit, 3, 1)
        grid.addWidget(self.buttonBox, 4, 1)

        self.setLayout(grid)

class FilterOptionsDialog(QDialog):

    def __init__(self, parent=None, titleString="Filter options",
                 filterOptions=None):
        """
        PyLoT widget FilterOptionsDialog is a QDialog object. It is an UI to
        adjust parameters for filtering seismic data.
        """
        super(FilterOptionsDialog, self).__init__()

        if parent is not None:
            self.filterOptions = parent.getFilterOptions()
        elif filterOptions is not None:
            self.filterOptions = FilterOptions(filterOptions)
        else:
            self.filterOptions = FilterOptions()

        _enable = True
        if self.getFilterOptions().getFilterType() is None:
            _enable = False

        self.freqminLabel = QLabel()
        self.freqminLabel.setText("minimum:")
        self.freqminSpinBox = QDoubleSpinBox()
        self.freqminSpinBox.setRange(5e-7, 1e6)
        self.freqminSpinBox.setDecimals(2)
        self.freqminSpinBox.setSuffix(' Hz')
        self.freqminSpinBox.setEnabled(_enable)

        self.freqmaxLabel = QLabel()
        self.freqmaxLabel.setText("maximum:")
        self.freqmaxSpinBox = QDoubleSpinBox()
        self.freqmaxSpinBox.setRange(5e-7, 1e6)
        self.freqmaxSpinBox.setDecimals(2)
        self.freqmaxSpinBox.setSuffix(' Hz')

        if _enable:
            self.freqminSpinBox.setValue(self.getFilterOptions().getFreq()[0])
            if self.getFilterOptions().getFilterType() in ['bandpass', 'bandstop']:
                self.freqmaxSpinBox.setValue(self.getFilterOptions().getFreq()[1])
        else:
            try:
                self.freqmaxSpinBox.setValue(self.getFilterOptions().getFreq())
                self.freqminSpinBox.setValue(self.getFilterOptions().getFreq())
            except TypeError, e:
                print e
                self.freqmaxSpinBox.setValue(1.)
                self.freqminSpinBox.setValue(.1)

        typeOptions = [None, "bandpass", "bandstop", "lowpass", "highpass"]

        self.orderLabel = QLabel()
        self.orderLabel.setText("Order:")
        self.orderSpinBox = QSpinBox()
        self.orderSpinBox.setRange(2, 10)
        self.orderSpinBox.setEnabled(_enable)
        self.selectTypeLabel = QLabel()
        self.selectTypeLabel.setText("Select filter type:")
        self.selectTypeCombo = QComboBox()
        self.selectTypeCombo.addItems(typeOptions)
        self.selectTypeLayout = QVBoxLayout()
        self.selectTypeLayout.addWidget(self.orderLabel)
        self.selectTypeLayout.addWidget(self.orderSpinBox)
        self.selectTypeLayout.addWidget(self.selectTypeLabel)
        self.selectTypeLayout.addWidget(self.selectTypeCombo)

        self.freqGroupBox = QGroupBox("Frequency range")
        self.freqGroupLayout = QGridLayout()
        self.freqGroupLayout.addWidget(self.freqminLabel, 0, 0)
        self.freqGroupLayout.addWidget(self.freqminSpinBox, 0, 1)
        self.freqGroupLayout.addWidget(self.freqmaxLabel, 1, 0)
        self.freqGroupLayout.addWidget(self.freqmaxSpinBox, 1, 1)
        self.freqGroupBox.setLayout(self.freqGroupLayout)

        self.freqmaxSpinBox.setEnabled(_enable)

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
                                          QDialogButtonBox.Cancel)

        grid = QGridLayout()
        grid.addWidget(self.freqGroupBox, 0, 2, 1, 2)
        grid.addLayout(self.selectTypeLayout, 1, 2, 1, 2)
        grid.addWidget(self.buttonBox, 2, 2, 1, 2)

        self.setLayout(grid)

        self.freqminSpinBox.valueChanged.connect(self.updateUi)
        self.freqmaxSpinBox.valueChanged.connect(self.updateUi)
        self.orderSpinBox.valueChanged.connect(self.updateUi)
        self.selectTypeCombo.currentIndexChanged.connect(self.updateUi)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)


    def updateUi(self):
        _enable = False
        if self.selectTypeCombo.currentText() not in ['bandpass', 'bandstop']:
            self.freqminLabel.setText("cutoff:")
            self.freqmaxSpinBox.setValue(self.freqminSpinBox.value())
        else:
            _enable = True
            self.freqminLabel.setText("minimum:")

        self.freqmaxLabel.setEnabled(_enable)
        self.freqmaxSpinBox.setEnabled(_enable)


        self.getFilterOptions().setFilterType(self.selectTypeCombo.currentText())
        freq = []
        freq.append(self.freqminSpinBox.value())
        if _enable:
            if self.freqminSpinBox.value() > self.freqmaxSpinBox.value():
                QMessageBox.warning(self, "Value error",
                                    "Maximum frequency must be at least the "
                                    "same value as minimum frequency (notch)!")
                self.freqmaxSpinBox.setValue(self.freqminSpinBox.value())
                self.freqmaxSpinBox.selectAll()
                self.freqmaxSpinBox.setFocus()
                return
            freq.append(self.freqmaxSpinBox.value())
        self.getFilterOptions().setFreq(freq)
        self.getFilterOptions().setOrder(self.orderSpinBox.value())

    def getFilterOptions(self):
        return self.filterOptions

    def accept(self):
        self.updateUi()
        QDialog.accept(self)


class LoadDataDlg(QDialog):

    def __init__(self, parent=None):
        super(LoadDataDlg, self).__init__(parent)

        pass


class HelpForm(QDialog):

    def __init__(self, page=QUrl('https://ariadne.geophysik.rub.de/trac/PyLoT'), parent=None):
        super(HelpForm, self).__init__(parent)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setAttribute(Qt.WA_GroupLeader)

        backAction = QAction(QIcon(":/back.png"), "&Back", self)
        backAction.setShortcut(QKeySequence.Back)
        homeAction = QAction(QIcon(":/home.png"), "&Home", self)
        homeAction.setShortcut("Home")
        self.pageLabel = QLabel()

        toolBar = QToolBar()
        toolBar.addAction(backAction)
        toolBar.addAction(homeAction)
        toolBar.addWidget(self.pageLabel)
        self.webBrowser = QWebView()
        self.webBrowser.load(page)

        layout = QVBoxLayout()
        layout.addWidget(toolBar)
        layout.addWidget(self.webBrowser, 1)
        self.setLayout(layout)

        self.connect(backAction, Signal("triggered()"),
                     self.webBrowser, Slot("backward()"))
        self.connect(homeAction, Signal("triggered()"),
                     self.webBrowser, Slot("home()"))
        self.connect(self.webBrowser, Signal("sourceChanged(QUrl)"),
                     self.updatePageTitle)

        self.resize(400, 600)
        self.setWindowTitle("{0} Help".format(QApplication.applicationName()))

    def updatePageTitle(self):
        self.pageLabel.setText(self.webBrowser.documentTitle())
