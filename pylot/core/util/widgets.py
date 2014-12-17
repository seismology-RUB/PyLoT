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
from PySide.QtGui import (QAction,
                          QApplication,
                          QComboBox,
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
                           SIGNAL,
                           SLOT)
from PySide.QtWebKit import QWebView
from pylot.core.read import FilterOptions
from pylot.core.util.defaults import OUTPUTFORMATS


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

        pass


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

        self.connect(self.buttonBox, SIGNAL("accepted()"), self,
                     SLOT("accept()"))
        self.connect(self.buttonBox.button(QDialogButtonBox.Apply),
                     SIGNAL("clicked()"), self.apply)
        self.connect(self.buttonBox, SIGNAL("rejected()"),
                     self, SLOT("reject()"))

    def apply(self):
        settings = QSettings()
        for widint in range(self.tabWidget.count()):
            curwid = self.tabWidget.widget(widint)
            values = self.getValues(curwid)
            settings.setValue()


class InputsTab(QWidget):

    def __init__(self, parent):
        super(InputsTab, self).__init__(parent)

        settings = QSettings()
        fulluser = settings.value("user/FullName")
        login = settings.value("user/Login")

        fullNameLabel = QLabel("Full name for user '{0}'".format(login))

        parent.fullNameEdit = QLineEdit()
        parent.fullNameEdit.setText(fulluser)

        dataroot = settings.value("data/dataRoot")
        dataDirLabel = QLabel("data directory:")
        parent.dataDirEdit = QLineEdit()
        parent.dataDirEdit.setText(dataroot)
        parent.dataDirEdit.selectAll()

        layout = QGridLayout()
        layout.addWidget(dataDirLabel, 0, 0)
        layout.addWidget(parent.dataDirEdit, 0, 1)
        layout.addWidget(fullNameLabel, 1, 0)
        layout.addWidget(parent.fullNameEdit, 1, 1)

        self.setLayout(layout)


class OutputsTab(QWidget):

    def __init__(self, parent=None):
        super(OutputsTab, self).__init__(parent)

        eventOutputLabel = QLabel("event ouput format")
        eventOutputComboBox = QComboBox()
        eventoutputformats = OUTPUTFORMATS.keys()
        eventOutputComboBox.addItems(eventoutputformats)

        layout = QGridLayout()
        layout.addWidget(eventOutputLabel, 0, 0)
        layout.addWidget(eventOutputComboBox, 0, 1)

        self.setLayout(layout)


class PhasesTab(QWidget):

    def __init__(self, parent=None):
        super(PhasesTab, self).__init__(parent)

        pass


class GraphicsTab(QWidget):

    def __init__(self, parent=None):
        super(GraphicsTab, self).__init__(parent)

        pass


class FilterOptionsDialog(QDialog):

    def __init__(self, parent=None, titleString="Filter options",
                 filterOptions=None):
        """
        PyLoT widget FilterOptionsDialog is a QDialog object. It is an UI to
        adjust parameters for filtering seismic data.
        """
        super(FilterOptionsDialog, self).__init__()

        if filterOptions is not None and parent.getSeismicPhase() != "P":
            self.filterOptions = filterOptions
        else:
            self.filterOptions = FilterOptions()

        self.freqminLabel = QLabel()
        self.freqminLabel.setText("minimum:")
        self.freqminSpinBox = QDoubleSpinBox()
        self.freqminSpinBox.setRange(5e-7, 1e6)
        self.freqminSpinBox.setDecimals(2)
        self.freqminSpinBox.setSuffix(' Hz')
        self.freqminSpinBox.setValue(self.getFilterOptions().getFreq[0])
        self.freqmaxLabel = QLabel()
        self.freqmaxLabel.setText("maximum:")
        self.freqmaxSpinBox = QDoubleSpinBox()
        self.freqmaxSpinBox.setRange(5e-7, 1e6)
        self.freqmaxSpinBox.setDecimals(2)
        self.freqmaxSpinBox.setSuffix(' Hz')
        if self.filterOptions.filterType in ['bandpass', 'bandstop']:
            self.freqmaxSpinBox.setValue(self.getFilterOptions().getFreq[1])

        typeOptions = ["bandpass", "bandstop", "lowpass", "highpass"]

        self.orderLabel = QLabel()
        self.orderLabel.setText("Order:")
        self.orderSpinBox = QSpinBox()
        self.orderSpinBox.setRange(2, 10)
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
        if self.selectTypeCombo.currentText() not in ['bandpass', 'bandstop']:
            self.freqminLabel.setText("cutoff:")
            self.freqmaxLabel.setEnabled(False)
            self.freqmaxSpinBox.setEnabled(False)
            self.freqmaxSpinBox.setValue(self.freqminSpinBox.value())
        else:
            self.freqminLabel.setText("minimum:")
            self.freqmaxLabel.setEnabled(True)
            self.freqmaxSpinBox.setEnabled(True)

        self.filterOptions.filterType = self.selectTypeCombo.currentText()
        freq = []
        freq.append(self.freqminSpinBox.value())
        if self.filterOptions.filterType in ['bandpass', 'bandstop']:
            if self.freqminSpinBox.value() > self.freqmaxSpinBox.value():
                QMessageBox.warning(self, "Value error",
                                    "Maximum frequency must be at least the "
                                    "same value as minimum frequency (notch)!")
                self.freqmaxSpinBox.setValue(self.freqminSpinBox.value())
                self.freqmaxSpinBox.selectAll()
                self.freqmaxSpinBox.setFocus()
                return
            freq.append(self.freqmaxSpinBox.value())
        self.filterOptions.freq = freq
        self.filterOptions.order = self.orderSpinBox.value()
        
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

        self.connect(backAction, SIGNAL("triggered()"),
                     self.webBrowser, SLOT("backward()"))
        self.connect(homeAction, SIGNAL("triggered()"),
                     self.webBrowser, SLOT("home()"))
        self.connect(self.webBrowser, SIGNAL("sourceChanged(QUrl)"),
                     self.updatePageTitle)

        self.resize(400, 600)
        self.setWindowTitle("{0} Help".format(QApplication.applicationName()))

    def updatePageTitle(self):
        self.pageLabel.setText(self.webBrowser.documentTitle())
