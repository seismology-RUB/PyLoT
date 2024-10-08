# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:27:35 2014
@author: sebastianw
"""

import copy
import datetime
import getpass
import glob
import logging
import multiprocessing
import os
import subprocess
import sys
import time
import traceback

import matplotlib
import numpy as np
from pylot.core.io.phases import getQualitiesfromxml

matplotlib.use('QT5Agg')

from matplotlib.figure import Figure

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvas
except ImportError:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.widgets import MultiCursor
from obspy import read

from PySide2 import QtCore, QtGui, QtWidgets
from PySide2.QtGui import QIcon, QPixmap, QKeySequence
from PySide2.QtWidgets import QAction, QApplication, QCheckBox, QComboBox, \
    QDateTimeEdit, QDialog, QDialogButtonBox, QDoubleSpinBox, QGroupBox, \
    QGridLayout, QLabel, QLineEdit, QMessageBox, \
    QTabWidget, QToolBar, QVBoxLayout, QHBoxLayout, QWidget, \
    QPushButton, QFileDialog, QInputDialog
from PySide2.QtCore import QSettings, Qt, QUrl, Signal, QTimer
from PySide2.QtWebEngineWidgets import QWebEngineView as QWebView
from obspy import Stream, Trace, UTCDateTime
from obspy.core.util import AttribDict
from obspy.taup import TauPyModel
from obspy.taup.utils import get_phase_names
from pylot.core.io.data import Data
from pylot.core.io.inputs import FilterOptions, PylotParameter
from pylot.core.pick.utils import getSNR, earllatepicker, getnoisewin, \
    getResolutionWindow, get_quality_class
from pylot.core.pick.compare import Comparison
from pylot.core.pick.autopick import fmpicker
from pylot.core.util.defaults import OUTPUTFORMATS, FILTERDEFAULTS
from pylot.core.util.utils import prep_time_axis, full_range, demeanTrace, isSorted, findComboBoxIndex, clims, \
    pick_linestyle_plt, pick_color_plt, \
    check4rotated, check4doubled, check_for_gaps_and_merge, check_for_nan, identifyPhase, \
    loopIdentifyPhase, trim_station_components, transformFilteroptions2String, \
    identifyPhaseID, get_bool, get_none, pick_color, getAutoFilteroptions, SetChannelComponents, \
    station_id_remove_channel, get_pylot_eventfile_with_extension, get_possible_pylot_eventfile_extensions
from autoPyLoT import autoPyLoT
from pylot.core.util.thread import Thread
from pylot.core.util.dataprocessing import Metadata

if sys.version_info.major == 3:
    import icons_rc_3 as icons_rc
else:
    raise ImportError(f'Python version {sys.version_info.major} of current interpreter not supported.'
                      f'\nPlease use Python 3+.')

# workaround to prevent PyCharm from deleting icons_rc import when optimizing imports
icons_rc = icons_rc


class QSpinBox(QtWidgets.QSpinBox):
    ''' Custom SpinBox, insensitive to Mousewheel (prevents accidental changes when scrolling through parameters) '''

    def wheelEvent(self, event):
        event.ignore()


class QDoubleSpinBox(QtWidgets.QDoubleSpinBox):
    ''' Custom DoubleSpinBox, insensitive to Mousewheel (prevents accidental changes when scrolling through parameters) '''

    def wheelEvent(self, event):
        event.ignore()


class TextLogWidget(QtWidgets.QTextEdit):
    highlight = Signal()

    def __init__(self, parent, highlight_input=False):
        super(TextLogWidget, self).__init__(parent)
        self.highlight_input = highlight_input
        self.append('DUMMY TEXT\n')

    def write(self, text):
        self.append(text)
        if self.highlight_input:
            self.highlight.emit()


class LogWidget(QtWidgets.QWidget):
    def __init__(self, parent):
        super(LogWidget, self).__init__(parent, Qt.Window)
        self.current_active_error = False
        self.qmb = None
        self.setWindowTitle('PyLoT Log')
        self.setMinimumWidth(800)
        self.setMinimumHeight(600)

        self.stdout = TextLogWidget(self)
        self.stdout.write('DUMMY TEXT 2\n')
        self.stderr = TextLogWidget(self, highlight_input=True)
        self.stderr.write('DUMMY TEXT 2\n')
        self.stderr.highlight.connect(self.active_error)

        self.tabs = QTabWidget()
        self.tabs.addTab(self.stdout, 'Log')
        self.tabs.addTab(self.stderr, 'Errors')

        self.layout = QtWidgets.QVBoxLayout()
        self.textfield = QtWidgets.QTextEdit()
        self.setLayout(self.layout)
        self.layout.addWidget(self.tabs)

    def active_error(self):
        if not self.current_active_error:
            self.current_active_error = True
            self.show()
            self.activateWindow()
            self.tabs.setCurrentWidget(self.stderr)
            self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Warning,
                                             'Error', 'Error occurred. Please check log!')
            self.qmb.buttonClicked.connect(self.reset_error)
            self.qmb.show()

    def reset_error(self):
        # used to make sure that write errors is finished before raising new Message box etc.
        self.current_active_error = False
        self.stderr.append(60 * '#' + '\n\n')


def getDataType(parent):
    type = QInputDialog().getItem(parent, "Select phases type", "Type:",
                                  ["manual", "automatic"])

    if type[0].startswith('auto'):
        type = 'auto'
    else:
        type = type[0]

    return type


def plot_pdf(_axes, x, y, annotation, bbox_props, xlabel=None, ylabel=None,
             title=None):
    # try method or data
    try:
        _axes.plot(x, y())  # y provided as method
    except:
        _axes.plot(x, y)  # y provided as data

    if title:
        _axes.set_title(title)
    if xlabel:
        _axes.set_xlabel(xlabel)
    if ylabel:
        _axes.set_ylabel(ylabel)
    _anno = _axes.annotate(annotation, xy=(.05, .5), xycoords='axes fraction')
    _anno.set_bbox(bbox_props)
    _anno.draggable()

    return _axes


def createAction(parent, text, slot=None, shortcut=None, icon=None,
                 tip=None, checkable=False):
    """
    :rtype : ~PyQt5.QtWidgets.QAction
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


class ProgressBarWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(ProgressBarWidget, self).__init__(parent)
        self.hlayout = QtWidgets.QHBoxLayout()
        self.pb = QtWidgets.QProgressBar()
        self.pb.setRange(0, 0)
        self.label = QLabel()
        self.hlayout.addWidget(self.pb)
        self.hlayout.addWidget(self.label)
        self.setLayout(self.hlayout)
        self.hide()


class AddMetadataWidget(QWidget):
    def __init__(self, parent=None, metadata=None, windowflag=Qt.Window):
        super(AddMetadataWidget, self).__init__(parent, windowflag)
        self.inventories_add = []
        self.inventories_delete = []

        self.setWindowTitle("Manage inventory files")
        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)
        self.setupUI()
        self.connect_signals()
        self.resize(600, 450)

        self.metadata = metadata if metadata else Metadata(verbosity=0)
        self.from_metadata()

        self.center()
        self.show()
        # self.__test__()

    # TODO: what is this for? Remove?
    def __test__(self):
        self.add_item(r'/rscratch/minos14/marcel/git/pylot/tests')
        self.add_item(r'/rscratch/minos14/marcel/git/pylot/inputs')
        self.add_item(r'/rscratch/minos14/marcel/git/pylot/pylot')
        self.add_item(r'/rscratch/minos14/marcel/git/pylot/pylot_not_existing')

    def center(self):
        fm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        fm.moveCenter(centerPoint)
        self.move(fm.topLeft())

    def setupUI(self):
        self.init_lineedit()
        self.init_list_buttons_layout()
        self.init_add_remove_buttons()
        self.init_list_layout()
        self.init_accept_cancel_buttons()

    def init_lineedit(self):
        self.selection_layout = QVBoxLayout()
        self.lineedit_title = QtWidgets.QLabel("Choose metadata file to add:")
        self.selection_layout.insertWidget(0, self.lineedit_title)
        self.lineedit_layout = QHBoxLayout()
        self.selection_box = QtWidgets.QLineEdit()
        self.browse_button = self.init_button("Browse", explanation="Browse the file explorer")
        self.lineedit_layout.addWidget(self.selection_box)
        self.lineedit_layout.addWidget(self.browse_button)
        self.selection_layout.insertLayout(1, self.lineedit_layout)
        self.main_layout.insertLayout(0, self.selection_layout)

    def init_list_buttons_layout(self):
        self.list_buttons_layout = QHBoxLayout()
        self.main_layout.insertLayout(1, self.list_buttons_layout, 0)

    def init_button(self, label, key=None, explanation=None):
        """
        generates a button with custom name, label and shortcut
        :param label: name of the button (str)
        :param key: shortcut key (PySide.QtCore.Qt.Key)
        :param explanation: text that shows when hovering over the button
        :return: QPushButton
        """
        button = QPushButton(label)
        if key is not None:
            button.setShortcut(QtGui.QKeySequence(key))
        if explanation is not None:
            button.setToolTip(explanation)
        return button

    def init_add_remove_buttons(self):
        self.add_remove_layout = QVBoxLayout()
        self.add_button = self.init_button("Add", Qt.Key_Insert, "Choose a file to add \n(Ins)")
        self.remove_button = self.init_button("Remove", Qt.Key_Delete, "Remove selected file(s) \n(Del)")
        self.add_remove_layout.addWidget(self.add_button, 0, QtCore.Qt.AlignBottom)
        self.add_remove_layout.addWidget(self.remove_button, 0, QtCore.Qt.AlignTop)
        self.list_buttons_layout.insertLayout(1, self.add_remove_layout, 0)

    def init_list_layout(self):
        self.list_layout = QVBoxLayout()
        self.list_buttons_layout.insertLayout(0, self.list_layout, 1)
        self.init_list_title()
        self.init_list_widget()

    def init_list_title(self):
        self.list_title = QtWidgets.QLabel("Current metadata files:")
        self.list_layout.insertWidget(0, self.list_title)

    def init_list_widget(self):
        self.list_view = QtWidgets.QListView()
        self.list_view.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.list_model = QtGui.QStandardItemModel(self.list_view)
        self.list_view.setModel(self.list_model)
        self.list_layout.insertWidget(1, self.list_view, 1)
        self.sel_model = self.list_view.selectionModel()
        self.sel_model.selectionChanged.connect(self.on_clicked)

    def init_accept_cancel_buttons(self):
        self.accept_cancel_layout = QHBoxLayout()
        self.accept_button = self.init_button("Accept", Qt.Key_Return, "Accept changes and close \n(Return)")
        self.cancel_button = self.init_button("Cancel", Qt.Key_Escape, "Disregard changes and close \n(Esc)")
        self.accept_cancel_layout.addWidget(self.accept_button)
        self.accept_cancel_layout.addWidget(self.cancel_button)
        self.main_layout.insertLayout(3, self.accept_cancel_layout)

    def connect_signals(self):
        self.browse_button.clicked.connect(self.directory_to_lineedit)
        self.add_button.clicked.connect(self.add_item_from_lineedit)
        self.remove_button.clicked.connect(self.remove_item)
        self.accept_button.clicked.connect(self.accept)
        self.cancel_button.clicked.connect(self.hide)

    def from_metadata(self):
        """
        already existing metadata folders are added to the visual list of metadata folders
        """
        for inventory_path in self.metadata.inventories:
            item = QtGui.QStandardItem(inventory_path)
            self.list_model.appendRow(item)
            if inventory_path in self.metadata.obspy_dmt_invs:
                item.setFlags(Qt.NoItemFlags)

    def directory_to_lineedit(self):
        inventory_path = self.open_directory()
        self.selection_box.setText(inventory_path)

    def open_directory(self):
        """
        choosing a directory in the dialog window and returning the path
        :return: directory path
        """
        fninv = QFileDialog.getExistingDirectory(self, self.tr("Select inventory..."), self.tr("Select folder"))
        if not fninv:
            return
        return fninv

    def on_clicked(self):
        indices = self.list_view.selectionModel().selectedIndexes()
        try:
            item = self.list_model.itemData(indices[-1])
            inventory_path = item[0]
        except IndexError:
            inventory_path = ""
        self.selection_box.setText(inventory_path)

    def add_item_from_lineedit(self):
        """
        checks if the inventory path is already in the list and if it leads to a directory;
        adds the path to the visible list;
        adds the path to the list of inventories that get added to the metadata upon accepting all changes
        :param inventory_path: path of the folder that contains the metadata (unicode string)
        """
        inventory_path = self.selection_box.text()
        if not inventory_path:
            return
        if inventory_path in self.metadata.inventories or inventory_path in self.inventories_add:
            if not inventory_path in self.inventories_delete:
                QMessageBox.warning(self, 'Info', 'Path already in list!')
                return
        if not os.path.isdir(inventory_path):
            QMessageBox.warning(self, 'Warning', 'Path is no directory!')
            return
        item = QtGui.QStandardItem(inventory_path)
        item.setEditable(False)
        self.inventories_add.append(inventory_path)
        self.list_model.appendRow(item)  # adds path to visible list
        self.selection_box.setText("")

    def remove_item(self):
        """
        removes marked path from the visible list;
        also adds the path to the list of inventories that get deleted upon accepting all changes or deletes the
        path from the list of of inventories that get added to the metadata upon accepting all changes
        """
        for index in reversed(sorted(self.list_view.selectionModel().selectedIndexes())):
            item = self.list_model.itemData(index)
            inventory_path = item[0]  # marked path
            self.list_model.removeRow(index.row())  # aus der Anzeige-Liste gelöscht
            if inventory_path in self.inventories_add:
                self.inventories_add.remove(inventory_path)
            else:
                self.inventories_delete.append(inventory_path)

    def accept(self):
        """
        all inventory pathes in the add list get added to the metadata (class) and all inventory pathes in the delete
        list get deleted from the metadata (class); the lists are emptied and the window closes
        """
        for inventory_path in self.inventories_add:
            self.metadata.add_inventory(inventory_path)
        for inventory_path in self.inventories_delete:
            self.metadata.remove_inventory(inventory_path)
        self.clear_lists()
        self.hide()

    def clear_lists(self):
        self.inventories_add = []
        self.inventories_delete = []


class ComparisonWidget(QWidget):
    def __init__(self, c, parent=None, windowflag=Qt.Window):
        self._data = c
        self._stats = c.stations
        self._canvas = PlotWidget(self)
        self._widgets = dict(stationsComboBox=None,
                             phasesComboBox=None,
                             histCheckBox=None)
        self._phases = 'PS'
        self._plotprops = dict(station=list(self.stations)[0], phase=list(self.phases)[0])
        super(ComparisonWidget, self).__init__(parent, windowflag)
        self.setupUI()
        self.resize(1280, 720)
        self.plotcomparison()
        self.center()

    def center(self):
        fm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        fm.moveCenter(centerPoint)
        self.move(fm.topLeft())

    def setupUI(self):

        _outerlayout = QVBoxLayout(self)
        _innerlayout = QVBoxLayout()

        _stats_combobox = QComboBox(self)
        _stats_combobox.setObjectName('stationsComboBox')
        _stats_combobox.setEditable(True)
        _stats_combobox.setInsertPolicy(QComboBox.NoInsert)
        _stats_combobox.addItems(sorted(self.stations))
        _stats_combobox.editTextChanged.connect(self.prepareplot)
        self.widgets = _stats_combobox

        _phases_combobox = QComboBox(self)
        _phases_combobox.setObjectName('phasesComboBox')
        _phases_combobox.addItems(['P', 'S'])
        _phases_combobox.currentIndexChanged.connect(self.prepareplot)
        self.widgets = _phases_combobox

        self._hist_checkbox = QCheckBox('Show histograms', self)
        self._hist_checkbox.setObjectName('histCheckBox')
        self._hist_checkbox.stateChanged.connect(self.plothist)
        self.widgets = self._hist_checkbox

        self._toolbar = QToolBar(self)
        self._toolbar.addWidget(_stats_combobox)
        self._toolbar.addWidget(_phases_combobox)
        self._toolbar.addWidget(self._hist_checkbox)

        _innerlayout.addWidget(self.canvas)

        _outerlayout.addWidget(self._toolbar)
        _outerlayout.addLayout(_innerlayout)

        # finally layout the entire widget
        self.setLayout(_outerlayout)

    @property
    def canvas(self):
        return self._canvas

    @canvas.setter
    def canvas(self, canvas_obj):
        self._canvas = canvas_obj

    @property
    def stations(self):
        return self._stats

    @stations.setter
    def stations(self, stations):
        self._stats = stations

    @property
    def phases(self):
        return self._phases

    @phases.setter
    def phases(self, value):
        self._phases = value

    @property
    def plotprops(self):
        return self._plotprops

    @plotprops.setter
    def plotprops(self, values):
        try:
            key, value = values
            if key not in self.plotprops.keys():
                raise KeyError("'key' {0} not found in "
                               "ComparisonDialog.plotprops keys.".format(key))
        except ValueError:
            raise ValueError("Pass an iterable with two items")
        else:
            self._plotprops[key] = value

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        assert not isinstance(data, Comparison)
        self.stations = data.stations
        self._data = data

    @property
    def widgets(self):
        return self._widgets

    @widgets.setter
    def widgets(self, widget):
        name = widget.objectName()
        if name in self.widgets.keys():
            self._widgets[name] = widget

    def showToolbar(self):
        self._toolbar.show()

    def hideToolbar(self):
        self._toolbar.hide()

    def setHistboxChecked(self, bool):
        self._hist_checkbox.setChecked(bool)

    def clf(self):
        self.canvas.figure.clf()

    @staticmethod
    def hasvalue(sender):
        text = sender.currentText()
        index = sender.findText(text.upper())
        return index

    def prepareplot(self):
        try:
            _widget = self.sender()
            name = _widget.objectName()
            text = _widget.currentText().upper()
            index = self.hasvalue(_widget)
            if name == 'stationsComboBox' and index is not -1:
                _widget.setCurrentIndex(index)
                self.plotprops = ('station', text)
            elif name == 'phasesComboBox':
                self.plotprops = ('phase', text)
        except ValueError:
            raise ValueError('No sender widget given!')
        finally:
            self.plotcomparison()

    def plotcomparison(self):
        from matplotlib import gridspec

        _gs = gridspec.GridSpec(3, 2)
        self.clf()
        self.canvas.figure._tight = True
        _axes = self.canvas.figure.add_subplot(_gs[0:2, :])
        _ax1 = self.canvas.figure.add_subplot(_gs[2, 0])
        _ax2 = self.canvas.figure.add_subplot(_gs[2, 1])
        self.canvas.figure.tight_layout()

        # _axes.cla()
        station = self.plotprops['station']
        phase = self.plotprops['phase']
        if not phase in self.data.comparison[station]:
            _axes.set_title('No pick found for phase {}.'.format(phase))
            self.canvas.draw()
            return
        pdf = self.data.comparison[station][phase]
        x, y, std, exp = pdf.axis, pdf.data, pdf.standard_deviation(), \
                         pdf.expectation()

        annotation = "%s difference on %s\n" \
                     "expectation: %7.4f s\n" \
                     "std: %7.4f s" % (phase, station,
                                       exp, std)
        bbox_props = dict(boxstyle='round', facecolor='lightgrey', alpha=.7)

        plot_pdf(_axes, x, y, annotation, bbox_props, 'time difference [s]',
                 'propability density [-]', phase)

        pdf_a = copy.deepcopy(self.data.get('auto')[station][phase])
        pdf_m = copy.deepcopy(self.data.get('manu')[station][phase])

        xauto, yauto, stdauto, expauto, alim = pdf_a.axis, pdf_a.data(), \
                                               pdf_a.standard_deviation(), \
                                               pdf_a.expectation(), \
                                               pdf_a.limits()
        xmanu, ymanu, stdmanu, expmanu, mlim = pdf_m.axis, pdf_m.data(), \
                                               pdf_m.standard_deviation(), \
                                               pdf_m.expectation(), \
                                               pdf_m.limits()
        # find common limits
        lims = clims(alim, mlim)
        # relative x axis
        x0 = lims[0]
        xmanu -= x0
        xauto -= x0
        lims = [lim - x0 for lim in lims]
        x0 = UTCDateTime(x0)

        # set annotation text
        mannotation = "probability density of manual pick\n" \
                      "expectation: %7.4f s\n" \
                      "std: %7.4f s" % (expmanu - x0.timestamp, stdmanu)

        aannotation = "probability density of automatic pick\n" \
                      "expectation: %7.4f s\n" \
                      "std: %7.4f s" % (expauto - x0.timestamp, stdauto)

        _ax1 = plot_pdf(_ax1, xmanu, ymanu, mannotation,
                        bbox_props=bbox_props, xlabel='seconds since '
                                                      '{0}'.format(x0),
                        ylabel='probability density [-]')
        _ax1.set_xlim(lims)

        _ax2 = plot_pdf(_ax2, xauto, yauto, aannotation,
                        bbox_props=bbox_props, xlabel='seconds since '
                                                      '{0}'.format(x0))
        _ax2.set_xlim(lims)

        _gs.update(wspace=0.5, hspace=0.5)

        self.canvas.draw()

    def plothist(self):
        name = self.sender().objectName()
        if self.widgets[name].isChecked():
            for wname, widget in self.widgets.items():
                if wname != name:
                    self.widgets[wname].setEnabled(False)
            self.canvas.figure.clf()
            self.canvas.figure._tight = True
            _axPstd, _axPexp = self.canvas.figure.add_subplot(221), self.canvas.figure.add_subplot(223)
            _axSstd, _axSexp = self.canvas.figure.add_subplot(222), self.canvas.figure.add_subplot(224)
            axes_dict = dict(P=dict(std=_axPstd, exp=_axPexp),
                             S=dict(std=_axSstd, exp=_axSexp))
            bbox_props = dict(boxstyle='round', facecolor='lightgrey', alpha=.7)
            for phase in self.phases:
                std = self.data.get_std_array(phase)
                std = std[np.isfinite(std)]
                stdxlims = [0., 1.2 * max(std)]
                exp = self.data.get_expectation_array(phase)
                exp = exp[np.isfinite(exp)]
                eps_exp = 0.05 * (max(exp) - min(exp))
                expxlims = [min(exp) - eps_exp, max(exp) + eps_exp]
                axes_dict[phase]['std'].hist(std, range=stdxlims, bins=20, normed=False)
                axes_dict[phase]['exp'].hist(exp, range=expxlims, bins=20,
                                             normed=False)
                std_annotation = "Distribution curve for {phase} differences'\n" \
                                 "standard deviations (all stations)\n" \
                                 "number of samples: {nsamples}".format(phase=phase, nsamples=len(std))
                _anno_std = axes_dict[phase]['std'].annotate(std_annotation, xy=(.05, .8), xycoords='axes fraction')
                _anno_std.set_bbox(bbox_props)
                _anno_std.draggable()
                exp_annotation = "Distribution curve for {phase} differences'\n" \
                                 "expectations (all stations)\n" \
                                 "number of samples: {nsamples}".format(phase=phase, nsamples=len(exp))
                _anno_exp = axes_dict[phase]['exp'].annotate(exp_annotation, xy=(.05, .8), xycoords='axes fraction')
                _anno_exp.set_bbox(bbox_props)
                _anno_exp.draggable()
                axes_dict[phase]['exp'].set_xlabel('Time [s]')

                # add colors (early, late) for expectation
                ax = axes_dict[phase]['exp']
                xlims = ax.get_xlim()
                ylims = ax.get_ylim()
                # ax.fill_between([xlims[0], 0], ylims[0], ylims[1], color=(0.9, 1.0, 0.9, 0.5), label='earlier than manual')
                # ax.fill_between([0, xlims[1]], ylims[0], ylims[1], color=(1.0, 0.9, 0.9, 0.5), label='later than manual')
            legend = ax.legend()
            legend.draggable()

            for ax in axes_dict['P'].values():
                ax.set_ylabel('number of picks [-]')

            self.canvas.draw()
        else:
            for wname, widget in self.widgets.items():
                if wname != name:
                    self.widgets[wname].setEnabled(True)
            self.canvas.figure.clf()
            self.plotcomparison()


class PlotWidget(FigureCanvas):
    def __init__(self, parent=None, xlabel='x', ylabel='y', title='Title'):
        self._parent = parent
        self._fig = Figure()
        self._xl = xlabel
        self._yl = ylabel
        self._title = title
        super(PlotWidget, self).__init__(self.figure)

    @property
    def figure(self):
        return self._fig

    @figure.setter
    def figure(self, fig):
        self._fig = fig

    @property
    def xlabel(self):
        return self._xl

    @xlabel.setter
    def xlabel(self, label):
        self._xl = label

    @property
    def ylabel(self):
        return self._yl

    @ylabel.setter
    def ylabel(self, label):
        self._yl = label

    @property
    def title(self):
        return self._title

    @title.setter
    def title(self, title):
        self._title = title

    @property
    def parent(self):
        return self._parent


class WaveformWidgetPG(QtWidgets.QWidget):
    def __init__(self, parent, title='Title'):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.pg = self.parent().pg
        # added because adding widget to scrollArea will set scrollArea to parent
        self.orig_parent = parent
        # attribute plotdict is a dictionary connecting position and a name
        self.plotdict = dict()
        # init labels
        self.xlabel = None
        self.ylabel = None
        self.title = None
        # create plot
        self.main_layout = QtWidgets.QVBoxLayout()
        self.label_layout = QtWidgets.QHBoxLayout()
        self.add_labels()
        self.connect_signals()
        self.plotWidget = self.pg.PlotWidget(self.parent(), title=title)
        self.main_layout.addWidget(self.plotWidget)
        self.main_layout.addLayout(self.label_layout)
        self.init_labels()
        self.activateObspyDMToptions(False)
        self.plotWidget.showGrid(x=False, y=True, alpha=0.3)
        self.wfstart, self.wfend = 0, 0
        self.pen_multicursor = self.pg.mkPen(self.parent()._style['multicursor']['rgba'])
        self.pen_linecolor = self.pg.mkPen(self.parent()._style['linecolor']['rgba'])
        self.pen_linecolor_highlight = self.pg.mkPen((255, 100, 100, 255))
        self.pen_linecolor_syn = self.pg.mkPen((100, 0, 255, 255))
        self.reinitMoveProxy()
        self._proxy = self.pg.SignalProxy(self.plotWidget.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
        # self.plotWidget.getPlotItem().setDownsampling(auto=True)

    def reinitMoveProxy(self):
        self.vLine = self.pg.InfiniteLine(angle=90, movable=False, pen=self.pen_multicursor)
        self.hLine = self.pg.InfiniteLine(angle=0, movable=False, pen=self.pen_multicursor)
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)

    def mouseMoved(self, evt):
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.plotWidget.sceneBoundingRect().contains(pos):
            mousePoint = self.plotWidget.getPlotItem().vb.mapSceneToView(pos)
            x, y, = (mousePoint.x(), mousePoint.y())
            # if x > 0:# and index < len(data1):
            wfID = self.orig_parent.getWFID(y)
            station = self.orig_parent.getTraceID(wfID)
            abstime = self.wfstart + x
            if self.orig_parent.get_current_event():
                self.status_label.setText("station = {}, T = {}, t = {} [s]".format(station, abstime, x))
            self.vLine.setPos(mousePoint.x())
            self.hLine.setPos(mousePoint.y())

    def connect_signals(self):
        self.qcombo_processed.activated.connect(self.parent().newWF)
        self.comp_checkbox.clicked.connect(self.parent().newWF)

    def init_labels(self):
        self.label_layout.addWidget(self.status_label)
        for label in self.perm_labels:
            self.label_layout.addWidget(label)
        mid_widget = QWidget()
        right_widget = QWidget()
        # use widgets as placeholder, so that child widgets keep position when others are hidden
        mid_layout = QHBoxLayout()
        right_layout = QHBoxLayout()
        mid_layout.addWidget(self.comp_checkbox)
        right_layout.addWidget(self.qcombo_processed)
        mid_widget.setLayout(mid_layout)
        right_widget.setLayout(right_layout)
        self.label_layout.addWidget(mid_widget)
        self.label_layout.addWidget(right_widget)
        self.comp_checkbox.setLayoutDirection(Qt.RightToLeft)
        self.label_layout.setStretch(0, 4)
        self.label_layout.setStretch(1, 0)
        self.label_layout.setStretch(2, 0)
        self.label_layout.setStretch(3, 0)
        self.label_layout.setStretch(4, 3)
        self.label_layout.setStretch(5, 1)

    def add_labels(self):
        self.status_label = QtWidgets.QLabel()
        self.perm_labels = []
        for index in range(3):
            label = QtWidgets.QLabel()
            self.perm_labels.append(label)
        self.qcombo_processed = QtWidgets.QComboBox()
        self.comp_checkbox = QtWidgets.QCheckBox('Load comparison data')
        self.addQCboxItem('processed', 'green')
        self.addQCboxItem('raw', 'black')
        # self.perm_qcbox_right.setAlignment(2)
        self.setLayout(self.main_layout)

    def getPlotDict(self):
        return self.plotdict

    def activateObspyDMToptions(self, activate: bool) -> None:
        self.qcombo_processed.setEnabled(activate)

    def activateCompareOptions(self, activate: bool) -> None:
        self.comp_checkbox.setEnabled(activate)

    def setPermText(self, number, text=None, color='black'):
        if not 0 <= number < len(self.perm_labels):
            raise ValueError('No label for number {}'.format(number))
        self.perm_labels[number].setText(text)
        self.perm_labels[number].setStyleSheet('color: {}'.format(color))

    def addQCboxItem(self, text=None, color='black'):
        item = QtGui.QStandardItem(text)
        model = self.qcombo_processed.model()
        model.appendRow(item)
        item.setForeground(QtGui.QColor('{}'.format(color)))

    def setQCboxItem(self, text):
        index = self.qcombo_processed.findText(text)
        self.qcombo_processed.setCurrentIndex(index)

    def setPlotDict(self, key, value):
        self.plotdict[key] = value

    def clearPlotDict(self):
        self.plotdict = dict()

    def plotWFData(self, wfdata, wfsyn=None, title=None, scaleddata=False, mapping=True,
                   component='*', nth_sample=1, verbosity=0, method='normal', gain=1., shift_syn=0.2):
        def station_sort(nslc):
            """Try to sort after station integer in case of a line array (e.g. active seismics)"""
            try:
                rval = sorted(nslc, key=lambda x: int(x.split('.')[1]))
                return rval
            except ValueError as e:
                # this is the standard case for seismological stations
                pass
            except Exception as e:
                print(f'Sorting by station integer failed with unknown exception: {e}')

            # fallback to default sorting
            return sorted(nslc)

        if not wfdata:
            print('Nothing to plot.')
            return

        self.title = title
        self.clearPlotDict()
        self.wfstart, self.wfend = full_range(wfdata)
        nmax = 0

        settings = QSettings()
        compclass = SetChannelComponents.from_qsettings(settings)

        if not component == '*':
            alter_comp = compclass.getCompPosition(component)
            # alter_comp = str(alter_comp[0])

            st_select = wfdata.select(component=component)
            st_select += wfdata.select(component=alter_comp)
        else:
            st_select = wfdata

        # st_select, gaps = check_for_gaps_and_merge(st_select) #MP MP commented because probably done twice

        # list containing tuples of network, station, channel (for sorting)
        nslc = []
        for trace in st_select:
            nslc.append(
                trace.get_id())  # (trace.stats.network, trace.stats.station, trace.stats.location trace.stats.channel))
        nslc = station_sort(nslc)
        nslc.reverse()
        plots = []

        try:
            self.plotWidget.getPlotItem().vb.setLimits(xMin=float(0),
                                                       xMax=float(self.wfend - self.wfstart),
                                                       yMin=.5,
                                                       yMax=len(nslc) + .5)
        except:
            print('Warning: Could not set zoom limits')

        for n, seed_id in enumerate(nslc):
            n += 1
            network, station, location, channel = seed_id.split('.')
            st = st_select.select(id=seed_id)
            trace = st[0].copy()
            st_syn = wfsyn.select(id=seed_id)
            if st_syn:
                trace_syn = st_syn[0].copy()
            else:
                trace_syn = Trace()
            if mapping:
                comp = channel[-1]
                n = compclass.getPlotPosition(str(comp))
                # n = n[0]
            if n > nmax:
                nmax = n
            if verbosity:
                msg = 'plotting %s channel of station %s' % (channel, station)
                print(msg)
            stime = trace.stats.starttime - self.wfstart
            time_ax = prep_time_axis(stime, trace)
            if st_syn:
                stime_syn = trace_syn.stats.starttime - self.wfstart
                time_ax_syn = prep_time_axis(stime_syn, trace_syn)

            if method == 'fast':
                trace.data, time_ax = self.minMax(trace, time_ax)
                if trace_syn:
                    trace_syn.data, time_ax_syn = self.minMax(trace_syn, time_ax_syn)

            if len(time_ax) > 0:
                if not scaleddata:
                    trace.detrend('constant')
                    trace.normalize(np.max(np.abs(trace.data)) * 2)
                    if st_syn:
                        trace_syn.detrend('constant')
                        trace_syn.normalize(np.max(np.abs(trace_syn.data)) * 2)
                # TODO: change this to numpy operations instead of lists?
                times = np.array([time for index, time in enumerate(time_ax) if not index % nth_sample])
                times_syn = np.array(
                    [time for index, time in enumerate(time_ax_syn) if not index % nth_sample] if st_syn else [])
                trace.data = np.array(
                    [datum * gain + n for index, datum in enumerate(trace.data) if not index % nth_sample])
                trace_syn.data = np.array([datum + n + shift_syn for index, datum in enumerate(trace_syn.data)
                                           if not index % nth_sample] if st_syn else [])
                plots.append((times, trace.data,
                              times_syn, trace_syn.data))
                self.setPlotDict(n, seed_id)
        self.xlabel = 'seconds since {0}'.format(self.wfstart)
        self.ylabel = ''
        self.setXLims([0, self.wfend - self.wfstart])
        self.setYLims([0.5, nmax + 0.5])
        return plots

    def minMax(self, trace, time_ax):
        '''
        create min/max array for fast plotting (approach based on obspy __plot_min_max function)
        :returns data, time_ax
        '''
        npixel = self.orig_parent.width()
        ndata = len(trace.data)
        pts_per_pixel = ndata // npixel
        if pts_per_pixel < 2:
            return trace.data, time_ax
        remaining_samples = ndata % pts_per_pixel
        npixel = ndata // pts_per_pixel
        if remaining_samples:
            data = trace.data[:-remaining_samples]
        else:
            data = trace.data
        data = data.reshape(npixel, pts_per_pixel)
        min_ = data.min(axis=1)
        max_ = data.max(axis=1)
        if remaining_samples:
            extreme_values = np.empty((npixel + 1, 2), dtype=float)
            extreme_values[:-1, 0] = min_
            extreme_values[:-1, 1] = max_
            extreme_values[-1, 0] = \
                trace.data[-remaining_samples:].min()
            extreme_values[-1, 1] = \
                trace.data[-remaining_samples:].max()
        else:
            extreme_values = np.empty((npixel, 2), dtype=float)
            extreme_values[:, 0] = min_
            extreme_values[:, 1] = max_
        data = extreme_values.flatten()
        time_ax = np.linspace(time_ax[0], time_ax[-1], num=len(data))
        return data, time_ax

    # def getAxes(self):
    #     return self.axes

    # def getXLims(self):
    #     return self.getAxes().get_xlim()

    # def getYLims(self):
    #     return self.getAxes().get_ylim()

    def setXLims(self, lims):
        vb = self.plotWidget.getPlotItem().getViewBox()
        vb.setXRange(float(lims[0]), float(lims[1]), padding=0)

    def setYLims(self, lims):
        vb = self.plotWidget.getPlotItem().getViewBox()
        vb.setYRange(float(lims[0]), float(lims[1]), padding=0)

    def setYTickLabels(self, pos, labels):
        pos = list(pos)
        ticks = list(zip(pos, labels))
        minorTicks = [(0, 0) for _ in labels]
        # leftAx.tickLength = 5
        # leftAx.orientation = 'right'
        self.getAxItem('left').setTicks([ticks, minorTicks])

    def updateXLabel(self, text):
        self.getAxItem('bottom').setLabel(text)
        self.draw()

    def updateYLabel(self, text):
        self.getAxItem('left').setLabel(text)
        self.draw()

    def getAxItem(self, position):
        return self.plotWidget.getPlotItem().axes[position]['item']

    def updateTitle(self, text):
        self.plotWidget.getPlotItem().setTitle(text)
        self.draw()

    def updateWidget(self):  # , xlabel, ylabel, title):
        self.updateXLabel(self.xlabel)
        self.updateYLabel(self.ylabel)
        self.updateTitle(self.title)

    def draw(self):
        pass


class PylotCanvas(FigureCanvas):
    def __init__(self, figure=None, parent=None, connect_events=True, multicursor=False,
                 panZoomX=True, panZoomY=True):
        if not figure:
            figure = Figure()
            # create axes
            self.ax = figure.add_subplot(111)

        self.axes = figure.axes
        self.figure = figure
        if hasattr(parent, '_style'):
            self.style = parent._style
        else:
            self.style = None
        if self.style:
            self.figure.set_facecolor(self.style['background']['rgba_mpl'])
        # attribute plotdict is a dictionary connecting position and a name
        self.plotdict = dict()
        # initialize super class
        super(PylotCanvas, self).__init__(self.figure)
        self.setParent(parent)
        self.orig_parent = parent

        if multicursor:
            color_cursor = (0., 0., 0., 1.) if not self.style else self.style['multicursor']['rgba_mpl']
            # add a cursor for station selection
            self.multiCursor = MultiCursor(self.figure.canvas, self.axes,
                                           horizOn=True, useblit=True,
                                           color=color_cursor, lw=1)

        # initialize panning attributes
        self.press = None
        self.xpress = None
        self.ypress = None
        self.cur_xlim = None
        self.cur_ylim = None

        # panZoom activated selection
        self.panZoomX = panZoomX
        self.panZoomY = panZoomY

        self.limits = {}

        for ax in self.axes:
            self.limits[ax] = {'x': (-np.inf, np.inf),
                               'y': (-np.inf, np.inf)}

        if connect_events:
            self.connectEvents()

        try:
            self.figure.tight_layout()
        except:
            pass

        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

    def panPress(self, gui_event):
        ax_check = False
        for ax in self.axes:
            if gui_event.inaxes == ax:
                ax_check = True
                break
        if not ax_check: return
        self.cur_xlim = ax.get_xlim()
        self.cur_ylim = ax.get_ylim()
        self.press = gui_event.xdata, gui_event.ydata
        self.press_rel = gui_event.x, gui_event.y
        self.xpress, self.ypress = self.press

    def pan(self, gui_event):
        if self.press is None:
            return
        if gui_event.button == 1:
            self.panMotion(gui_event)
        elif gui_event.button == 3:
            if self.panZoomX or self.panZoomY:
                self.panZoom(gui_event)

    def panMotion(self, gui_event):
        ax_check = False
        for ax in self.axes:
            if gui_event.inaxes == ax:
                ax_check = True
                break
        if not ax_check: return
        dx = gui_event.xdata - self.xpress
        dy = gui_event.ydata - self.ypress
        self.cur_xlim -= dx
        self.cur_ylim -= dy
        ax.set_xlim(self.cur_xlim)
        ax.set_ylim(self.cur_ylim)
        self.refreshPickDlgText()
        ax.figure.canvas.draw_idle()

    def panRelease(self, gui_event):
        self.press = None
        self.press_rel = None
        self.figure.canvas.draw_idle()

    def panZoom(self, gui_event, threshold=2., factor=1.1):
        if not gui_event.x and not gui_event.y:
            return
        if not gui_event.button == 3:
            return
        ax_check = False
        for ax in self.axes:
            if gui_event.inaxes == ax:
                ax_check = True
                break
        if not ax_check: return

        # self.updateCurrentLimits() #maybe put this down to else:

        # calculate delta (relative values in axis)
        old_x, old_y = self.press_rel
        xdiff = gui_event.x - old_x
        ydiff = gui_event.y - old_y

        # threshold check
        if abs(xdiff) < threshold and abs(ydiff) < threshold:
            return

        # refresh press positions to new position
        self.press = gui_event.xdata, gui_event.ydata
        self.press_rel = gui_event.x, gui_event.y
        self.xpress, self.ypress = self.press

        if abs(xdiff) >= threshold and self.panZoomX:
            x_left, x_right = self.getXLims(ax)
            new_xlim = self.calcPanZoom(self.xpress, x_left, x_right, factor, (xdiff > 0))
            self.setXLims(ax, new_xlim)
        if abs(ydiff) >= threshold and self.panZoomY:
            y_bot, y_top = self.getYLims(ax)
            new_ylim = self.calcPanZoom(self.ypress, y_bot, y_top, factor, (ydiff > 0))
            self.setYLims(ax, new_ylim)

        self.refreshPickDlgText()
        self.draw()

    def set_frame_color(self, color='k'):
        for ax in self.axes:
            for spine in ax.spines.values():
                spine.set_edgecolor(color)

    def set_frame_linewidth(self, linewidth=1.):
        for ax in self.axes:
            for spine in ax.spines.values():
                spine.set_linewidth(linewidth)

    def saveFigure(self):
        if self.figure:
            fd = QtWidgets.QFileDialog()
            fname, filter = fd.getSaveFileName(self.parent(), filter='Images (*.png *.svg *.jpg)')
            if not fname:
                return
            if not any([fname.endswith(item) for item in ['.png', '.svg', '.jpg']]):
                fname += '.png'
            self.figure.savefig(fname)

    @staticmethod
    def calcPanZoom(origin, lower_b, upper_b, factor, positive):
        d_lower = abs(origin - lower_b)
        d_upper = abs(origin - upper_b)

        if positive:
            d_lower *= 1 - 1 / factor
            d_upper *= 1 - 1 / factor
            lower_b += d_lower
            upper_b -= d_upper
        else:
            d_lower /= 1 + 1 / factor
            d_upper /= 1 + 1 / factor
            lower_b -= d_lower
            upper_b += d_upper

        new_lim = [lower_b, upper_b]
        new_lim.sort()

        return new_lim

    def scrollZoom(self, gui_event, factor=2.):
        if not gui_event.xdata or not gui_event.ydata:
            return
        ax_check = False
        for ax in self.axes:
            if gui_event.inaxes == ax:
                ax_check = True
                break
        if not ax_check: return

        self.updateCurrentLimits()

        if gui_event.button == 'up':
            scale_factor = 1 / factor
        elif gui_event.button == 'down':
            # deal with zoom out
            scale_factor = factor
        else:
            # deal with something that should never happen
            scale_factor = 1
            print(gui_event.button)

        new_xlim = gui_event.xdata - \
                   scale_factor * (gui_event.xdata - self.getXLims(ax))
        new_ylim = gui_event.ydata - \
                   scale_factor * (gui_event.ydata - self.getYLims(ax))

        new_xlim.sort()
        global_x = self.getGlobalLimits(ax, 'x')
        global_y = self.getGlobalLimits(ax, 'y')
        new_xlim[0] = max(new_xlim[0], global_x[0])
        new_xlim[1] = min(new_xlim[1], global_x[1])
        new_ylim.sort()
        new_ylim[0] = max(new_ylim[0], global_y[0])
        new_ylim[1] = min(new_ylim[1], global_y[1])

        self.setXLims(ax, new_xlim)
        self.setYLims(ax, new_ylim)

        self.refreshPickDlgText()
        self.draw()

    def refreshPickDlgText(self):
        # TODO: Maybe decreasing performance if activated too often on move event
        # refresh text for pickdlg if given
        parent = self.parent()
        if hasattr(parent, 'refreshArrivalsText'):
            parent.refreshArrivalsText()
        if hasattr(parent, 'refreshPhaseText'):
            parent.refreshPhaseText()

    def keyPressHandler(self, gui_event):
        if gui_event.key == 'ctrl+p':
            self.saveFigure()

    def connectEvents(self):
        self.cidscroll = self.connectScrollEvent(self.scrollZoom)
        self.cidpress = self.connectPressEvent(self.panPress)
        self.cidmotion = self.connectMotionEvent(self.pan)
        self.cidrelease = self.connectReleaseEvent(self.panRelease)
        self.cidkpress = self.connectKeyPressEvent(self.keyPressHandler)

    def disconnectEvents(self):
        self.disconnectScrollEvent(self.cidscroll)
        self.disconnectMotionEvent(self.cidmotion)
        self.disconnectPressEvent(self.cidpress)
        self.disconnectReleaseEvent(self.cidrelease)
        self.disconnectKeyPressEvent(self.cidkpress)

        self.cidscroll = None
        self.cidrelease = None
        self.cidpress = None
        self.cidmotion = None
        self.cidkpress = None

    def disconnectPressEvent(self, cid):
        self.mpl_disconnect(cid)

    def connectPressEvent(self, slot):
        return self.mpl_connect('button_press_event', slot)

    def disconnectMotionEvent(self, cid):
        self.mpl_disconnect(cid)

    def connectMotionEvent(self, slot):
        return self.mpl_connect('motion_notify_event', slot)

    def disconnectReleaseEvent(self, cid):
        self.mpl_disconnect(cid)

    def connectReleaseEvent(self, slot):
        return self.mpl_connect('button_release_event', slot)

    def disconnectScrollEvent(self, cid):
        self.mpl_disconnect(cid)

    def connectScrollEvent(self, slot):
        return self.mpl_connect('scroll_event', slot)

    def disconnectKeyPressEvent(self, cid):
        self.mpl_disconnect(cid)

    def connectKeyPressEvent(self, slot):
        return self.mpl_connect('key_press_event', slot)

    def getPlotDict(self):
        return self.plotdict

    def setPlotDict(self, key, value):
        self.plotdict[key] = value

    def clearPlotDict(self):
        self.plotdict = dict()

    @staticmethod
    def calcPlotPositions(wfdata, compclass):
        possible_plot_pos = list(range(len(wfdata)))
        plot_positions = {}
        for trace in wfdata:
            comp = trace.stats.channel[-1]
            try:
                position = compclass.getPlotPosition(str(comp))
            except ValueError as e:
                continue
            plot_positions[trace.stats.channel] = position
        for channel, plot_pos in plot_positions.items():
            while not plot_pos in possible_plot_pos or not plot_pos - 1 in plot_positions.values():
                if plot_pos == 0:
                    break
                plot_pos -= 1
                if plot_pos < 0:
                    raise Exception('Plot position lower zero. This should not happen.')
            plot_positions[channel] = plot_pos
        return plot_positions

    def plotWFData(self, wfdata, wfdata_compare=None, title=None, zoomx=None, zoomy=None,
                   noiselevel=None, scaleddata=False, mapping=True,
                   component='*', nth_sample=1, iniPick=None, verbosity=0,
                   plot_additional=False, additional_channel=None, scaleToChannel=None,
                   snr=None):
        def get_wf_dict(data: Stream = Stream(), linecolor = 'k', offset: float = 0., **plot_kwargs):
            return dict(data=data, linecolor=linecolor, offset=offset, plot_kwargs=plot_kwargs)


        ax = self.axes[0]
        ax.cla()

        self.clearPlotDict()
        wfstart, wfend = full_range(wfdata)
        nmax = 0

        settings = QSettings()
        compclass = SetChannelComponents.from_qsettings(settings)

        linecolor = (0., 0., 0., 1.) if not self.style else self.style['linecolor']['rgba_mpl']

        plot_streams = dict(wfdata=get_wf_dict(linecolor=linecolor, linewidth=0.7),
                            wfdata_comp=get_wf_dict(offset=0.1, linecolor='b', alpha=0.7, linewidth=0.5))

        if not component == '*':
            alter_comp = compclass.getCompPosition(component)
            # alter_comp = str(alter_comp[0])

            plot_streams['wfdata']['data'] = wfdata.select(component=component)
            plot_streams['wfdata']['data'] += wfdata.select(component=alter_comp)
            if wfdata_compare:
                plot_streams['wfdata_comp']['data'] = wfdata_compare.select(component=component)
                plot_streams['wfdata_comp']['data'] += wfdata_compare.select(component=alter_comp)
        else:
            plot_streams['wfdata']['data'] = wfdata
            if wfdata_compare:
                plot_streams['wfdata_comp']['data'] = wfdata_compare

        st_main = plot_streams['wfdata']['data']

        if mapping:
            plot_positions = self.calcPlotPositions(st_main, compclass)

        # list containing tuples of network, station, channel and plot position (for sorting)
        nslc = []
        for plot_pos, trace in enumerate(st_main):
            if not trace.stats.channel[-1] in ['Z', 'N', 'E', '1', '2', '3']:
                print('Warning: Unrecognized channel {}'.format(trace.stats.channel))
                continue
            nslc.append(trace.get_id())
        nslc.sort()
        nslc.reverse()

        for n, seed_id in enumerate(nslc):
            network, station, location, channel = seed_id.split('.')
            for wf_name, wf_dict in plot_streams.items():
                st_select = wf_dict.get('data')
                if not st_select:
                    continue
                st = st_select.select(id=seed_id)
                trace = st[0].copy()
                if mapping:
                    n = plot_positions[trace.stats.channel]
                if n > nmax:
                    nmax = n
                if verbosity:
                    msg = 'plotting %s channel of station %s' % (channel, station)
                    print(msg)
                stime = trace.stats.starttime - wfstart
                time_ax = prep_time_axis(stime, trace)
                if time_ax is not None:
                    if scaleToChannel:
                        st_scale = wfdata.select(channel=scaleToChannel)
                        if st_scale:
                            tr = st_scale[0]
                            trace.detrend('constant')
                            trace.normalize(np.max(np.abs(tr.data)) * 2)
                            scaleddata = True
                    if not scaleddata:
                        trace.detrend('constant')
                        trace.normalize(np.max(np.abs(trace.data)) * 2)

                    offset = wf_dict.get('offset')

                    times = [time for index, time in enumerate(time_ax) if not index % nth_sample]
                    data = [datum + n + offset for index, datum in enumerate(trace.data) if not index % nth_sample]
                    ax.axhline(n, color="0.5", lw=0.5)
                    ax.plot(times, data, color=wf_dict.get('linecolor'), **wf_dict.get('plot_kwargs'))
                    if noiselevel is not None:
                        for level in [-noiselevel[channel], noiselevel[channel]]:
                            ax.plot([time_ax[0], time_ax[-1]],
                                    [n + level, n + level],
                                    color=wf_dict.get('linecolor'),
                                    linestyle='dashed')
                    self.setPlotDict(n, seed_id)
        if plot_additional and additional_channel:
            compare_stream = wfdata.select(channel=additional_channel)
            if compare_stream:
                trace = compare_stream[0]
                if scaleToChannel:
                    st_scale = wfdata.select(channel=scaleToChannel)
                    if st_scale:
                        tr = st_scale[0]
                        trace.detrend('constant')
                        trace.normalize(np.max(np.abs(tr.data)) * 2)
                        scaleddata = True
                if not scaleddata:
                    trace.detrend('constant')
                    trace.normalize(np.max(np.abs(trace.data)) * 2)
                time_ax = prep_time_axis(stime, trace)
                times = [time for index, time in enumerate(time_ax) if not index % nth_sample]
                p_data = compare_stream[0].data
                # #normalize
                # p_max = max(abs(p_data))
                # p_data /= p_max
                for index in range(3):
                    ax.plot(times, p_data, color='red', alpha=0.5, linewidth=0.7)
                    p_data += 1

        if iniPick:
            ax.vlines(iniPick, ax.get_ylim()[0], ax.get_ylim()[1],
                      colors='m', linestyles='dashed',
                      linewidth=2)
        xlabel = 'seconds since {0}'.format(wfstart)
        ylabel = ''
        self.updateWidget(xlabel, ylabel, title)
        self.setXLims(ax, [0, wfend - wfstart])
        self.setYLims(ax, [-0.5, nmax + 0.5])
        if zoomx is not None:
            self.setXLims(ax, zoomx)
        if zoomy is not None:
            self.setYLims(ax, zoomy)
        if snr is not None:
            if snr < 2:
                warning = 'LOW SNR'
                if snr < 1.5:
                    warning = 'VERY LOW SNR'
                ax.text(0.1, 0.9, 'WARNING - {}'.format(warning), ha='center', va='center', transform=ax.transAxes,
                        color='red')

        self.draw()

    @staticmethod
    def getXLims(ax):
        return ax.get_xlim()

    @staticmethod
    def getYLims(ax):
        return ax.get_ylim()

    @staticmethod
    def setXLims(ax, lims):
        ax.set_xlim(lims)

    @staticmethod
    def setYLims(ax, lims):
        ax.set_ylim(lims)

    def setYTickLabels(self, pos, labels):
        self.axes[0].set_yticks(list(pos))
        self.axes[0].set_yticklabels(labels)
        self.draw()

    def updateXLabel(self, text):
        self.axes[0].set_xlabel(text)
        self.draw()

    def updateYLabel(self, text):
        self.axes[0].set_ylabel(text)
        self.draw()

    def updateTitle(self, text):
        self.axes[0].set_title(text, verticalalignment='bottom')
        self.draw()

    def updateWidget(self, xlabel, ylabel, title):
        self.updateXLabel(xlabel)
        self.updateYLabel(ylabel)
        self.updateTitle(title)

    def insertLabel(self, pos, text):
        pos = pos / max(self.axes[0].ylim)
        axann = self.axes[0].annotate(text, xy=(.03, pos),
                                      xycoords='axes fraction')
        axann.set_bbox(dict(facecolor='lightgrey', alpha=.6))

    def setZoomBorders2content(self):
        if not self.axes:
            return
        for ax in self.limits.keys():
            xlims = self.getXLims(ax)
            ylims = self.getYLims(ax)

            self.limits[ax] = {'x': xlims,
                               'y': ylims}

            for axis, limit in self.limits[ax].items():
                self.setGlobalLimits(ax, axis, limit)

    def updateCurrentLimits(self):
        for ax in self.limits.keys():
            self.setXLims(ax, self.getXLims(ax))
            self.setYLims(ax, self.getYLims(ax))

    def getGlobalLimits(self, ax, axis):
        return self.limits[ax][axis]

    def setGlobalLimits(self, ax, axis, limits):
        self.limits[ax][axis] = limits

    def resetZoom(self):
        for ax in self.figure.axes:
            self.setXLims(ax, self.getGlobalLimits(ax, 'x'))
            self.setYLims(ax, self.getGlobalLimits(ax, 'y'))
        self.draw()


class SearchFileByExtensionDialog(QtWidgets.QDialog):
    def __init__(self, parent=None, label='Text: ', default_text='.xml', events=None):
        super(SearchFileByExtensionDialog, self).__init__(parent)
        self.events = events
        self.filepaths = []
        self.file_extensions = []
        self.check_all_state = True
        self.merge_strategy = None
        self.default_text = default_text
        self.label = label
        self.setButtons()
        self.setupUi()
        self.connectSignals()
        self.showPaths()
        self.refreshSelectionBox()
        # self.refresh_timer = QTimer(self)
        # self.refresh_timer.timeout.connect(self.showPaths)
        # self.refresh_timer.start(10000)

        self.resize(800, 450)

    def setupUi(self):
        ncol = 4
        self.main_layout = QtWidgets.QVBoxLayout()
        self.header_layout = QtWidgets.QHBoxLayout()
        self.footer_layout = QtWidgets.QHBoxLayout()
        #
        self.setLayout(self.main_layout)

        # widgets inside the dialog
        self.textLabel = QtWidgets.QLabel(self.label)
        self.comboBox = QtWidgets.QComboBox()
        self.comboBox.addItem(self.default_text)
        self.comboBox.setEditable(True)

        # optional search button, currently disabled. List refreshed when text changes
        self.searchButton = QtWidgets.QPushButton('Search')
        self.searchButton.setVisible(False)

        # check/uncheck button for table
        self.checkAllButton = QtWidgets.QPushButton('Check/Uncheck all')

        # radiobutton for merge selection
        self.mergeRadioButtonGroup = QtWidgets.QButtonGroup()
        self.merge_button = QtWidgets.QRadioButton('Merge')
        self.overwrite_button = QtWidgets.QRadioButton('Overwrite')
        self.mergeRadioButtonGroup.addButton(self.merge_button)
        self.mergeRadioButtonGroup.addButton(self.overwrite_button)
        self.merge_button.setChecked(True)
        self.merge_strategy = self.merge_button.text()

        # table
        self.tableWidget = QtWidgets.QTableWidget()
        tableWidget = self.tableWidget
        tableWidget.setColumnCount(ncol)
        tableWidget.setRowCount(len(self.events))
        tableWidget.setHorizontalHeaderLabels(('', 'Event ID', 'Filename', 'Last modified'))
        tableWidget.setEditTriggers(tableWidget.NoEditTriggers)
        tableWidget.setSortingEnabled(True)
        header = tableWidget.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)

        self.statusText = QtWidgets.QLabel()

        self.header_layout.addWidget(self.textLabel)
        self.header_layout.addWidget(self.comboBox)
        self.header_layout.addWidget(self.searchButton)

        self.footer_layout.addWidget(self.checkAllButton)
        self.footer_layout.addWidget(self.statusText)
        self.footer_layout.addWidget(self.merge_button)
        self.footer_layout.addWidget(self.overwrite_button)

        self.footer_layout.setStretch(0, 0)
        self.footer_layout.setStretch(1, 1)

        self.main_layout.addLayout(self.header_layout)
        self.main_layout.addWidget(self.tableWidget)
        self.main_layout.addLayout(self.footer_layout)
        self.main_layout.addWidget(self._buttonbox)

    def showPaths(self):
        self.filepaths = []
        fext = self.comboBox.currentText()
        self.tableWidget.clearContents()
        for index, event in enumerate(self.events):
            filename = get_pylot_eventfile_with_extension(event, fext)
            pf_selected_item = QtWidgets.QTableWidgetItem()
            check_state = QtCore.Qt.Checked if filename else QtCore.Qt.Unchecked
            pf_selected_item.setCheckState(check_state)
            self.tableWidget.setItem(index, 0, pf_selected_item)
            self.tableWidget.setItem(index, 1, QtWidgets.QTableWidgetItem(f'{event.pylot_id}'))
            if filename:
                self.filepaths.append(filename)
                ts = int(os.path.getmtime(filename))

                # create QTableWidgetItems of filepath and last modification time
                fname_item = QtWidgets.QTableWidgetItem(f'{os.path.split(filename)[-1]}')
                fname_item.setData(3, filename)
                ts_item = QtWidgets.QTableWidgetItem(f'{datetime.datetime.fromtimestamp(ts)}')
                self.tableWidget.setItem(index, 2, fname_item)
                self.tableWidget.setItem(index, 3, ts_item)

        self.update_status()

    def refreshSelectionBox(self):
        fext = self.comboBox.currentText()
        self.file_extensions = [fext]

        for event in self.events:
            extensions = get_possible_pylot_eventfile_extensions(event, '*.xml')
            for ext in extensions:
                if not ext in self.file_extensions:
                    self.file_extensions.append(ext)

        self.comboBox.clear()
        for ext in sorted(self.file_extensions):
            self.comboBox.addItem(ext)

    def setButtons(self):
        self._buttonbox = QDialogButtonBox(QDialogButtonBox.Ok |
                                           QDialogButtonBox.Cancel)

    def toggleCheckAll(self):
        self.check_all_state = not self.check_all_state
        self.checkAll(self.check_all_state)

    def checkAll(self, state):
        state = QtCore.Qt.Checked if state else QtCore.Qt.Unchecked
        for row_ind in range(self.tableWidget.rowCount()):
            item = self.tableWidget.item(row_ind, 0)
            item.setCheckState(state)

    def getChecked(self):
        filepaths = []
        for row_ind in range(self.tableWidget.rowCount()):
            item_check = self.tableWidget.item(row_ind, 0)
            if item_check.checkState() == QtCore.Qt.Checked:
                item_fname = self.tableWidget.item(row_ind, 2)
                if item_fname:
                    filepath = item_fname.data(3)
                    filepaths.append(filepath)
        return filepaths

    def update_status(self, row=None, col=None):
        if col is not None and col != 0:
            return
        filepaths = self.getChecked()
        if len(filepaths) > 0:
            status_text = f"Found {len(filepaths)} eventfile{'s' if len(filepaths) > 1 else ''}. Do you want to load them?"
        else:
            status_text = 'Did not find/select any files for specified file mask.'
        self.statusText.setText(status_text)

    def update_merge_strategy(self):
        for button in (self.merge_button, self.overwrite_button):
            if button.isChecked():
                self.merge_strategy = button.text()

    def connectSignals(self):
        self._buttonbox.accepted.connect(self.accept)
        self._buttonbox.rejected.connect(self.reject)
        self.comboBox.editTextChanged.connect(self.showPaths)
        self.searchButton.clicked.connect(self.showPaths)
        self.checkAllButton.clicked.connect(self.toggleCheckAll)
        self.checkAllButton.clicked.connect(self.update_status)
        self.tableWidget.cellClicked.connect(self.update_status)
        self.merge_button.clicked.connect(self.update_merge_strategy)
        self.overwrite_button.clicked.connect(self.update_merge_strategy)


class SingleTextLineDialog(QtWidgets.QDialog):
    def __init__(self, parent=None, label='Text: ', default_text='.xml'):
        super(SingleTextLineDialog, self).__init__(parent)
        self.default_text = default_text
        self.label = label
        self.setButtons()
        self.connectSignals()
        self.setupUi()

    def setupUi(self):
        self.main_layout = QtWidgets.QVBoxLayout()
        self.sub_layout = QtWidgets.QHBoxLayout()
        #
        self.setLayout(self.main_layout)
        self.textLabel = QtWidgets.QLabel(self.label)
        self.lineEdit = QtWidgets.QLineEdit(self.default_text)

        self.sub_layout.addWidget(self.textLabel)
        self.sub_layout.addWidget(self.lineEdit)

        self.main_layout.addLayout(self.sub_layout)
        self.main_layout.addWidget(self._buttonbox)

    def setButtons(self):
        self._buttonbox = QDialogButtonBox(QDialogButtonBox.Ok |
                                           QDialogButtonBox.Cancel)

    def connectSignals(self):
        self._buttonbox.accepted.connect(self.accept)
        self._buttonbox.rejected.connect(self.reject)


class PhaseDefaults(QtWidgets.QDialog):
    def __init__(self, parent=None, nrow=10,
                 phase_defaults=['P', 'S'],
                 current_phases=[]):
        super(PhaseDefaults, self).__init__(parent)
        self.nrow = nrow
        self.checktoggle = True
        self.main_layout = QtWidgets.QVBoxLayout()
        self.sub_layout = QtWidgets.QGridLayout()
        self.phase_names = phase_defaults
        self.current_phases = current_phases
        self.setButtons()
        self.setupUi()
        self.connectSignals()
        self.setWindowTitle('Default Phases')
        self.selected_phases = []

    def setButtons(self):
        self._check_all_button = QtWidgets.QPushButton('Check/Uncheck all')
        self._buttonbox = QDialogButtonBox(QDialogButtonBox.Ok |
                                           QDialogButtonBox.Cancel)

    def setupUi(self):
        self.setLayout(self.main_layout)
        self.main_layout.addWidget(self._check_all_button)
        self.main_layout.addLayout(self.sub_layout)
        self.init_phase_layout()
        self.main_layout.addWidget(self._buttonbox)

    def connectSignals(self):
        self._check_all_button.clicked.connect(self.toggleAllChecked)
        self._buttonbox.accepted.connect(self.accept)
        self._buttonbox.rejected.connect(self.reject)

    def toggleAllChecked(self):
        for box in self._checkboxes.values():
            box.setChecked(self.checktoggle)
        self.checktoggle = not self.checktoggle

    def init_phase_layout(self):
        self._checkboxes = {}
        row = 0
        column = 0
        for index, phase in enumerate(self.phase_names):
            if row > self.nrow:
                column += 1
                row = 0
            checkbox = self.create_phase_box(phase)
            self.sub_layout.addWidget(checkbox,
                                      row, column)
            self._checkboxes[phase] = checkbox
            checkbox.setChecked(bool(phase in self.current_phases))
            row += 1

    def create_phase_box(self, phase_name):
        checkbox = QtWidgets.QCheckBox(phase_name)
        return checkbox

    def update_selected_phases(self):
        self.selected_phases = []
        for phase in self.phase_names:
            if self._checkboxes[phase].isChecked():
                self.selected_phases.append(phase)

    def accept(self):
        self.update_selected_phases()
        QtWidgets.QDialog.accept(self)


class PickDlg(QDialog):
    update_picks = QtCore.Signal(dict)

    def __init__(self, parent=None, data=None, data_compare=None, station=None, network=None, location=None, picks=None,
                 autopicks=None, rotate=False, parameter=None, embedded=False, metadata=None, show_comp_data=False,
                 event=None, filteroptions=None, wftype=None):
        super(PickDlg, self).__init__(parent, Qt.Window)
        self.orig_parent = parent
        self.setAttribute(Qt.WA_DeleteOnClose)

        # initialize attributes
        self.parameter = parameter
        model = self.parameter.get('taup_model')
        self._embedded = embedded
        self.showCompData = show_comp_data
        self.station = station
        self.network = network
        self.location = location
        self.rotate = rotate
        self.metadata = metadata
        self.wftype = wftype
        self.pylot_event = event
        self.components = 'ZNE'
        self.currentPhase = None
        self.phaseText = []
        self.phaseLines = {}
        self.arrivals = []
        self.arrivalsText = []
        self.cidpick = []
        self.cidpress = None
        self.removed_picks = []
        settings = QSettings()
        pylot_user = getpass.getuser()
        self._user = settings.value('user/Login', pylot_user)
        self._dirty = False
        self._style = parent._style
        if picks:
            self.picks = copy.deepcopy(picks)
            self._init_picks = picks
        else:
            self.picks = {}
            self._init_picks = {}
        if autopicks:
            self.autopicks = copy.deepcopy(autopicks)
            self._init_autopicks = autopicks
        else:
            self.autopicks = {}
            self._init_autopicks = {}
        if filteroptions:
            self.filteroptions = filteroptions
        else:
            self.filteroptions = FILTERDEFAULTS
        self.pick_block = False

        # set attribute holding data
        if data is None or not data:
            try:
                data = parent.get_data().getWFData().copy()
                self.data = data.select(station=station)
            except AttributeError as e:
                errmsg = 'You either have to put in a data or an appropriate ' \
                         'parent (PyLoT MainWindow) object: {0}'.format(e)
                raise Exception(errmsg)
        else:
            self.data = data
        self.data_compare = data_compare

        self.nextStation = QtWidgets.QCheckBox('Continue with next station ')

        # comparison channel
        self.referenceChannel = QtWidgets.QComboBox()
        self.referenceChannel.activated.connect(self.resetPlot)

        # comparison channel
        self.compareCB = QtWidgets.QCheckBox()
        self.compareCB.setChecked(self.showCompData)
        self.compareCB.clicked.connect(self.switchCompData)
        self.compareCB.clicked.connect(self.resetPlot)
        self.compareCB.setVisible(bool(self.data_compare))

        # scale channel
        self.scaleChannel = QtWidgets.QComboBox()
        self.scaleChannel.activated.connect(self.resetPlot)

        # initialize panning attributes
        self.press = None
        self.xpress = None
        self.ypress = None
        self.cur_xlim = None
        self.cur_ylim = None

        self.stime, self.etime = full_range(self.getWFData())

        # initialize plotting widget
        self.multicompfig = PylotCanvas(parent=self, multicursor=True)
        self.phaseplot = PhasePlotWidget(self)
        self.phaseplot.hide()

        # setup ui
        self.setupUi()

        # fill compare and scale channels
        self.referenceChannel.addItem('-', None)
        self.scaleChannel.addItem('individual', None)

        for trace in self.getWFData():
            channel = trace.stats.channel
            self.referenceChannel.addItem(channel, trace)
            if not channel[-1] in ['Z', 'N', 'E', '1', '2', '3']:
                print('Skipping unknown channel for scaling: {}'.format(channel))
                continue
            self.scaleChannel.addItem(channel, trace)
            actionP = self.pChannels.addAction(str(channel))
            actionS = self.sChannels.addAction(str(channel))
            actionP.setCheckable(True)
            actionS.setCheckable(True)
            actionP.setChecked(self.getChannelSettingsP(channel))
            actionS.setChecked(self.getChannelSettingsS(channel))

        # plot data
        title = self.getStation()
        if self.wftype is not None:
            title += ' | ({})'.format(self.wftype)

        self.multicompfig.plotWFData(wfdata=self.getWFData(), wfdata_compare=self.getWFDataComp(),
                                     title=title)

        self.multicompfig.setZoomBorders2content()

        self.multicompfig.updateCurrentLimits()
        self.multicompfig.draw()
        self.multicompfig.setFocus()

        # set plot labels
        self.setPlotLabels()

        # draw picks if present
        self.drawAllPicks()

        # init expected picks using obspy Taup
        try:
            if self.metadata and model != "None":
                self.model = TauPyModel(model)
                self.get_arrivals()
                self.drawArrivals()
                self.activateArrivalsButton(True)
            else:
                self.activateArrivalsButton(False)
        except Exception as e:
            print('Warning: Could not init expected picks from taup: {}'.format(e))
            print(traceback.format_exc())
            self.activateArrivalsButton(False)

        # init pick delete (with middle mouse click)
        self.connect_pick_delete()
        self.connect_mouse_motion()
        self.setWindowTitle('Pickwindow on station: {}'.format(self.getStation()))
        self.setWindowState(QtCore.Qt.WindowMaximized)

    def setupUi(self):
        menuBar = QtWidgets.QMenuBar(self)
        if not self._embedded:
            exitMenu = menuBar.addMenu('File')
            exitAction = QtWidgets.QAction('Close', self)
            exitAction.triggered.connect(self.close)
            exitMenu.addAction(exitAction)

        # create matplotlib toolbar to inherit functionality
        self.figToolBar = NavigationToolbar2QT(self.multicompfig, self)
        self.figToolBar.hide()

        # create icons
        filter_icon_p = QIcon()
        filter_icon_p.addPixmap(QPixmap(':/icons/filter_p.png'))
        filter_icon_s = QIcon()
        filter_icon_s.addPixmap(QPixmap(':/icons/filter_s.png'))
        key_a_icon = QIcon()
        key_a_icon.addPixmap(QPixmap(':/icons/key_A.png'))
        zoom_icon = QIcon()
        zoom_icon.addPixmap(QPixmap(':/icons/zoom_in.png'))
        home_icon = QIcon()
        home_icon.addPixmap(QPixmap(':/icons/zoom_0.png'))
        del_icon = QIcon()
        del_icon.addPixmap(QPixmap(':/icons/delete.png'))
        sync_icon = QIcon()
        sync_icon.addPixmap(QPixmap(':/icons/sync.png'))

        # create actions
        self.filterActionP = createAction(parent=self, text='Apply P Filter',
                                          slot=self.filterP,
                                          icon=filter_icon_p,
                                          tip='Toggle filtered/original'
                                              ' waveforms',
                                          checkable=True,
                                          shortcut='P')
        self.filterActionS = createAction(parent=self, text='Apply S Filter',
                                          slot=self.filterS,
                                          icon=filter_icon_s,
                                          tip='Toggle filtered/original'
                                              ' waveforms',
                                          checkable=True,
                                          shortcut='S')
        self.autoFilterAction = createAction(parent=self, text='Automatic Filtering',
                                             slot=self.toggleAutoFilter,
                                             icon=key_a_icon,
                                             tip='Filter automatically before initial pick',
                                             checkable=True,
                                             shortcut='Ctrl+A')
        self.zoomAction = createAction(parent=self, text='Zoom',
                                       slot=self.zoom, icon=zoom_icon,
                                       tip='Zoom into waveform',
                                       checkable=True)
        self.resetZoomAction = createAction(parent=self, text='Home',
                                            slot=self.resetPlot, icon=home_icon,
                                            tip='Reset zoom to original limits')
        self.resetPicksAction = createAction(parent=self, text='Delete Picks',
                                             slot=self.delPicks, icon=del_icon,
                                             tip='Delete current picks.')
        self.renamePhaseAction = createAction(parent=self, text='Rename Phase',
                                              slot=self.initRenamePhase, icon=sync_icon,
                                              tip='Rename a Phase.', checkable=True,
                                              shortcut='R')

        self.addPickPhases(menuBar)

        self.pChannels = menuBar.addMenu('P-Channels')
        self.sChannels = menuBar.addMenu('S-Channels')
        self.pChannels.triggered.connect(self.updateChannelSettingsP)
        self.sChannels.triggered.connect(self.updateChannelSettingsS)

        settings = QSettings()
        self.autoFilterAction.setChecked(get_bool(settings.value('autoFilter')))

        # create other widget elements
        phaseitems = [None] + list(FILTERDEFAULTS.keys())

        # create buttons for P and S filter and picking
        self.p_button = QPushButton('P', self)
        self.s_button = QPushButton('S', self)
        self.p_button.setMinimumWidth(100)
        self.s_button.setMinimumWidth(100)
        self.p_button.setCheckable(True)
        self.s_button.setCheckable(True)
        # set button tooltips
        # self.p_button.setToolTip('Hotkey: "1"')
        # self.s_button.setToolTip('Hotkey: "2"')

        self.plot_arrivals_button = QPushButton('Plot phases')
        self.plot_arrivals_button.setCheckable(True)

        # create accept/reject button
        self.accept_button = QPushButton('&Accept')
        self.reject_button = QPushButton('&Reject')
        self.disable_ar_buttons()

        self.statusbar = QtWidgets.QStatusBar(self)

        # add hotkeys
        self._shortcut_space = QtWidgets.QShortcut(QtGui.QKeySequence(' '), self)
        self._shortcut_space.activated.connect(self.accept_button.clicked)
        # button shortcuts (1 for P-button, 2 for S-button)
        # self.p_button.setShortcut(QKeySequence('1'))
        # self.s_button.setShortcut(QKeySequence('2'))

        # layout the outermost appearance of the Pick Dialog
        _outerlayout = QVBoxLayout()
        _dialtoolbar = QToolBar()
        _dialtoolbar.setStyleSheet('QToolBar{spacing:5px;}')

        # fill toolbar with content
        _dialtoolbar.addAction(self.filterActionP)
        _dialtoolbar.addAction(self.filterActionS)
        _dialtoolbar.addAction(self.autoFilterAction)
        _dialtoolbar.addWidget(self.p_button)
        _dialtoolbar.addWidget(self.s_button)
        _dialtoolbar.addAction(self.zoomAction)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addAction(self.resetZoomAction)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addAction(self.resetPicksAction)
        _dialtoolbar.addAction(self.renamePhaseAction)
        _dialtoolbar.addSeparator()
        if self._embedded:
            manu_label = QLabel('Manual Onsets:')
            manu_label.setStyleSheet('QLabel {'
                                     'padding:2px;'
                                     'padding-left:5px}')
            _dialtoolbar.addWidget(manu_label)
            _dialtoolbar.addWidget(self.accept_button)
            _dialtoolbar.addWidget(self.reject_button)
        else:
            _dialtoolbar.addWidget(self.nextStation)
            _dialtoolbar.addSeparator()
        est_label = QLabel('Estimated onsets:')
        est_label.setStyleSheet('QLabel {'
                                'padding:2px;'
                                'padding-left:5px}')
        _dialtoolbar.addWidget(est_label)
        _dialtoolbar.addWidget(self.plot_arrivals_button)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addWidget(QtWidgets.QLabel('Plot reference channel: '))
        _dialtoolbar.addWidget(self.referenceChannel)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addWidget(QtWidgets.QLabel('Compare: '))
        _dialtoolbar.addWidget(self.compareCB)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addWidget(QtWidgets.QLabel('Scaling: '))
        _dialtoolbar.addWidget(self.scaleChannel)

        # layout the innermost widget
        _innerlayout = QVBoxLayout()
        _innerinnerlayout = QtWidgets.QHBoxLayout()
        _lowerlayout = QHBoxLayout()
        _innerinnerlayout.addWidget(self.multicompfig)
        _innerinnerlayout.addWidget(self.phaseplot)
        _innerlayout.addLayout(_innerinnerlayout)
        _innerlayout.addLayout(_lowerlayout)
        _lowerlayout.addWidget(self.statusbar)

        # add button box to the dialog
        _buttonbox = QDialogButtonBox(QDialogButtonBox.Ok |
                                      QDialogButtonBox.Cancel)

        # merge widgets and layouts to establish the dialog
        if not self._embedded:
            _lowerlayout.addWidget(_buttonbox)
        _outerlayout.addWidget(menuBar)
        _outerlayout.addWidget(_dialtoolbar)
        _outerlayout.addLayout(_innerlayout)
        _outerlayout.setStretch(0, 0)
        _outerlayout.setStretch(1, 0)
        _outerlayout.setStretch(2, 1)
        _lowerlayout.setStretch(0, 5)
        _lowerlayout.setStretch(1, 1)
        _innerlayout.setStretch(0, 1)
        _innerlayout.setStretch(1, 0)

        # connect widget element signals with slots (methods to the dialog
        # object
        self.p_button.clicked.connect(self.p_clicked)
        self.s_button.clicked.connect(self.s_clicked)
        self.accept_button.clicked.connect(self.accept)
        self.reject_button.clicked.connect(self.reject)
        self.accept_button.clicked.connect(self.disable_ar_buttons)
        self.reject_button.clicked.connect(self.disable_ar_buttons)
        self.plot_arrivals_button.clicked.connect(self.toggle_arrivals_plot)
        _buttonbox.accepted.connect(self.accept)
        _buttonbox.rejected.connect(self.reject)

        # finally layout the entire dialog
        self.setLayout(_outerlayout)
        self.resize(1280, 720)

    def activateArrivalsButton(self, val=True):
        self.plot_arrivals_button.setEnabled(val)

    def toggle_arrivals_plot(self):
        if self.plot_arrivals_button.isChecked():
            self.plot_arrivals()
        else:
            self.hide_arrivals_plot()

    def hide_arrivals_plot(self):
        self.phaseplot.hide()

    def plot_arrivals(self):
        if self.phaseplot.new:
            self.get_arrivals(True)
            ax = self.phaseplot.ax
            self.arrivals.plot(ax=ax, show=False)
            ax.legend(loc=1)
            self.phaseplot.new = False
            self.phaseplot.draw()
        self.phaseplot.show()

    def setDirty(self, bool):
        self._dirty = bool

    def get_arrivals(self, plot=False):
        if not self.metadata:
            print('get_arrivals: No metadata given. Return!')
            return
        func = {True: self.model.get_ray_paths_geo,
                False: self.model.get_travel_times_geo}
        phases = self.prepare_phases()
        trace = self.data.traces[0]
        station_id = trace.get_id()
        starttime = trace.stats.starttime
        station_coords = self.metadata.get_coordinates(station_id, starttime)
        if not station_coords:
            print('get_arrivals: No station coordinates found. Return!')
            return
        origins = self.pylot_event.origins
        if phases == ['None', 'None']:
            print("get_arrivals: Creation info (manual or auto) not available! Return!")
            return
        if origins:
            source_origin = origins[0]
        else:
            print('get_arrivals: No source origin given. Return!')
            return
            # raise ValueError('No source origin given.')
        arrivals = func[plot](source_origin.depth,
                              source_origin.latitude,
                              source_origin.longitude,
                              station_coords.get('latitude'),
                              station_coords.get('longitude'),
                              phases)
        self.arrivals = arrivals

    @staticmethod
    def prepare_phases():
        settings = QtCore.QSettings()
        p_phases = settings.value('p_phases')
        s_phases = settings.value('s_phases')
        phases = ''
        if not p_phases and not s_phases:
            print('No phases for TauPy selected in Preferences.')
        if p_phases and s_phases:
            phases = p_phases + ',' + s_phases
        elif p_phases and not s_phases:
            phases = p_phases
        elif s_phases and not p_phases:
            phases = s_phases
        phases = phases.split(',')
        phases = [phase.strip() for phase in phases]
        return phases

    def drawArrivals(self, textOnly=False):
        if not self.arrivals:
            return
        ax = self.multicompfig.axes[0]
        if not textOnly:
            ylims = self.getGlobalLimits(ax, 'y')
        else:
            ylims = self.multicompfig.getYLims(ax)
        stime = self.getStartTime()
        source_origin = self.pylot_event.origins[0]
        source_time = source_origin.time
        for arrival in self.arrivals:
            arrival_time_abs = source_time + arrival.time
            time_rel = arrival_time_abs - stime
            if not textOnly:
                ax.plot([time_rel, time_rel], ylims, '0.3', linestyle='dashed')
            self.arrivalsText.append(ax.text(time_rel, ylims[0], arrival.name, color='0.5'))

    def drawArrivalsText(self):
        return self.drawArrivals(textOnly=True)

    def refreshArrivalsText(self, event=None):
        self.removeArrivalsText()
        self.drawArrivalsText()

    def removeArrivalsText(self):
        for textItem in self.arrivalsText:
            try:
                textItem.remove()
            except:
                pass
        self.arrivalsText = []

    def addPickPhases(self, menuBar):
        settings = QtCore.QSettings()
        p_phases = settings.value('p_phases')
        s_phases = settings.value('s_phases')

        if p_phases:
            p_phases = p_phases.split(',')
        else:
            p_phases = []
        if s_phases:
            s_phases = s_phases.split(',')
        else:
            s_phases = []

        phases = {'P': p_phases,
                  'S': s_phases}
        if not 'P' in phases['P'] and not 'p' in phases['P']:
            phases['P'] = ['P'] + phases['P']
        if not 'S' in phases['S'] and not 's' in phases['S']:
            phases['S'] = ['S'] + phases['S']

        picksMenu = menuBar.addMenu('Picks')
        self.picksActions = {}

        # dictionary points on corresponding phase_select function
        phaseSelect = {'P': self.p_phase_select,
                       'S': self.s_phase_select}

        hotkeys = {'P': [1, 2, 3, 4, 'q', 'w', 'e', 'r'],
                   'S': [5, 6, 7, 8, 't', 'z', 'u', 'i']}

        # loop over P and S (use explicit list instead of iter over dict.keys to keep order)
        for phaseIndex, phaseID in enumerate(['P', 'S']):
            # loop through phases in list
            for index, phase in enumerate(phases[phaseID]):
                # remove zeros
                phase = phase.strip()
                # add hotkeys
                try:
                    shortcut = str(hotkeys[phaseID][index])
                except IndexError:
                    shortcut = None

                # create action and add to menu
                # phase name transferred using lambda function
                slot = lambda ph=phase, phID=phaseID: phaseSelect[phID](ph)
                picksAction = createAction(parent=self, text=phase,
                                           slot=slot,
                                           shortcut=shortcut)
                picksMenu.addAction(picksAction)
                self.picksActions[str(phase)] = picksAction  # save action in dictionary
            if phaseIndex == 0:
                picksMenu.addSeparator()

        filterOptionsAction = createAction(parent=self, text="&Filter parameter ...",
                                           slot=self.filterOptions,
                                           shortcut='Ctrl+F',
                                           icon=self.orig_parent.filter_icon)
        filterMenu = menuBar.addMenu('Filter')
        filterMenu.addAction(self.filterActionP)
        filterMenu.addAction(self.filterActionS)
        filterMenu.addAction(self.autoFilterAction)
        filterMenu.addAction(filterOptionsAction)

    def filterOptions(self):
        if self.orig_parent.adjustFilterOptions():
            phase = None
            if self.filterActionP.isChecked():
                phase = 'P'
            elif self.filterActionS.isChecked():
                phase = 'S'
            if phase:
                self.plotWFData(phase=phase, filter=True)

    def disable_ar_buttons(self):
        self.enable_ar_buttons(False)

    def enable_ar_buttons(self, bool=True):
        self.accept_button.setEnabled(bool)
        self.reject_button.setEnabled(bool)

    def p_phase_select(self, phase):
        if not self.p_button.isChecked():
            self.p_button.setEnabled(True)
            self.p_button.setChecked(True)
            self.p_button.setText(phase)
        else:
            if str(phase) == str(self.p_button.text()):
                self.reset_p_button()
            else:
                self.p_button.setText(phase)
        self.p_clicked()

    def s_phase_select(self, phase):
        if not self.s_button.isChecked():
            self.s_button.setEnabled(True)
            self.s_button.setChecked(True)
            self.s_button.setText(phase)
        else:
            if str(phase) == str(self.s_button.text()):
                self.reset_s_button()
            else:
                self.s_button.setText(phase)
        self.s_clicked()

    def p_clicked(self):
        if self.p_button.isChecked():
            self.reset_s_button()
            self.s_button.setEnabled(False)
            self.init_p_pick()
        else:
            self.leave_picking_mode()

    def s_clicked(self):
        if self.s_button.isChecked():
            self.reset_p_button()
            self.p_button.setEnabled(False)
            self.init_s_pick()
        else:
            self.leave_picking_mode()

    def init_p_pick(self):
        color = pick_color('manual', 'P')
        self.set_button_border_color(self.p_button, color)
        self.currentPhase = str(self.p_button.text())
        self.activatePicking()

    def init_s_pick(self):
        color = pick_color('manual', 'S')
        self.set_button_border_color(self.s_button, color)
        self.currentPhase = str(self.s_button.text())
        self.activatePicking()

    @staticmethod
    def getPhaseID(phase):
        return identifyPhaseID(phase)

    def set_button_border_color(self, button, color=None):
        '''
        Set background color of a button.
        button: type = QtGui.QAbstractButton
        color: type = QtGui.QColor or type = str (RGBA)
        '''
        if type(color) == QtGui.QColor:
            button.setStyleSheet({'QPushButton{background-color:transparent}'})
            palette = button.palette()
            role = button.backgroundRole()
            palette.setColor(role, color)
            button.setPalette(palette)
            button.setAutoFillBackground(True)
        elif type(color) == str:
            button.setStyleSheet('QPushButton{border-color: %s}' % color)
        elif type(color) == tuple:
            button.setStyleSheet('QPushButton{border-color: rgba%s}' % str(color))
        elif not color:
            button.setStyleSheet(self.orig_parent._style['stylesheet'])

    def reset_p_button(self):
        self.set_button_border_color(self.p_button)
        self.p_button.setEnabled(True)
        self.p_button.setChecked(False)
        self.p_button.setText('P')

    def reset_s_button(self):
        self.set_button_border_color(self.s_button)
        self.s_button.setEnabled(True)
        self.s_button.setChecked(False)
        self.s_button.setText('S')

    def leave_picking_mode(self):
        self.currentPhase = None
        self.reset_p_button()
        self.reset_s_button()
        self.resetZoom()
        self.refreshPlot()
        self.deactivatePicking()

    def activatePicking(self):
        self.leave_rename_phase()
        self.renamePhaseAction.setEnabled(False)
        self.referenceChannel.setEnabled(False)
        self.scaleChannel.setEnabled(False)
        phase = self.currentPhase
        phaseID = self.getPhaseID(phase)
        if not phaseID:
            self.warn_unknown_phase(phase)
            self.leave_picking_mode()
            return
        color = pick_color_plt('manual', phaseID)
        self.multicompfig.set_frame_color(color)
        self.multicompfig.set_frame_linewidth(1.5)
        if self.zoomAction.isChecked():
            self.zoomAction.trigger()
        self.multicompfig.disconnectEvents()
        self.cidpress = self.multicompfig.connectPressEvent(self.setIniPick)
        if self.autoFilterAction.isChecked():
            for filteraction in [self.filterActionP, self.filterActionS]:
                filteraction.setChecked(False)
            self.filterWFData()
            self.draw()
        else:
            self.draw()
        self.disconnect_pick_delete()

    def deactivatePicking(self):
        defaultcolor = self.orig_parent._style['linecolor']['rgba_mpl']
        self.multicompfig.set_frame_color(defaultcolor)
        self.multicompfig.set_frame_linewidth(1)

        self.disconnectPressEvent()
        self.multicompfig.connectEvents()
        self.renamePhaseAction.setEnabled(True)
        self.referenceChannel.setEnabled(True)
        self.scaleChannel.setEnabled(True)
        self.connect_pick_delete()
        self.draw()

    def getParameter(self):
        return self.parameter

    def getStartTime(self):
        return self.stime

    def getEndTime(self):
        return self.etime

    def getComponents(self):
        return self.components

    def getStation(self):
        return '{}.{}.{}'.format(self.network, self.station, self.location)

    def getChannelID(self, key):
        if key < 0: key = 0
        return self.multicompfig.getPlotDict()[int(key)].split('.')[-1]

    def getTraceID(self, channels):
        plotDict = self.multicompfig.getPlotDict()
        traceIDs = []
        for channel in channels:
            channel = channel.upper()
            for seed_id in plotDict.items():
                channelID = seed_id.split('.')[-1]
                if channelID[1].upper().endswith(channel):
                    traceIDs.append(seed_id)
        return traceIDs

    def getUser(self):
        return self._user

    def getFilterOptions(self, phase, gui_filter=False):
        settings = QSettings()
        phaseID = self.getPhaseID(phase)

        if get_bool(settings.value('useGuiFilter')) or gui_filter:
            filteroptions = self.filteroptions[phaseID]
        else:
            filteroptions = getAutoFilteroptions(phaseID, self.parameter)

        if type(filteroptions) == dict:
            freq = [filteroptions['freqmin'], filteroptions['freqmax']]
            order = filteroptions['corners']
            filtertype = filteroptions['type']
            return FilterOptions(filtertype=filtertype, freq=freq, order=order)
        else:
            return filteroptions

    def getXLims(self):
        return self.cur_xlim

    def getYLims(self):
        return self.cur_ylim

    def setXLims(self, limits):
        self.cur_xlim = limits

    def setYLims(self, limits):
        self.cur_ylim = limits

    def getGlobalLimits(self, ax, axis):
        return self.multicompfig.getGlobalLimits(ax, axis)

    def getWFData(self):
        return self.data

    def getWFDataComp(self):
        if self.showCompData:
            return self.data_compare
        else:
            return Stream()

    def selectWFData(self, channel):
        component = channel[-1].upper()
        wfdata = Stream()

        def selectTrace(tr, components):
            if tr.stats.channel[-1].upper() in components:
                return tr

        if component == 'E' or component == 'N':
            for trace in self.getWFData():
                trace = selectTrace(trace, 'NE')
                if trace:
                    wfdata.append(trace)
        elif component == '1' or component == '2':
            for trace in self.getWFData():
                trace = selectTrace(trace, '12')
                if trace:
                    wfdata.append(trace)
        else:
            wfdata = self.getWFData().select(component=component)
        return wfdata

    def getPicks(self, picktype='manual'):
        if picktype == 'manual':
            return self.picks
        elif picktype == 'auto':
            return self.autopicks
        else:
            raise TypeError('Unknown picktype {0}'.format(picktype))

    def resetPicks(self):
        self.picks = {}

    def delPicks(self):
        self.resetPicks()
        self.refreshPlot()

    def initRenamePhase(self):
        if self.renamePhaseAction.isChecked():
            self.multicompfig.disconnectEvents()
            self.multicompfig.set_frame_color('orange')
            self.draw()
            self.statusbar.showMessage('Click on a phase you want to rename.')
        else:
            self.multicompfig.set_frame_color()
            self.multicompfig.connectEvents()
            self.draw()

    def getPickPhases(self, data, phase):
        st = Stream()
        phases = {'P': self.pChannels,
                  'S': self.sChannels}
        if not phase in phases.keys():
            raise ValueError('Unknown phase ID {}'.format(phase))
        for action in phases[phase].actions():
            if action.isChecked():
                st += data.select(channel=action.text())
        return st

    @staticmethod
    def calcNoiseScaleFactor(noiselevel, zoomfactor=5., norm=1):
        # calculate factor to upscale a trace normed to 'norm' in a way that all values
        # zoomfactor*noiselevel are found within -0.5*norm and 0.5*norm
        scaleFactor = (norm / 2.) / (zoomfactor * noiselevel)
        return scaleFactor

    def setIniPick(self, gui_event):
        self.multicompfig.set_frame_color('green')
        trace_number = round(gui_event.ydata)

        self.multicompfig.disconnectEvents()
        self.disconnectPressEvent()
        self.cidpress = self.multicompfig.connectPressEvent(self.setPick)

        if self.getPhaseID(self.currentPhase) == 'P':
            self.set_button_border_color(self.p_button, 'green')
            self.setIniPickP(gui_event)
        elif self.getPhaseID(self.currentPhase) == 'S':
            self.set_button_border_color(self.s_button, 'green')
            self.setIniPickS(gui_event)

        self.zoomAction.setEnabled(False)

        # reset labels
        self.setPlotLabels()
        self.draw()

    def currentFilterPhase(self):
        filterphase = None
        if self.filterActionP.isChecked():
            filterphase = 'P'
        elif self.filterActionS.isChecked():
            filterphase = 'S'
        return filterphase

    def getNoiseWin(self, phase):
        twins_phase = {'P': 'tsnrz',
                       'S': 'tsnrh'}

        return self.parameter.get(twins_phase[phase])[:3]

    def setIniPickP(self, gui_event):
        self.setIniPickPS(gui_event, phase='P')

    def setIniPickS(self, gui_event):
        self.setIniPickPS(gui_event, phase='S')

    def setIniPickPS(self, gui_event, phase):
        phase = self.getPhaseID(phase)

        nfac_phase = {'P': 'nfacP',
                      'S': 'nfacS'}

        parameter = self.parameter
        ini_pick = gui_event.xdata

        nfac = parameter.get(nfac_phase[phase])
        noise_win, gap_win, signal_win = self.getNoiseWin(phase)

        stime = self.getStartTime()

        # copy wfdata for plotting
        wfdata = self.getWFData().copy()
        wfdata_comp = self.getWFDataComp().copy()
        wfdata = self.getPickPhases(wfdata, phase)
        wfdata_comp = self.getPickPhases(wfdata_comp, phase)
        for wfd in [wfdata, wfdata_comp]:
            if wfd:
                wfd.normalize()

        if not wfdata:
            QtWidgets.QMessageBox.warning(self, 'No channel to plot',
                                          'No channel to plot for phase: {}. '
                                          'Make sure to select the correct channels for P and S '
                                          'in the menu in the top panel.'.format(phase))
            self.leave_picking_mode()
            return

        # filter wfdata and trace on which is picked prior to determination of SNR
        filterphase = self.currentFilterPhase()
        if filterphase:
            filteroptions = self.getFilterOptions(filterphase).parseFilterOptions()
            try:
                for wfd in [wfdata, wfdata_comp]:
                    if wfd:
                        wfd.detrend('linear')
                        wfd.filter(**filteroptions)
                # wfdata.filter(**filteroptions)# MP MP removed filtering of original wfdata
            except ValueError as e:
                self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Information,
                                                 'Denied',
                                                 'setIniPick{}: Could not filter waveform: {}'.format(phase, e))
                self.qmb.show()

        snr = []
        noiselevels = {}
        # determine SNR and noiselevel
        for trace in wfdata.traces:
            st = wfdata.select(channel=trace.stats.channel)
            stime_diff = trace.stats.starttime - stime
            result = getSNR(st, (noise_win, gap_win, signal_win), ini_pick - stime_diff)
            snr.append(result[0])
            noiselevel = result[2]
            if noiselevel:
                noiselevel *= nfac
            else:
                noiselevel = nfac
            noiselevels[trace.stats.channel] = noiselevel

        # prepare plotting of wfdata
        for wfd in [wfdata, wfdata_comp]:
            if wfd:
                for trace in wfd:
                    t = prep_time_axis(trace.stats.starttime - stime, trace)
                    inoise = getnoisewin(t, ini_pick, noise_win, gap_win)
                    trace = demeanTrace(trace, inoise)
                    # upscale trace wfdata in a way that each trace is vertically zoomed to noiselevel*factor
                    channel = trace.stats.channel
                    noiselevel = noiselevels[channel]
                    noiseScaleFactor = self.calcNoiseScaleFactor(noiselevel, zoomfactor=5.)
                    trace.data *= noiseScaleFactor
                    noiselevels[channel] *= noiseScaleFactor

        mean_snr = np.mean(snr)
        x_res = getResolutionWindow(mean_snr, parameter.get('extent'))

        xlims = [ini_pick - x_res, ini_pick + x_res]
        ylims = list(np.array([-.5, .5]) + [0, len(wfdata) - 1])

        title = self.getStation() + ' picking mode'
        title += ' | SNR: {}'.format(mean_snr)
        if filterphase:
            filtops_str = transformFilteroptions2String(filteroptions)
            title += ' | Filteroptions: {}'.format(filtops_str)

        plot_additional = bool(self.referenceChannel.currentText())
        additional_channel = self.referenceChannel.currentText()
        self.multicompfig.plotWFData(wfdata=wfdata,
                                     wfdata_compare=wfdata_comp,
                                     title=title,
                                     zoomx=xlims,
                                     zoomy=ylims,
                                     noiselevel=noiselevels,
                                     scaleddata=True,
                                     iniPick=ini_pick,
                                     plot_additional=plot_additional,
                                     additional_channel=additional_channel,
                                     snr=mean_snr)

    def setPick(self, gui_event):

        parameter = self.parameter

        # get axes limits
        self.multicompfig.updateCurrentLimits()

        # setting pick
        pick = gui_event.xdata  # get pick time relative to the traces timeaxis not to the global
        channel = self.getChannelID(round(gui_event.ydata))
        # TODO: channel ID not correct when calcPlotPositions altered positions?

        # get name of phase actually picked
        phase = self.currentPhase

        # get filter parameter for the phase to be picked
        filterphase = self.currentFilterPhase()
        if filterphase:
            filteroptions = self.getFilterOptions(filterphase).parseFilterOptions()
        else:
            filteroptions = None

        # copy and filter data for earliest and latest possible picks
        wfdata = self.getWFData().copy().select(channel=channel)
        if filteroptions:
            try:
                wfdata.detrend('linear')
                wfdata.filter(**filteroptions)
            except ValueError as e:
                self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Information,
                                                 'Denied', 'setPick: Could not filter waveform: {}'.format(e))
                self.qmb.show()

        # get earliest and latest possible pick and symmetric pick error
        # TODO: Hardcoded channel 3 for Z!
        if wfdata[0].stats.channel[2] == 'Z' or wfdata[0].stats.channel[2] == '3':
            nfac = parameter.get('nfacP')
            TSNR = parameter.get('tsnrz')
        else:
            nfac = parameter.get('nfacS')
            TSNR = parameter.get('tsnrh')

        # return absolute time values for phases
        stime = self.getStartTime()
        stime_diff = wfdata[0].stats.starttime - stime

        [epp, lpp, spe] = earllatepicker(wfdata, nfac, (TSNR[0], TSNR[1], TSNR[2]),
                                         pick - stime_diff, verbosity=1)

        mpp = stime + pick

        if epp:
            epp = stime + epp + stime_diff
        if lpp:
            lpp = stime + lpp + stime_diff

        noise_win, gap_win, signal_win = self.getNoiseWin(phase)
        snr, snrDB, noiselevel = getSNR(wfdata, (noise_win, gap_win, signal_win), pick - stime_diff)
        print('SNR of final pick: {}'.format(snr))
        if snr < 1.5:
            QMessageBox.warning(self, 'SNR too low', 'WARNING! SNR of final pick below 1.5! SNR = {}'.format(snr))

        # get first motion and quality classes
        FM = ''
        if self.getPhaseID(phase) == 'P':
            # get first motion quality of P onset is sufficeint
            minFMweight = parameter.get('minfmweight')
            minFMSNR = parameter.get('minFMSNR')
            quality = get_quality_class(spe, parameter.get('timeerrorsP'))
            if quality <= minFMweight and snr >= minFMSNR:
                FM = fmpicker(self.getWFData().select(channel=channel).copy(), wfdata.copy(), parameter.get('fmpickwin'),
                              pick - stime_diff)

        # save pick times for actual phase
        phasepicks = dict(epp=epp, lpp=lpp, mpp=mpp, spe=spe, fm=FM,
                          picker='manual', channel=channel,
                          network=wfdata[0].stats.network,
                          filteroptions=transformFilteroptions2String(filteroptions))

        saved = self.savePick(phase, phasepicks)
        if saved:
            self.setDirty(True)

        self.disconnectPressEvent()
        self.enable_ar_buttons()
        self.zoomAction.setEnabled(True)

        self.leave_picking_mode()

    def savePick(self, phase, phasepicks):
        if not self.getPhaseID(phase):
            self.warn_unknown_phase(phase)
            return

        self.picks[phase] = phasepicks
        return True

    def warn_unknown_phase(self, phase=None):
        QtWidgets.QMessageBox.warning(self, 'Unknown phase ID',
                                      'Could not identify phase ID: {}.'.format(phase))

    def disconnectPressEvent(self):
        self.multicompfig.mpl_disconnect(self.cidpress)
        self.cidpress = None

    def drawAllPicks(self):
        self.removePhaseText()
        self.drawPicks(picktype='manual')
        self.drawPicks(picktype='auto')
        self.draw()

    def drawPicks(self, phase=None, picktype='manual', textOnly=False, picks=None):
        # plotting picks
        ax = self.multicompfig.axes[0]
        if not textOnly:
            ylims = self.multicompfig.getGlobalLimits(ax, 'y')
        else:
            ylims = ax.get_ylim()
        if not picks:
            if self.getPicks(picktype):
                if phase is not None and not phase == 'SPt':
                    if (type(self.getPicks(picktype)[phase]) is dict
                            or type(self.getPicks(picktype)[phase]) is AttribDict):
                        picks = self.getPicks(picktype)[phase]
                elif phase is None:
                    for phase, picks in self.getPicks(picktype).items():
                        if phase == 'SPt':
                            continue
                        self.drawPicks(phase, picktype, textOnly, picks)
                    return
                else:
                    return
            else:
                return

        # get quality classes
        if self.getPhaseID(phase) == 'P':
            quality = get_quality_class(picks['spe'], self.parameter['timeerrorsP'])
            phaseID = 'P'
        elif self.getPhaseID(phase) == 'S':
            quality = get_quality_class(picks['spe'], self.parameter['timeerrorsS'])
            phaseID = 'S'

        # if no mpp is there, return
        if not picks['mpp']:
            return

        mpp = picks['mpp'] - self.getStartTime()
        if picks['epp'] and picks['lpp'] and not textOnly:
            epp = picks['epp'] - self.getStartTime()
            lpp = picks['lpp'] - self.getStartTime()
        spe = picks['spe']

        if not picktype in ['auto', 'manual']:
            raise TypeError('Unknown picktype {0}'.format(picktype))

        if picktype == 'manual':
            baseorder = 5
        elif picktype == 'auto':
            baseorder = 0

        color = pick_color_plt(picktype, phaseID, quality)
        if not textOnly:
            linestyle_mpp, width_mpp = pick_linestyle_plt(picktype, 'mpp')
            vl = ax.axvline(mpp, ylims[0], ylims[1], color=color, linestyle=linestyle_mpp, linewidth=width_mpp,
                            label='{}-{}-Pick (quality: {})'.format(phase, picktype, quality), picker=5,
                            zorder=baseorder + 9)
            phaseLineKey = '{}-{}'.format(phase, picktype)
            self.phaseLines[phaseLineKey] = vl
            if spe:
                ax.fill_between([mpp - spe, mpp + spe], ylims[0], ylims[1],
                                alpha=.25, color=color, label='{}-{}-SPE'.format(phase, picktype), zorder=baseorder + 1)
            if picks['epp']:
                linestyle_epp, width_epp = pick_linestyle_plt(picktype, 'epp')
                ax.axvline(epp, ylims[0], ylims[1], color=color, linestyle=linestyle_epp,
                           linewidth=width_epp, label='{}-{}-EPP'.format(phase, picktype), zorder=baseorder + 2)
            if picks['lpp']:
                linestyle_lpp, width_lpp = pick_linestyle_plt(picktype, 'lpp')
                ax.axvline(lpp, ylims[0], ylims[1], color=color, linestyle=linestyle_lpp,
                           linewidth=width_lpp, label='{}-{}-LPP'.format(phase, picktype), zorder=baseorder + 2)
            if picktype == 'auto':
                ax.plot(mpp, ylims[1], color=color, marker='v', zorder=baseorder + 3)
                ax.plot(mpp, ylims[0], color=color, marker='^', zorder=baseorder + 3)
        # append phase text (if textOnly: draw with current ylims)
        self.phaseText.append(ax.text(mpp, ylims[1], phase, color=color, zorder=baseorder + 10))
        # indicate first motion
        fm = picks.get('fm')
        if fm:
            self.phaseText.append(ax.text(mpp - 0.03 * mpp, ylims[1] - ylims[1] / 12, fm, color=color,
                                          zorder=baseorder + 10))
        ax.legend(loc=1)

    def connect_mouse_motion(self):
        self.cidmotion = self.multicompfig.mpl_connect(
            'motion_notify_event', self.on_motion)

    def connect_pick_delete(self):
        self.cidpick = self.multicompfig.mpl_connect('pick_event', self.onpick)
        self.cidpick = self.multicompfig.mpl_connect('motion_notify_event', self.on_hover_info)

    def disconnect_pick_delete(self):
        if hasattr(self, 'cidpick'):
            self.multicompfig.mpl_disconnect(self.cidpick)

    def on_motion(self, event):
        x = event.xdata
        if x is not None:
            time_code = 'T = {}, t = {} [s]'.format(self.stime + x, x)
            user_help = ' - Left-Click to Drag | Right-Click to Pan-Zoom |' \
                        ' Mousewheel to Zoom | Middle-Click to Delete Pick'
            self.statusbar.showMessage(time_code + user_help)

    def onpick(self, event):
        if event.mouseevent.button == 1:
            self.onpick_info(event)
        elif event.mouseevent.button == 2:
            self.onpick_delete(event)

    def on_hover_info(self, event):
        if not any([phase.contains(event)[0] for phase in self.phaseLines.values()]):
            return
        x = event.xdata
        if not x:
            return
        allpicks, pick_rel, phase, picktype = self.identify_selected_picks(x)
        if pick_rel is None:
            return
        pick = allpicks[picktype][phase]
        message = '{} {}-pick'.format(picktype, phase)
        if 'mpp' in pick:
            message += ', MPP: {}'.format(pick['mpp'])
        if 'spe' in pick:
            message += ', SPE: {} [s]'.format(pick['spe'])
        if 'filteroptions' in pick:
            message += ', FILTER: {}'.format(pick['filteroptions'])
        x = event.x
        y = event.y
        y = self.size().height() - y
        pt = self.mapToGlobal(QtCore.QPoint(x, y))
        QtWidgets.QToolTip.showText(pt, message)

    def onpick_info(self, event):
        if not event.mouseevent.button == 1:
            return
        x = event.mouseevent.xdata
        allpicks, pick_rel, phase, picktype = self.identify_selected_picks(x)
        if pick_rel is None:
            return
        pick = allpicks[picktype][phase]
        message = '{} {}-pick'.format(picktype, phase)
        if 'mpp' in pick:
            message += ', MPP: {}'.format(pick['mpp'])
        if 'spe' in pick:
            message += ', SPE: {}'.format(pick['spe'])
        if 'filteroptions' in pick:
            message += ', FILTER: {}'.format(pick['filteroptions'])

        if self.renamePhaseAction.isChecked():
            self.renamePhase(picktype, phase)

        self.statusbar.showMessage(message, 10e3)

    def onpick_delete(self, event):
        if not event.mouseevent.button == 2:
            return
        x = event.mouseevent.xdata
        self.remove_pick_by_x(x)
        self.refreshPlot()

    def renamePhase(self, picktype, phase):
        allpicks = {'manual': self.picks,
                    'auto': self.autopicks}
        picks = allpicks[picktype]
        dialog = QtWidgets.QInputDialog(parent=self)
        new_phase, executed = dialog.getText(self, 'Rename phase', 'Rename phase {} to:'.format(phase))
        if executed:
            try:
                self.renamePhaseInDict(picks, phase, new_phase)
            except KeyError as e:
                QtWidgets.QMessageBox.warning(self, 'Could not rename phase',
                                              'Could not rename phase {} to {}: {}'.format(phase, new_phase, e))
        self.leave_rename_phase()
        self.refreshPlot()

    def renamePhaseInDict(self, picks, phase_old, phase_new):
        if phase_new in picks:
            raise KeyError('New phase ID already assigned.')
        picks_new = picks[phase_old].copy()
        saved = self.savePick(phase_new, picks_new)
        if saved:
            self.phaseLines[phase_new] = self.phaseLines.pop(phase_old)
            picks.pop(phase_old)
            self.setDirty(True)

    def leave_rename_phase(self):
        self.renamePhaseAction.setChecked(False)
        self.multicompfig.set_frame_color()
        self.multicompfig.connectEvents()

    def remove_pick_by_x(self, x):
        if not self.picks and not self.autopicks:
            return
        allpicks, pick_rel, phase, picktype = self.identify_selected_picks(x)
        if pick_rel is None:
            return
        # delete the value from corresponding dictionary
        removed_pick = allpicks[picktype].pop(phase)
        # delete line from vlines dictionary
        if phase in self.phaseLines.keys():
            del (self.phaseLines[phase])
        # information output
        msg = 'Deleted {} pick for phase {}, station {} at timestamp {} (relative time: {} s)'
        print(msg.format(picktype, phase, '{}.{}'.format(self.network, self.station),
                         self.getStartTime() + pick_rel, pick_rel))
        removed_pick['station'] = self.station
        self.removed_picks.append(removed_pick)
        self.setDirty(True)

    def identify_selected_picks(self, x):
        # init empty list and get stat5ion starttime
        X = []
        pick_rel = phase = picktype = None
        starttime = self.getStartTime()
        # init dictionaries to iterate through and iterate over them
        allpicks = {'manual': self.picks,
                    'auto': self.autopicks}
        for picktype in allpicks.keys():
            picks = allpicks[picktype]
            for phase in picks:
                if not type(picks[phase]) in [dict, AttribDict]:
                    continue
                if not picks[phase]['mpp']:
                    continue
                pick_rel = picks[phase]['mpp'] - starttime
                # add relative pick time, phaseID and picktype index
                X.append((pick_rel, phase, picktype))
        if len(X) > 0:
            # find index and value closest to x
            index, value = min(enumerate([val[0] for val in X]), key=lambda y: abs(y[1] - x))
            # unpack the found value
            pick_rel, phase, picktype = X[index]
        return allpicks, pick_rel, phase, picktype

    def drawPhaseText(self):
        self.drawPicks(picktype='manual', textOnly=True)
        self.drawPicks(picktype='auto', textOnly=True)

    def removePhaseText(self):
        for textItem in self.phaseText:
            try:
                textItem.remove()
            except:
                pass
        self.phaseText = []

    def refreshPhaseText(self, event=None):
        self.removePhaseText()
        self.drawPhaseText()

    def togglePickBlocker(self):
        return not self.pick_block

    def filterWFData(self, phase=None):
        if not phase:
            phase = self.currentPhase
            if self.getPhaseID(phase) == 'P':
                self.filterActionP.setChecked(True)
                self.filterP()
            elif self.getPhaseID(phase) == 'S':
                self.filterActionS.setChecked(True)
                self.filterS()
            return
        self.plotWFData(phase=phase, filter=True)

    def plotWFData(self, phase=None, filter=False):
        if self.pick_block:
            return
        self.cur_xlim = self.multicompfig.axes[0].get_xlim()
        self.cur_ylim = self.multicompfig.axes[0].get_ylim()
        # self.multicompfig.updateCurrentLimits()
        wfdata = self.getWFData().copy()
        wfdata_comp = self.getWFDataComp().copy()
        title = self.getStation()
        if filter:
            filtoptions = None
            if phase:
                filtoptions = self.getFilterOptions(self.getPhaseID(phase), gui_filter=True).parseFilterOptions()

            if filtoptions is not None:
                for wfd in [wfdata, wfdata_comp]:
                    if wfd:
                        wfd.detrend('linear')
                        wfd.taper(0.02, type='cosine')
                        wfd.filter(**filtoptions)
                filtops_str = transformFilteroptions2String(filtoptions)
                title += ' | Filteroptions: {}'.format(filtops_str)

        if self.wftype is not None:
            title += ' | ({})'.format(self.wftype)

        plot_additional = bool(self.referenceChannel.currentText())
        additional_channel = self.referenceChannel.currentText()
        scale_channel = self.scaleChannel.currentText()
        self.multicompfig.plotWFData(wfdata=wfdata, wfdata_compare=wfdata_comp,
                                     title=title,
                                     zoomx=self.getXLims(),
                                     zoomy=self.getYLims(),
                                     plot_additional=plot_additional,
                                     additional_channel=additional_channel,
                                     scaleToChannel=scale_channel)
        self.setPlotLabels()
        self.drawAllPicks()
        self.drawArrivals()
        self.draw()

    def filterP(self):
        self.filterActionS.setChecked(False)
        if self.filterActionP.isChecked():
            self.filterWFData('P')
        else:
            self.refreshPlot()

    def filterS(self):
        self.filterActionP.setChecked(False)
        if self.filterActionS.isChecked():
            self.filterWFData('S')
        else:
            self.refreshPlot()

    def toggleAutoFilter(self):
        settings = QSettings()
        settings.setValue('autoFilter', self.autoFilterAction.isChecked())

    @staticmethod
    def updateChannelSettingsP(action):
        settings = QSettings()
        settings.setValue('p_channel_{}'.format(action.text()), action.isChecked())

    @staticmethod
    def updateChannelSettingsS(action):
        settings = QSettings()
        settings.setValue('s_channel_{}'.format(action.text()), action.isChecked())

    @staticmethod
    def getChannelSettingsP(channel):
        settings = QSettings()
        rval = get_bool(settings.value('p_channel_{}'.format(channel)))
        compclass = SetChannelComponents.from_qsettings(settings)
        components = ['Z']
        for component in components[:]:
            components.append(compclass.getCompPosition(component))
        if not rval in [True, False]:
            if any([channel.endswith(component) for component in components]):
                rval = True
            else:
                rval = False
        return rval

    @staticmethod
    def getChannelSettingsS(channel):
        settings = QSettings()
        rval = get_bool(settings.value('s_channel_{}'.format(channel)))
        compclass = SetChannelComponents.from_qsettings(settings)
        components = ['N', 'E']
        for component in components[:]:
            components.append(compclass.getCompPosition(component))
        if not rval in [True, False]:
            if any([channel.endswith(component) for component in components]):
                rval = True
            else:
                rval = False
        return rval

    def resetPlot(self):
        self.resetZoom()
        self.refreshPlot()

    def switchCompData(self):
        self.showCompData = self.compareCB.isChecked()

    def refreshPlot(self):
        if self.autoFilterAction.isChecked():
            self.filterActionP.setChecked(False)
            self.filterActionS.setChecked(False)
        # data = self.getWFData().copy()
        # title = self.getStation()
        filter = False
        phase = None
        if self.filterActionP.isChecked():
            phase = 'P'
            filter = True
        if self.filterActionS.isChecked():
            phase = 'S'
            filter = True
        self.plotWFData(phase=phase, filter=filter)

    def resetZoom(self):
        ax = self.multicompfig.axes[0]
        self.multicompfig.setXLims(ax, self.multicompfig.getGlobalLimits(ax, 'x'))
        self.multicompfig.setYLims(ax, self.multicompfig.getGlobalLimits(ax, 'y'))
        if not self.zoomAction.isChecked():
            self.multicompfig.connectEvents()

    def setPlotLabels(self):
        # get channel labels
        pos = self.multicompfig.getPlotDict().keys()
        labels = [self.multicompfig.getPlotDict()[key].split('.')[-1] for key in pos]

        ax = self.multicompfig.figure.axes[0]

        # set channel labels
        self.multicompfig.setYTickLabels(pos, labels)

    def zoom(self):
        if self.zoomAction.isChecked() and self.pick_block:
            self.zoomAction.setChecked(False)
        elif self.zoomAction.isChecked():
            self.multicompfig.disconnectEvents()
            self.figToolBar.zoom()
        else:
            self.figToolBar.zoom()
            self.multicompfig.connectEvents()

    def draw(self):
        self.multicompfig.draw()

    def apply(self):
        picks = self.getPicks()
        self.update_picks.emit(picks)
        # for pick in picks:
        #     print(pick, picks[pick])

    def discard(self):
        self.picks = self._init_picks
        self.autopicks = self._init_autopicks
        self.update_picks.emit(self.picks)
        # for pick in picks:
        #     print(pick, picks[pick])

    def reject(self):
        self.discard()
        if not self._embedded:
            QDialog.reject(self)
        else:
            self.refreshPlot()
            self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Information,
                                             'Denied', 'New picks rejected!')
            self.qmb.show()

    def accept(self):
        self.apply()
        if not self._embedded:
            QDialog.accept(self)
        else:
            self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Information,
                                             'Accepted', 'New picks applied!')
            self.qmb.show()


class PhasePlotWidget(FigureCanvas):
    def __init__(self, parent=None):
        self._parent = parent
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111, projection='polar')
        self.new = True
        super(PhasePlotWidget, self).__init__(self.fig)


class CanvasWidget(QWidget):
    '''
    '''

    def __init__(self, parent, canvas):
        QtWidgets.QWidget.__init__(self, parent)  # , 1)
        canvas = canvas
        self.main_layout = QtWidgets.QVBoxLayout()
        self.setLayout(self.main_layout)
        self.main_layout.addWidget(canvas)
        canvas.setZoomBorders2content()


class MultiEventWidget(QWidget):
    start = Signal()
    '''
    '''

    def __init__(self, options=None, parent=None, windowflag=Qt.Window):
        QtWidgets.QWidget.__init__(self, parent, windowflag)

        self.options = options
        self.setupUi()
        # set initial size
        self.resize(1280, 720)
        self.center()

    def center(self):
        fm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        fm.moveCenter(centerPoint)
        self.move(fm.topLeft())

    def setupUi(self):
        # init main layout
        self.main_layout = QtWidgets.QVBoxLayout()
        self.setLayout(self.main_layout)
        # init main splitter
        self.main_splitter = QtWidgets.QSplitter()
        self.main_splitter.setChildrenCollapsible(False)

        self.init_checkboxes()

        self.eventbox = QtWidgets.QComboBox()
        self.button_clear = QtWidgets.QPushButton('Clear')

        self.main_layout.insertWidget(1, self.main_splitter)

        self.main_layout.setStretch(0, 0)
        self.main_layout.setStretch(1, 1)
        self.main_splitter.setStretchFactor(0, 1)
        self.main_splitter.setStretchFactor(1, 2)

    def init_checkboxes(self):
        self.rb_layout = QtWidgets.QHBoxLayout()

        self.rb_dict = {}

        self.start_button = QtWidgets.QPushButton('Start')

        for index, (key, func, color) in enumerate(self.options):
            rb = QtWidgets.QRadioButton(key)
            rb.toggled.connect(self.check_rb_selection)
            if color:
                color = 'rgba{}'.format(color)
            else:
                color = 'transparent'
            rb.setStyleSheet('QRadioButton{'
                             'background-color: %s;'
                             'border-style:outset;'
                             'border-width:1px;'
                             'border-radius:5px;'
                             'padding:5px;'
                             '}' % str(color))
            if index == 0:
                rb.setChecked(True)
            self.rb_dict[key] = rb
            self.rb_layout.insertWidget(index, rb)
            self.rb_layout.setStretch(index, 0)

        self.pb = QtWidgets.QProgressBar()
        self.pb.setRange(0, 0)
        self.pb.setVisible(False)

        # space holder for progressbar
        self._pb_space = QtWidgets.QWidget()

        self.rb_layout.addWidget(self.start_button)

        self.rb_layout.addWidget(self.pb)
        self.rb_layout.addWidget(self._pb_space)

        self.rb_layout.setStretch(len(self.options) + 1, 1)
        self.rb_layout.setStretch(len(self.options) + 2, 1)

        self.main_layout.insertLayout(0, self.rb_layout)

    def refresh_tooltips(self):
        for key, func, color in self.options:
            eventlist = func()
            if not type(eventlist) == list:
                eventlist = [eventlist]
            tooltip = ''
            for index, event in enumerate(eventlist):
                if not event:
                    continue
                tooltip += '{}'.format(event.pylot_id)
                if not index + 1 == len(eventlist):
                    tooltip += '\n'
            if not tooltip:
                tooltip = 'No events for this selection'
            self.rb_dict[key].setToolTip(tooltip)
        self.check_rb_selection()

    def check_rb_selection(self):
        for rb in self.rb_dict.values():
            if rb.isChecked():
                check_events = (rb.toolTip() == 'No events for this selection')
                self.start_button.setEnabled(not check_events)

    def enable(self, bool):
        for rb in self.rb_dict.values():
            rb.setEnabled(bool)
        self.start_button.setEnabled(bool)
        self.pb.setVisible(not bool)
        self._pb_space.setVisible(bool)
        self.eventbox.setEnabled(bool)
        self.button_clear.setEnabled(bool)


class AutoPickWidget(MultiEventWidget):
    '''
    '''

    def __init__(self, parent, options):
        MultiEventWidget.__init__(self, options, parent, Qt.Window)
        self.events2plot = {}
        self.connect_buttons()
        self.init_plot_layout()
        self.init_log_layout()
        self.reinitEvents2plot()
        self.setWindowTitle('Autopick events interactively')
        self.set_main_stretch()

    def connect_buttons(self):
        self.start_button.clicked.connect(self.run)
        self.button_clear.clicked.connect(self.reinitEvents2plot)

    def init_plot_layout(self):
        # init tab widget
        self.tab_plots = QtWidgets.QTabWidget()
        self.gb_plots = QtWidgets.QGroupBox('Plots')
        self.gb_plots.setMinimumSize(100, 100)
        self.main_splitter.insertWidget(1, self.gb_plots)
        self.plot_layout = QtWidgets.QVBoxLayout()
        self.plot_layout.insertWidget(1, self.tab_plots)
        self.gb_plots.setLayout(self.plot_layout)

    def init_log_layout(self):
        self.gb_log = QtWidgets.QGroupBox('Log')
        self.gb_log.setMinimumSize(100, 100)
        self.main_splitter.insertWidget(0, self.gb_log)

    def insert_log_widget(self, widget):
        vl = QtWidgets.QVBoxLayout()
        vl.addWidget(widget)
        self.gb_log.setLayout(vl)

    def add_plot_widget(self, widget, key, eventID):
        eventID += ' [picked: {}]'.format(time.strftime('%X %x %z'))
        if not eventID in self.events2plot.keys():
            self.events2plot[eventID] = {}
        self.events2plot[eventID][key] = widget

    def generate_combobox(self):
        self.eventbox.clear()
        for eventID, widgets in self.events2plot.items():
            self.eventbox.addItem(str(eventID), widgets)
        self.eventbox.currentIndexChanged.connect(self.draw_plots)
        self.draw_plots()

    def draw_plots(self, index=0):
        self.refresh_plot_tabs()
        widgets = self.eventbox.itemData(index)
        if not widgets:
            return
        for key, widget in widgets.items():
            self.tab_plots.addTab(widget, str(key))

    def update_plots(self):
        self.refresh_plot_tabs()
        if len(self.events2plot) > 0:
            self.eventbox_layout = QtWidgets.QHBoxLayout()
            self.generate_combobox()
            self.eventbox_layout.addWidget(self.eventbox)
            self.eventbox_layout.addWidget(self.button_clear)
            self.eventbox_layout.setStretch(0, 1)
            self.plot_layout.insertLayout(0, self.eventbox_layout)

    def set_main_stretch(self):
        self.main_layout.setStretch(0, 0)
        self.main_layout.setStretch(1, 1)

    def reinitEvents2plot(self):
        for eventID, eventDict in self.events2plot.items():
            for widget_key, widget in eventDict.items():
                widget.setParent(None)
        self.events2plot = {}
        self.eventbox.clear()
        self.refresh_plot_tabs()

    def refresh_plot_tabs(self):
        self.tab_plots.clear()

    def run(self):
        self.refresh_plot_tabs()
        self.start.emit()


class CompareEventsWidget(MultiEventWidget):
    '''
    '''

    def __init__(self, parent, options, eventdict, comparisons):
        MultiEventWidget.__init__(self, options, parent, Qt.Window)
        self.eventdict = eventdict
        self.comparisons = comparisons
        self.compare_widget = QtWidgets.QWidget()
        self.init_eventbox()
        self.init_event_area()
        self.fill_eventbox()
        self.connect_buttons()
        self.setWindowTitle('Compare events')
        self.set_main_stretch()
        self.center()

    def center(self):
        fm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        fm.moveCenter(centerPoint)
        self.move(fm.topLeft())

    def connect_buttons(self):
        self.start_button.clicked.connect(self.run)
        self.start_button.setText('Show Histograms')

    def init_event_area(self):
        self.event_layout = QVBoxLayout()
        self.event_layout.insertWidget(0, self.eventbox)
        self.event_area = QGroupBox('Single Event')
        self.event_area.setLayout(self.event_layout)
        self.main_layout.insertWidget(1, self.event_area)

    def init_eventbox(self):
        self.eventbox_layout = QtWidgets.QHBoxLayout()
        self.eventbox_layout.addWidget(self.eventbox)
        self.eventbox.currentIndexChanged.connect(self.update_comparison)

    def fill_eventbox(self):
        event_ids = list(self.eventdict.keys())
        for event_id in sorted(event_ids):
            self.eventbox.addItem(str(event_id))
        self.update_comparison()

    def update_eventbox(self):
        self.eventbox.clear()
        self.fill_eventbox()

    def update_comparison(self):
        self.compare_widget.setParent(None)
        self.compare_widget = ComparisonWidget(
            self.comparisons[self.eventbox.currentText()], self, Qt.Widget)
        self.event_layout.insertWidget(1, self.compare_widget)
        self.set_main_stretch()

    def set_main_stretch(self):
        self.main_layout.setStretch(0, 0)
        self.main_layout.setStretch(1, 1)
        self.main_layout.setStretch(2, 0)
        self.event_layout.setStretch(0, 0)
        self.event_layout.setStretch(1, 1)

    def run(self):
        self.start.emit()


class TuneAutopicker(QWidget):
    update = QtCore.Signal(str)
    '''
    QWidget used to modifiy and test picking parameters for autopicking algorithm.
    :param: parent
    :type: PyLoT Mainwindow
    '''

    def __init__(self, parent, obspy_dmt=False):
        QtWidgets.QWidget.__init__(self, parent, Qt.Window)
        self._style = parent._style
        self.setWindowTitle('PyLoT - Tune Autopicker')
        self.parameter = self.parent()._inputs
        self.fig_dict = self.parent().fig_dict
        self.data = Data()
        self.obspy_dmt = obspy_dmt
        self.wftype = None
        self.pdlg_widget = None
        self.pylot_picks = None
        self.init_main_layouts()
        self.init_eventlist()
        self.init_figure_tabs()
        self.init_stationlist()
        self.init_pbwidget()
        self.connect_signals()
        self.add_parameters()
        self.add_buttons()
        self.add_log()
        self.set_stretch()
        self.resize(1280, 720)
        self.center()
        self._manual_pick_plots = []
        self.fnames = None
        if hasattr(self.parent(), 'metadata'):
            self.metadata = self.parent().metadata
        else:
            self.metadata = None
            # self.setWindowModality(QtCore.Qt.WindowModality.ApplicationModal)
            # self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)

    def center(self):
        fm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        fm.moveCenter(centerPoint)
        self.move(fm.topLeft())

    def set_fig_dict(self, fig_dict):
        for key, value in fig_dict.items():
            if key is not 'mainFig':
                value._tight = True
        self.fig_dict = fig_dict

    def init_main_layouts(self):
        self.main_layout = QtWidgets.QVBoxLayout()
        self.tune_layout = QtWidgets.QHBoxLayout()
        self.trace_layout = QtWidgets.QHBoxLayout()
        self.parameter_layout = QtWidgets.QVBoxLayout()

        self.main_layout.addLayout(self.trace_layout)
        self.main_layout.addLayout(self.tune_layout)
        self.setLayout(self.main_layout)

    def init_eventlist(self):
        self.eventBox = self.parent().createEventBox()
        self.eventBox.setMaxVisibleItems(20)
        self.fill_eventbox()
        self.trace_layout.addWidget(self.eventBox)

    def init_stationlist(self):
        self.stationBox = QtWidgets.QComboBox()
        self.stationBox.setMaxVisibleItems(42)
        self.trace_layout.addWidget(self.stationBox)
        self.fill_stationbox()
        self.figure_tabs.setCurrentIndex(0)

    def connect_signals(self):
        self.eventBox.activated.connect(self.fill_stationbox)
        self.eventBox.activated.connect(self.update_eventID)
        self.eventBox.activated.connect(self.fill_tabs)
        self.stationBox.activated.connect(self.fill_tabs)

    def catch_station_ids(self):
        """
        Fill self.station_ids dictionary.
        The dict keys are the strings of the station ids in format network.station.location, the values are lists of
        filenames, which contain traces of that station.
        E.g. if file A.mseed and B.mseed contain traces from station X.Y.Z, then the dict would contain
        {"X.Y.Z" : ["/full/path/to/A.mseed", "/full/path/to/B.mseed"]}
        """
        self.station_ids = {}
        eventpath = self.get_current_event_fp()
        self.wftype = 'processed' if self.obspy_dmt else ''
        wf_path = os.path.join(eventpath, self.wftype)
        if not os.path.exists(wf_path) and self.obspy_dmt:
            self.wftype = 'raw'
            wf_path = os.path.join(eventpath, self.wftype)
        for filename in os.listdir(wf_path):
            filename = os.path.join(eventpath, self.wftype, filename)
            try:
                st = read(filename, headonly=True)
            except Exception as e:
                print('Warning: Could not read file {} as a stream object: {}'.format(filename, e))
                continue
            for trace in st:
                station_id = trace.get_id()
                station_id = station_id_remove_channel(station_id)
                if not station_id in self.station_ids:
                    self.station_ids[station_id] = []
                self.station_ids[station_id].append(filename)

    def fill_stationbox(self):
        self.stationBox.clear()
        model = self.stationBox.model()

        self.catch_station_ids()
        st_ids_list = list(self.station_ids.keys())
        st_ids_list.sort()
        for station_id in st_ids_list:
            item = QtGui.QStandardItem(station_id)
            if station_id.split('.')[1] in self.get_current_event().pylot_picks:
                item.setBackground(self.parent()._ref_test_colors['ref'])
            model.appendRow(item)

    def load_wf_data(self):
        fnames = self.station_ids[self.get_current_station_id()]
        if not fnames == self.fnames:
            self.fnames = fnames
            self.data.setWFData(fnames)
            wfdat = self.data.getWFData()  # all available streams
            # remove possible underscores in station names
            # wfdat = remove_underscores(wfdat)
            # rotate misaligned stations to ZNE
            # check for gaps and doubled channels
            wfdat, _ = check_for_gaps_and_merge(wfdat)
            # check for nans
            check_for_nan(wfdat)
            check4doubled(wfdat)
            wfdat = check4rotated(wfdat, self.parent().metadata, verbosity=0)
            # trim station components to same start value
            trim_station_components(wfdat, trim_start=True, trim_end=False)

    def init_figure_tabs(self):
        self.figure_tabs = QtWidgets.QTabWidget()
        self.fill_figure_tabs()

    def init_pbwidget(self):
        self.pb_widget = ProgressBarWidget()

    def init_tab_names(self):
        self.ptb_names = ['aicFig', 'slength', 'checkZ4s', 'refPpick', 'el_Ppick', 'fm_picker']
        self.stb_names = ['aicARHfig', 'refSpick', 'el_S1pick', 'el_S2pick']

    def add_parameters(self):
        self.paraBox = PylotParameterWidget(self.parameter, parent=self, windowflag=Qt.Widget)
        self.paraBox.set_tune_mode(True)
        self.update_eventID()
        self.parameter_layout.addWidget(self.paraBox)
        self.parameter_layout.addWidget(self.pb_widget)
        self.tune_layout.insertLayout(1, self.parameter_layout)

    def add_buttons(self):
        self.pick_button = QtWidgets.QPushButton('Pick Trace')
        self.pick_button.setStyleSheet('QPushButton{border-color: rgba(110, 200, 0, 255)}')
        self.pick_button.clicked.connect(self.call_picker)
        self.close_button = QtWidgets.QPushButton('Close')
        self.close_button.clicked.connect(self.hide)
        self.trace_layout.addWidget(self.pick_button)
        self.trace_layout.setStretch(0, 1)
        self.parameter_layout.addWidget(self.close_button)

    def add_log(self):
        self.listWidget = QtWidgets.QListWidget()
        self.figure_tabs.insertTab(4, self.listWidget, 'log')

    def add_log_item(self, text):
        self.listWidget.addItem(text)
        self.listWidget.scrollToBottom()

    def get_current_event(self):
        path = self.get_current_event_fp()
        return self.parent().project.getEventFromPath(path)

    def get_current_event_name(self):
        return self.eventBox.currentText().split('/')[-1].rstrip('*')

    def get_current_event_fp(self):
        return self.eventBox.currentText().rstrip('*')

    def get_current_event_picks(self, station):
        event = self.get_current_event()
        if station in event.pylot_picks.keys():
            return event.pylot_picks[station]

    def get_current_event_autopicks(self, station):
        event = self.get_current_event()
        if event.pylot_autopicks:
            if station in event.pylot_autopicks:
                return event.pylot_autopicks[station]

    def get_current_station(self):
        return str(self.stationBox.currentText()).split('.')[1]

    def get_current_station_id(self):
        return str(self.stationBox.currentText())

    @staticmethod
    def gen_tab_widget(name, canvas):
        widget = QtWidgets.QWidget()
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addWidget(canvas)
        widget.setLayout(v_layout)
        return widget

    def gen_pick_dlg(self):
        if not self.get_current_event():
            if self.pdlg_widget:
                self.pdlg_widget.setParent(None)
            self.pdlg_widget = None
            return
        self.load_wf_data()
        station_id_list = self.get_current_station_id().split(".")
        # this sometimes only contains network.station.
        if len(station_id_list) == 3:
            # location is empty string
            network, station, location = station_id_list
        elif len(station_id_list) == 4:
            # location contains text
            network, station, location, _ = station_id_list
        else:
            station = self.get_current_station()
            network = None
            location = None

        wfdata = self.data.getWFData()
        wfdata_comp = self.data.getWFDataComp()
        metadata = self.parent().metadata
        event = self.get_current_event()
        filteroptions = self.parent().filteroptions
        wftype = self.wftype if self.obspy_dmt else ''
        self.pickDlg = PickDlg(self.parent(), data=wfdata.select(station=station).copy(),
                               data_comp=wfdata_comp.select(station=station).copy(),
                               station=station, network=network,
                               location=location, parameter=self.parameter,
                               picks=self.get_current_event_picks(station),
                               autopicks=self.get_current_event_autopicks(station),
                               metadata=metadata, event=event, filteroptions=filteroptions,
                               model=self.parameter.get('taup_model'),
                               embedded=True, wftype=wftype)
        self.pickDlg.update_picks.connect(self.picks_from_pickdlg)
        self.pickDlg.update_picks.connect(self.fill_eventbox)
        self.pickDlg.update_picks.connect(self.fill_stationbox)
        self.pickDlg.update_picks.connect(lambda: self.parent().setDirty(True))
        self.pickDlg.update_picks.connect(self.parent().enableSaveEventAction)
        self.pickDlg.update_picks.connect(self.plot_manual_picks_to_figs)
        self.pdlg_widget = QtWidgets.QWidget(self)
        hl = QtWidgets.QHBoxLayout()
        self.pdlg_widget.setLayout(hl)
        hl.addWidget(self.pickDlg)

    def picks_from_pickdlg(self, picks=None):
        seed_id = self.get_current_station_id()
        station = self.get_current_station()
        replot = self.parent().addPicks(station, picks)
        self.get_current_event().setPick(station, picks)
        if self.get_current_event() == self.parent().get_current_event():
            if replot:
                self.parent().plotWaveformDataThread()
                self.parent().drawPicks()
            else:
                self.parent().drawPicks(seed_id)
            self.parent().draw()

    def clear_plotitem(self, plotitem):
        if type(plotitem) == list:
            for item in plotitem:
                self.clear_plotitem(item)
            return
        try:
            plotitem.remove()
        except Exception as e:
            print('Warning could not remove item {}: {}'.format(plotitem, e))

    def plot_manual_picks_to_figs(self):
        picks = self.get_current_event_picks(self.get_current_station())
        if not picks or not 'P' in picks or not 'S' in picks:
            return
        for plotitem in self._manual_pick_plots:
            self.clear_plotitem(plotitem)
        self._manual_pick_plots = []
        st = self.data.getWFData()
        tr = st.select(station=self.get_current_station())[0]
        starttime = tr.stats.starttime
        # create two lists with figure names and subindices (for subplots) to get the correct axes
        p_axes = [
            ('mainFig', 0),
            ('aicFig', 0),
            ('slength', 0),
            ('refPpick', 0),
            ('el_Ppick', 0),
            ('fm_picker', 0),
            ('fm_picker', 1)]
        s_axes = [
            ('mainFig', 1),
            ('mainFig', 2),
            ('aicARHfig', 0),
            ('refSpick', 0),
            ('el_S1pick', 0),
            ('el_S2pick', 0)]
        qualityPpick = get_quality_class(picks['P']['spe'], self.parameter['timeerrorsP'])
        qualitySpick = get_quality_class(picks['S']['spe'], self.parameter['timeerrorsS'])
        for p_ax in p_axes:
            axes = self.parent().fig_dict[p_ax[0]].axes
            if not axes:
                continue
            ax = axes[p_ax[1]]
            self.plot_manual_pick_to_ax(ax=ax, picks=picks, phase='P',
                                        starttime=starttime, quality=qualityPpick)
        for s_ax in s_axes:
            axes = self.parent().fig_dict[s_ax[0]].axes
            if not axes:
                continue
            ax = axes[s_ax[1]]
            self.plot_manual_pick_to_ax(ax=ax, picks=picks, phase='S',
                                        starttime=starttime, quality=qualitySpick)
        for canvas in self.parent().canvas_dict.values():
            canvas.draw_idle()

    def plot_manual_pick_to_ax(self, ax, picks, phase, starttime, quality):
        mpp = picks[phase]['mpp'] - starttime
        color = pick_color_plt('manual', phase, quality)

        y_top = 0.9 * ax.get_ylim()[1]
        y_bot = 0.9 * ax.get_ylim()[0]
        self._manual_pick_plots.append(ax.vlines(mpp, y_bot, y_top,
                                                 color=color, linewidth=2,
                                                 label='manual {} Onset (quality: {})'.format(phase, quality)))
        self._manual_pick_plots.append(ax.plot([mpp - 0.5, mpp + 0.5],
                                               [y_bot, y_bot], linewidth=2,
                                               color=color))
        self._manual_pick_plots.append(ax.plot([mpp - 0.5, mpp + 0.5],
                                               [y_top, y_top], linewidth=2,
                                               color=color))
        ax.legend(loc=1)

    def fill_tabs(self, event=None, picked=False):
        self.clear_all()
        self.gen_pick_dlg()
        canvas_dict = self.parent().canvas_dict
        self.overview = self.gen_tab_widget('Overview', canvas_dict['mainFig'])
        id0 = self.figure_tabs.insertTab(0, self.pdlg_widget, 'Traces Plot')
        id1 = self.figure_tabs.insertTab(1, self.overview, 'Overview')
        id2 = self.figure_tabs.insertTab(2, self.p_tabs, 'P')
        id3 = self.figure_tabs.insertTab(3, self.s_tabs, 'S')
        if picked and self.get_current_event():
            self.fill_p_tabs(canvas_dict)
            self.fill_s_tabs(canvas_dict)
            self.toggle_autopickTabs(bool(self.fig_dict['mainFig'].axes))
            self.plot_manual_picks_to_figs()
        else:
            self.disable_autopickTabs()
        try:
            self.fig_dict['main_fig'].tight_layout()
        except:
            pass
        self.figure_tabs.setCurrentIndex(0)

    def fill_p_tabs(self, canvas_dict):
        for name in self.ptb_names:
            id = self.p_tabs.addTab(self.gen_tab_widget(name, canvas_dict[name]), name)
            self.p_tabs.setTabEnabled(id, bool(self.fig_dict[name].axes))
            try:
                self.fig_dict[name].tight_layout()
            except:
                pass

    def fill_s_tabs(self, canvas_dict):
        for name in self.stb_names:
            id = self.s_tabs.addTab(self.gen_tab_widget(name, canvas_dict[name]), name)
            self.s_tabs.setTabEnabled(id, bool(self.fig_dict[name].axes))
            try:
                self.fig_dict[name].tight_layout()
            except:
                pass

    def fill_figure_tabs(self):
        self.clear_all()
        self.p_tabs = QtWidgets.QTabWidget()
        self.s_tabs = QtWidgets.QTabWidget()
        self.tune_layout.insertWidget(0, self.figure_tabs)
        self.init_tab_names()

    def fill_eventbox(self):
        project = self.parent().project
        if not project:
            return
        # update own list
        self.parent().fill_eventbox(eventBox=self.eventBox, select_events='ref')
        index_start = self.parent().eventBox.currentIndex()
        index = index_start
        if index == -1:
            index += 1
        nevents = self.eventBox.model().rowCount()
        path = self.eventBox.itemText(index).split('*')[0]
        if project.getEventFromPath(path).isTestEvent():
            for index in range(nevents):
                path = self.eventBox.itemText(index)
                if project.getEventFromPath(index):
                    if not project.getEventFromPath(index).isTestEvent():
                        break
                # in case all events are marked as test events and last event is reached
                if index == nevents - 1:
                    index = -1
        self.eventBox.setCurrentIndex(index)
        if not index == index_start:
            self.eventBox.activated.emit(index)
        # update parent
        self.parent().fill_eventbox()

    def update_eventID(self):
        self.paraBox.boxes['eventID'].setText(
            self.get_current_event_name())
        self.figure_tabs.setCurrentIndex(0)

    def call_picker(self):
        self.parameter = self.params_from_gui()
        station = self.get_current_station()
        if not station:
            self._warn('No station selected')
            return
        wfpath = self.wftype if self.obspy_dmt else ''
        args = {'parameter': self.parameter,
                'station': station,
                'fnames': 'None',
                'eventid': [self.get_current_event_fp()],
                'iplot': 2,
                'fig_dict': self.fig_dict,
                'savexml': False,
                'obspyDMT_wfpath': wfpath}
        event = self.get_current_event()
        # self.parent().saveData(event, event.path, '.xml') MP MP uncommented because overwriting pick files in tune mode
        for key in self.fig_dict.keys():
            if not key == 'plot_style':
                self.fig_dict[key].clear()
        self.ap_thread = Thread(self, autoPyLoT, arg=args,
                                progressText='Picking trace...',
                                pb_widget=self.pb_widget,
                                redirect_stdout=True)
        self.enable(False)
        self.ap_thread.message.connect(self.add_log_item)
        self.ap_thread.finished.connect(self.finish_picker)
        self.figure_tabs.setCurrentIndex(4)
        self.ap_thread.start()
        # picks = autoPyLoT(self.parameter, fnames='None', iplot=2, fig_dict=self.fig_dict)

    def finish_picker(self):
        self.enable(True)
        if not self.ap_thread._executed:
            msg = 'Could not execute picker:\n{}'.format(
                self.ap_thread._executedError)
            info = self.ap_thread._executedErrorInfo
            self._warn(msg, info)
            return
        self.pylot_picks = self.ap_thread.data[self.get_current_event_name()]
        if not self.pylot_picks:
            self._warn('No picks found. See terminal output.')
            return
        # renew tabs
        # self.fill_figure_tabs()
        self.set_stretch()
        self.update.emit('Update')
        self.figure_tabs.setCurrentIndex(1)

    def enable(self, bool):
        self.pick_button.setEnabled(bool)
        self.paraBox.setEnabled(bool)
        self.eventBox.setEnabled(bool)
        self.stationBox.setEnabled(bool)
        self.overview.setEnabled(bool)
        self.p_tabs.setEnabled(bool)
        self.s_tabs.setEnabled(bool)

    def params_from_gui(self):
        parameters = self.paraBox.params_from_gui()
        if self.parent():
            self.parent()._inputs = parameters
        return parameters

    def set_stretch(self):
        self.tune_layout.setStretch(0, 3)
        self.tune_layout.setStretch(1, 1)

    def clear_all(self):
        if hasattr(self, 'pdlg_widget'):
            if self.pdlg_widget:
                self.pdlg_widget.setParent(None)
                del self.pdlg_widget
        if hasattr(self, 'overview'):
            self.overview.setParent(None)
        if hasattr(self, 'p_tabs'):
            self.p_tabs.clear()
            self.p_tabs.setParent(None)
        if hasattr(self, 's_tabs'):
            self.s_tabs.clear()
            self.s_tabs.setParent(None)

    def disable_autopickTabs(self):
        self.toggle_autopickTabs(False)

    def toggle_autopickTabs(self, bool):
        self.figure_tabs.setTabEnabled(1, bool)
        self.figure_tabs.setTabEnabled(2, bool)
        self.figure_tabs.setTabEnabled(3, bool)

    def _warn(self, message, info=None):
        self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Warning,
                                         'Warning', message)
        self.qmb.setDetailedText(str(info))
        self.qmb.show()


class PylotParameterWidget(QtWidgets.QWidget):
    accepted = QtCore.Signal(str)
    rejected = QtCore.Signal(str)

    def __init__(self, parameter, parent=None, windowflag=Qt.Window):
        '''
        Generate Widget containing parameters for PyLoT.
        :param: parameter
        :type: PylotParameter (object)
        '''
        QtWidgets.QWidget.__init__(self, parent, windowflag)
        self.parameter = parameter
        self.tabs = QtWidgets.QTabWidget()
        self.layout = QtWidgets.QVBoxLayout()
        self._init_save_buttons()
        self._init_tabs()
        self._init_dialog_buttons()
        self.labels = {}
        self.boxes = {}
        self.groupboxes = {}
        self._exclusive_widgets = []
        self._init_sublayouts()
        self.setLayout(self.layout)
        self.add_main_parameters_tab()
        self.add_special_pick_parameters_tab()
        self.params_to_gui()
        self._toggle_advanced_settings()
        self.resize(720, 860)
        self.center()
        self.setWindowModality(QtCore.Qt.WindowModality.ApplicationModal)
        self.accepted.connect(self.params_from_gui)
        self.rejected.connect(self.params_to_gui)

    def center(self):
        fm = self.frameGeometry()
        screen = QtWidgets.QApplication.desktop().screenNumber(QtWidgets.QApplication.desktop().cursor().pos())
        centerPoint = QtWidgets.QApplication.desktop().screenGeometry(screen).center()
        fm.moveCenter(centerPoint)
        self.move(fm.topLeft())

    def _init_sublayouts(self):
        self._main_layout = QtWidgets.QVBoxLayout()
        self._advanced_layout = QtWidgets.QVBoxLayout()
        self._create_advanced_cb()

    def _init_save_buttons(self):
        self._buttons_layout = QtWidgets.QHBoxLayout()
        self.loadButton = QtWidgets.QPushButton('&Load settings')
        self.saveButton = QtWidgets.QPushButton('&Save settings')
        self.defaultsButton = QtWidgets.QPushButton('&Defaults')
        self._buttons_layout.addWidget(self.loadButton)
        self._buttons_layout.addWidget(self.saveButton)
        self._buttons_layout.addWidget(self.defaultsButton)
        self.layout.addLayout(self._buttons_layout)
        self.loadButton.clicked.connect(self.openFile)
        self.saveButton.clicked.connect(self.saveFile)
        self.defaultsButton.clicked.connect(self.restoreDefaults)

    def _init_tabs(self):
        self.layout.addWidget(self.tabs)

    def _init_dialog_buttons(self):
        self._dialog_buttons = QtWidgets.QHBoxLayout()
        self._okay = QtWidgets.QPushButton('Ok')
        self._close = QtWidgets.QPushButton('Close')
        self._apply = QtWidgets.QPushButton('Apply')
        self._dialog_buttons.addWidget(self._okay)
        self._dialog_buttons.addWidget(self._close)
        self._dialog_buttons.addWidget(self._apply)
        self._okay.clicked.connect(self.accept)
        self._okay.clicked.connect(self.close)
        self._apply.clicked.connect(self.accept)
        self._close.clicked.connect(self.close)
        self.layout.addLayout(self._dialog_buttons)

    def _create_advanced_cb(self):
        self._advanced_cb = QtWidgets.QCheckBox('Enable Advanced Settings')
        self._advanced_layout.insertWidget(0, self._advanced_cb)
        self._advanced_cb.toggled.connect(self._toggle_advanced_settings)

    def _toggle_advanced_settings(self):
        if self._advanced_cb.isChecked():
            self._enable_advanced(True)
        else:
            self._enable_advanced(False)

    def _enable_advanced(self, enable):
        for lst in self.parameter.get_special_para_names().values():
            for param in lst:
                box = self.boxes[param]
                if type(box) is not list:
                    box.setEnabled(enable)
                else:
                    for b in box:
                        b.setEnabled(enable)

    def set_tune_mode(self, bool):
        names = ['Directories', 'NLLoc',
                 'Seismic Moment']
        for name in names:
            self.hide_groupbox(name)
        if bool:
            self._apply.hide()
            self._okay.hide()
            self._close.hide()
        else:
            self._apply.show()
            self._okay.show()
            self._close.show()

    def init_boxes(self, parameter_names):
        grid = QtWidgets.QGridLayout()

        for index1, name in enumerate(parameter_names):
            if name in ['rootpath', 'database']:
                logging.warning(
                    f'Deprecated parameter loaded: {name}. Check if datapath is still correct in parameter widget.'
                )
                continue
            default_item = self.parameter.get_defaults()[name]
            tooltip = default_item['tooltip']
            tooltip += ' | type: {}'.format(default_item['type'])
            if not type(default_item['type']) == tuple:
                typ = default_item['type']
                box = self.create_box(typ, tooltip)
                self.boxes[name] = box
                namestring = default_item['namestring']
            elif type(default_item['type']) == tuple:
                boxes = []
                values = self.parameter[name]
                for index2, val in enumerate(values):
                    typ = default_item['type'][index2]
                    boxes.append(self.create_box(typ, tooltip))
                headline = default_item['namestring'][1:]
                box, lower = self.create_multi_box(boxes, headline)
                self.boxes[name] = boxes
                namestring = default_item['namestring'][0]
            text = namestring + ' [?]'
            label = QtWidgets.QLabel(text)
            self.labels[name] = label
            label.setToolTip(tooltip)
            grid.addWidget(label, index1, 1)
            grid.addWidget(box, index1, 2)
        return grid

    def create_box(self, typ, tooltip):
        if typ == str:
            box = QtWidgets.QLineEdit()
        elif typ == float:
            box = QDoubleSpinBox()
            box.setDecimals(4)
            box.setRange(-10e4, 10e4)
        elif typ == int:
            box = QSpinBox()
        elif typ == bool:
            box = QtWidgets.QCheckBox()
        else:
            raise TypeError('Unrecognized type {}'.format(typ))
        return box

    def create_multi_box(self, boxes, headline=None):
        box = QtWidgets.QWidget()
        gl = QtWidgets.QGridLayout()
        column = 0
        if headline:
            for index, item in enumerate(headline):
                if not item:
                    continue
                gl.addWidget(QtWidgets.QLabel(item), index, 0, Qt.Alignment(2))
                column = 1
        for index, b in enumerate(boxes):
            gl.addWidget(b, index, column)
        box.setLayout(gl)
        return box, column

    def add_tab(self, layout, name):
        widget = QtWidgets.QWidget()
        scrollA = QtWidgets.QScrollArea()
        scrollA.setWidgetResizable(True)
        scrollA.setWidget(widget)
        widget.setLayout(layout)
        self.tabs.addTab(scrollA, name)

    def add_main_parameters_tab(self):
        self.add_to_layout(self._main_layout, 'Directories',
                           self.parameter.get_main_para_names()['dirs'], 0)
        self.add_to_layout(self._main_layout, 'NLLoc',
                           self.parameter.get_main_para_names()['nlloc'], 1)
        self.add_to_layout(self._main_layout, 'Seismic Moment',
                           self.parameter.get_main_para_names()['smoment'], 2)
        self.add_to_layout(self._main_layout, 'Local Magnitude',
                           self.parameter.get_main_para_names()['localmag'], 3)
        self.add_to_layout(self._main_layout, 'Filter Settings',
                           self.parameter.get_main_para_names()['filter'], 4)
        self.add_to_layout(self._main_layout, 'Common Settings Characteristic Function',
                           self.parameter.get_main_para_names()['pick'], 5)
        self.add_tab(self._main_layout, 'Main Settings')

    def add_special_pick_parameters_tab(self):
        self.add_to_layout(self._advanced_layout, 'Z-component',
                           self.parameter.get_special_para_names()['z'], 1)
        self.add_to_layout(self._advanced_layout, 'H-components',
                           self.parameter.get_special_para_names()['h'], 2)
        self.add_to_layout(self._advanced_layout, 'First-motion picker',
                           self.parameter.get_special_para_names()['fm'], 3)
        self.add_to_layout(self._advanced_layout, 'Quality assessment',
                           self.parameter.get_special_para_names()['quality'], 4)
        self.add_tab(self._advanced_layout, 'Advanced Settings')

    # def gen_h_separator(self):
    #     separator = QtGui.QFrame()
    #     separator.setFrameShape(QtGui.QFrame.HLine)
    #     return separator

    # def gen_headline(self, text):
    #     label=QtGui.QLabel(text)
    #     font=QtGui.QFont()
    #     font.setBold(True)
    #     label.setFont(font)
    #     return label

    def refresh(self):
        for groupbox in self.groupboxes.values():
            layout = groupbox._parentLayout
            position = groupbox._position
            layout.insertWidget(position, groupbox)

    def get_groupbox_exclusive(self, name):
        widget = QtWidgets.QWidget(self, Qt.WindowFlags(1))
        layout = QtWidgets.QVBoxLayout()
        widget.setLayout(layout)
        layout.addWidget(self.groupboxes[name])
        self._exclusive_widgets.append(widget)
        return widget

    def get_groupbox_dialog(self, name):
        widget = self.get_groupbox_exclusive(name)
        dialog = QtWidgets.QDialog(self.parent())
        layout = QtWidgets.QVBoxLayout()
        dialog.setLayout(layout)
        buttonbox = QtWidgets.QDialogButtonBox(QDialogButtonBox.Ok |
                                               QDialogButtonBox.Cancel)
        buttonbox.accepted.connect(dialog.accept)
        buttonbox.accepted.connect(self.refresh)
        buttonbox.accepted.connect(self.params_from_gui)
        buttonbox.rejected.connect(dialog.reject)
        buttonbox.rejected.connect(self.refresh)
        buttonbox.rejected.connect(self.params_to_gui)
        layout.addWidget(widget)
        layout.addWidget(buttonbox)
        self._exclusive_dialog = dialog
        return dialog

    def add_to_layout(self, layout, name, items, position):
        groupbox = QtWidgets.QGroupBox(name)
        groupbox._position = position
        groupbox._parentLayout = layout
        self.groupboxes[name] = groupbox
        groupbox.setLayout(self.init_boxes(items))
        layout.insertWidget(position, groupbox)

    def show_groupboxes(self):
        for name in self.groupboxes.keys():
            self.show_groupbox(name)
        self._advanced_cb.show()

    def hide_groupboxes(self):
        for name in self.groupboxes.keys():
            self.hide_groupbox(name)
        self._advanced_cb.hide()

    def show_groupbox(self, name):
        if name in self.groupboxes.keys():
            self.groupboxes[name].show()
        else:
            print('Groupbox {} not part of object.'.format(name))

    def hide_groupbox(self, name):
        if name in self.groupboxes.keys():
            self.groupboxes[name].hide()
        else:
            print('Groupbox {} not part of object.'.format(name))

    def show_file_buttons(self):
        self.saveButton.show()
        self.loadButton.show()
        self.defaultsButton.show()

    def hide_file_buttons(self):
        self.saveButton.hide()
        self.loadButton.hide()
        self.defaultsButton.hide()

    def show_parameter(self, name=None):
        if not name:
            for name in self.boxes.keys():
                self.show_parameter(name)
            return
        if name in self.boxes.keys() and name in self.labels.keys():
            # comprising case type(self.boxes[name]) == list
            boxes = self.boxes[name]
            if not type(boxes) == list:
                boxes = [boxes]
            for box in boxes:
                box.show()
            self.labels[name].show()
        else:
            print('Parameter {} not part of object.'.format(name))

    def hide_parameter(self, name=None):
        if not name:
            for name in self.boxes.keys():
                self.hide_parameter(name)
            return
        if name in self.boxes.keys() and name in self.labels.keys():
            # comprising case type(self.boxes[name]) == list
            boxes = self.boxes[name]
            if not type(boxes) == list:
                boxes = [boxes]
            for box in boxes:
                box.hide()
            self.labels[name].hide()
        else:
            print('Parameter {} not part of object.'.format(name))

    def params_from_gui(self):
        for param in self.parameter.get_all_para_names():
            box = self.boxes[param]
            value = self.getValue(box)
            self.parameter.checkValue(param, value)
            self.parameter.setParamKV(param, value)
        return self.parameter

    def params_to_gui(self, tuneMode=False):
        for param in self.parameter.get_all_para_names():
            if param == 'eventID':
                if tuneMode:
                    continue
            box = self.boxes[param]
            value = self.parameter[param]
            # self.parameter.checkValue(param, value)
            self.setValue(box, value)

    def setValue(self, box, value):
        if type(box) == QtWidgets.QLineEdit:
            box.setText(str(value))
        elif type(box) == QSpinBox or type(box) == QDoubleSpinBox:
            if not value:
                value = 0.
            box.setValue(value)
        elif type(box) == QtWidgets.QCheckBox:
            if value == 'True':
                value = True
            if value == 'False' or value is None:
                value = False
            box.setChecked(value)
        elif type(box) == list:
            for index, b in enumerate(box):
                self.setValue(b, value[index])

    def getValue(self, box):
        if type(box) == QtWidgets.QLineEdit:
            value = str(box.text())
        elif type(box) == QSpinBox or type(box) == QDoubleSpinBox:
            value = box.value()
        elif type(box) == QtWidgets.QCheckBox:
            value = box.isChecked()
        elif type(box) == list:
            value = []
            for b in box:
                value.append(self.getValue(b))
            value = tuple(value)
        return value

    def openFile(self):
        fd = QtWidgets.QFileDialog()
        fname = fd.getOpenFileName(self, 'Browse for settings file.',
                                   filter='PyLoT input file (*.in)')
        if fname[0]:
            try:
                self.parameter.from_file(fname[0])
                self.params_to_gui(tuneMode=True)
            except Exception as e:
                self._warn('Could not open file {}:\n{}'.format(fname[0], e))
                return

    def saveFile(self):
        fd = QtWidgets.QFileDialog()
        fname = fd.getSaveFileName(self, 'Browse for settings file.',
                                   filter='PyLoT input file (*.in)')[0]
        if fname:
            if not fname.endswith('.in'):
                fname += '.in'
            try:
                self.params_from_gui()
                self.parameter.export2File(fname)
            except Exception as e:
                self._warn('Could not save file {}:\n{}'.format(fname, e))
                return

    def restoreDefaults(self):
        try:
            self.parameter.reset_defaults()
            self.params_to_gui(tuneMode=True)
        except Exception as e:
            self._warn('Could not restore defaults:\n{}'.format(e))
            return

    def show(self):
        self.refresh()
        self.show_parameter()
        if hasattr(self, '_exclusive_dialog'):
            self._exclusive_dialog.close()
        self._exclusive_widgets = []
        QtWidgets.QWidget.show(self)

    def close(self):
        self.rejected.emit('reject')
        QtWidgets.QWidget.close(self)

    def accept(self):
        self.accepted.emit('accept')

    def _warn(self, message):
        self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Warning,
                                         'Warning', message)
        self.qmb.show()


class AutoPickDlg(QDialog):
    def __init__(self, parent=None, sge=False):
        super(AutoPickDlg, self).__init__(parent)

        if not self.parent():
            print('WARNING: No parent given! PyLoT parameter can not be read!')

        self.pp_export = os.path.join(os.path.expanduser('~'), '.pylot', '.pylot_exported.in')

        self.sge = sge
        self.layout = QVBoxLayout()
        self.setWindowTitle('Pick whole project...')
        self.setLayout(self.layout)

        self.setupUI()
        self.resize(720, 0)

    def setupUI(self):
        self.init_sb_ncores()
        self.init_gb()

        self.layout.addLayout(self.layout_ncores, 2)
        self.layout.addWidget(self.gb)
        self.addJobWidget()

    def init_sb_ncores(self):
        self.layout_ncores = QHBoxLayout()
        self.label_ncores = QLabel('Number of CPU cores:')
        self.sb_ncores = QSpinBox()
        if not self.sge:
            maxcores = multiprocessing.cpu_count()
            self.sb_ncores.setMinimum(1)
            self.sb_ncores.setMaximum(maxcores)
            self.layout_ncores.addWidget(self.label_ncores)
            self.layout_ncores.addWidget(self.sb_ncores)
            self.layout_ncores.setStretch(0, 1)
            self.layout_ncores.setStretch(1, 0)

    def init_gb(self):
        title_sge = {
            True: 'Submit autoPyLoT process to grid engine',
            False: 'Spawn autoPyLot process on local machine'
        }

        self.gb = QGroupBox(title_sge[self.sge])

    def addJobWidget(self):
        widget_sge = {
            True: Submit2Grid(),
            False: SubmitLocal()
        }

        self.job_widget = widget_sge[self.sge]
        self.job_widget.button.clicked.connect(self.accept)

        self.jobLayout = QVBoxLayout()
        self.jobLayout.addWidget(self.job_widget)

        self.gb.setLayout(self.jobLayout)

    def exportParameter(self):
        self.parent().exportEvents()
        pylot_params = self.parent()._inputs
        self.addEvents2pp(pylot_params)
        pylot_params.export2File(self.pp_export)

    def addEvents2pp(self, pylot_parameter):
        eventIDs = []
        for event in self.parent().project.eventlist:
            eventIDs.append(event.pylot_id)
        pylot_parameter['eventID'] = eventIDs

    def accept(self):
        self.exportParameter()
        self.job_widget.start(self.pp_export, self.sb_ncores.value())
        QDialog.accept(self)


class Submit2Grid(QWidget):
    def __init__(self, parent=None):
        super(Submit2Grid, self).__init__(parent)
        self.main_layout = QVBoxLayout()
        self.sub_layout = QVBoxLayout()
        self.setLayout(self.main_layout)
        self.label = QLabel('Full GE command (without script name):')
        self.textedit = QLineEdit()
        self.button = QPushButton('Run')

        self.script_fn = '.autoPyLot.sh'

        self.sub_layout.addWidget(self.label)
        self.sub_layout.addWidget(self.textedit)

        self.main_layout.addLayout(self.sub_layout)
        self.main_layout.addWidget(self.button)

        self.setDefaultCommand()

    def setDefaultCommand(self):
        default_command = 'qsub -l low -cwd -q TARGET_MACHINE -pe mpi-fu NCORES'
        self.textedit.setText(default_command)

    def start(self, pp_export, ncores=None):
        self.genShellScript(pp_export)
        self.execute_script()

    def genShellScript(self, pp_export):
        outfile = open(self.script_fn, 'w')
        outfile.write('#!/bin/sh\n\n')
        try:
            ncores = int(self.textedit.text().split()[-1])
            ncores = '--ncores {}'.format(ncores)
        except:
            ncores = None
        outfile.write('python autoPyLoT.py -i {} {}\n'.format(pp_export, ncores))
        outfile.close()

    def execute_script(self):
        command = self.textedit.text().strip().split(' ')
        command.append(self.script_fn)
        p = subprocess.Popen(command)
        cmd_str = str()
        for item in command:
            cmd_str += item + ' '
        print('exec. command: {}'.format(cmd_str))
        print('Spawned autoPyLoT process with pid {}'.format(p.pid))


class SubmitLocal(QWidget):
    def __init__(self, parent=None):
        super(SubmitLocal, self).__init__(parent)
        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)
        self.button = QPushButton('Run')

        self.script_fn = ['python', 'autoPyLoT.py', '-i']

        self.main_layout.addWidget(self.button)

    def start(self, pp_export, ncores):
        self.execute_command(pp_export, ncores)

    def execute_command(self, pp_export, ncores):
        command = self.script_fn[:]
        command.append(pp_export)
        command.append('--ncores')
        command.append(str(ncores))
        cmd_str = str()
        for item in command:
            cmd_str += item + ' '
        print('exec. command: {}'.format(cmd_str))
        p = subprocess.Popen(command)
        print('Spawned autoPyLoT process with pid {}'.format(p.pid))


class PropertiesDlg(QDialog):
    def __init__(self, parent=None, infile=None, inputs=None):
        super(PropertiesDlg, self).__init__(parent)
        self._pylot_mainwindow = self.parent()

        self.infile = infile
        self.inputs = inputs
        self.dirty = False

        self.setWindowTitle("PyLoT Properties")
        self.tabWidget = QTabWidget()
        self.tabWidget.addTab(InputsTab(self), "Main")
        # self.tabWidget.addTab(OutputsTab(self), "Outputs")
        self.tabWidget.addTab(PhasesTab(self, inputs), "Phases")
        self.tabWidget.addTab(GraphicsTab(self), "Graphics")
        # self.tabWidget.addTab(LocalisationTab(self), "Loc. Tools")
        self.tabWidget.addTab(LocalisationTab(self), "NonLinLoc")
        self.tabWidget.addTab(ChannelOrderTab(self), "Channel Order")
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok |
                                          QDialogButtonBox.Apply |
                                          QDialogButtonBox.Close |
                                          QDialogButtonBox.RestoreDefaults)

        layout = QVBoxLayout()
        layout.addWidget(self.tabWidget)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)
        self.setFixedWidth(700)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(self.apply)
        self.buttonBox.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.restore)

    def setDirty(self, level=1):
        if type(self.dirty) == int:
            if level < self.dirty:
                return
        self.dirty = level

    def getinfile(self):
        return self.infile

    def accept(self, *args, **kwargs):
        self.apply()
        QDialog.accept(self)

    def apply(self):
        for widint in range(self.tabWidget.count()):
            curwid = self.tabWidget.widget(widint)
            if curwid.dirty:
                self.setDirty(curwid.dirty)
            values = curwid.getValues()
            if values is not None:
                self.setValues(values)

    def close(self):
        self.reset_current()
        QDialog.close(self)

    def show(self):
        self.keep_current()
        QDialog.show(self)

    def restore(self):
        for widint in range(self.tabWidget.count()):
            curwid = self.tabWidget.widget(widint)
            values = curwid.resetValues(self.getinfile())
            if values is not None:
                self.setValues(values)

    def keep_current(self):
        self._current_values = []
        for widint in range(self.tabWidget.count()):
            curwid = self.tabWidget.widget(widint)
            values = curwid.getValues()
            if values is not None:
                self._current_values.append(values)

    def reset_current(self):
        for values in self._current_values:
            self.setValues(values)

    @staticmethod
    def setValues(tabValues):
        settings = QSettings()
        compclass = SetChannelComponents.from_qsettings(settings)

        for setting, value in tabValues.items():
            settings.setValue(setting, value)
            if value is not None:
                if setting in ['Z', 'N', 'E']:
                    compclass.setCompPosition(value, setting, False)

        settings.sync()


class PropTab(QWidget):
    def __init__(self, parent=None):
        super(PropTab, self).__init__(parent)

    def getValues(self):
        return None

    def resetValues(self, infile):
        return None


class InputsTab(PropTab):
    def __init__(self, parent, infile=None):
        super(InputsTab, self).__init__(parent)

        self.dirty = False
        settings = QSettings()
        pylot_user = getpass.getuser()
        fulluser = settings.value("user/FullName")
        login = settings.value("user/Login")

        fullNameLabel = QLabel("Full name for user '{0}': ".format(pylot_user))

        # get the full name of the actual user
        self.fullNameEdit = QLineEdit()
        try:
            self.fullNameEdit.setText(fulluser)
        except TypeError as e:
            self.fullNameEdit.setText(fulluser[0])
        # TODO: check settings and substitute datapath for root
        # information about data structure
        dataroot = settings.value("data/dataRoot")
        curstructure = settings.value("data/Structure")
        self.dataDirEdit = QLineEdit()
        self.dataDirEdit.setText(dataroot)
        self.dataDirEdit.selectAll()
        self.structureSelect = QComboBox()
        self.cutLabel = QLabel
        self.cuttimesLayout = QHBoxLayout()
        self.autosaveXML_checkbox = QCheckBox()
        self.tstartBox = QSpinBox()
        self.tstopBox = QSpinBox()
        for spinbox in [self.tstartBox, self.tstopBox]:
            spinbox.setRange(-99999, 99999)
        self.tstartBox.setValue(float(settings.value('tstart')) if get_none(settings.value('tstart')) else 0)
        self.tstopBox.setValue(float(settings.value('tstop')) if get_none(settings.value('tstop')) else 1e6)
        self.cuttimesLayout.addWidget(self.tstartBox, 10)
        self.cuttimesLayout.addWidget(QLabel('[s] and'), 0)
        self.cuttimesLayout.addWidget(self.tstopBox, 10)
        self.cuttimesLayout.addWidget(QLabel('[s]'), 0)

        from pylot.core.util.structure import DATASTRUCTURE

        self.structureSelect.addItems(list(DATASTRUCTURE.keys()))

        dsind = findComboBoxIndex(self.structureSelect, curstructure)

        self.structureSelect.setCurrentIndex(dsind)

        layout = QtWidgets.QFormLayout()
        layout.addRow("Data root directory: ", self.dataDirEdit)
        layout.addRow("Full name for user '{0}': ".format(pylot_user), self.fullNameEdit)
        layout.addRow("Data structure: ", self.structureSelect)
        layout.addRow('Cut waveforms between (relative to source origin time if given): ', self.cuttimesLayout)
        layout.addRow('Automatically save events as XML when saving Project: ', self.autosaveXML_checkbox)

        self.setLayout(layout)
        self.connectSignals()

    def connectSignals(self):
        self.tstartBox.valueChanged.connect(self.checkSpinboxStop)
        self.tstopBox.valueChanged.connect(self.checkSpinboxStart)
        self.tstartBox.valueChanged.connect(self.setDirty)
        self.tstopBox.valueChanged.connect(self.setDirty)
        self.autosaveXML_checkbox.stateChanged.connect(self.setDirty)

    def setDirty(self):
        self.dirty = 2

    def checkSpinboxStop(self):
        if self.tstopBox.value() < self.tstartBox.value():
            self.tstopBox.setValue(self.tstartBox.value())

    def checkSpinboxStart(self):
        if self.tstopBox.value() < self.tstartBox.value():
            self.tstartBox.setValue(self.tstopBox.value())

    def getValues(self):
        values = {"data/dataRoot": self.dataDirEdit.text(),
                  "user/FullName": self.fullNameEdit.text(),
                  "data/Structure": self.structureSelect.currentText(),
                  "tstart": self.tstartBox.value(),
                  "tstop": self.tstopBox.value(),
                  "autosaveXML": self.autosaveXML_checkbox.isChecked()}
        return values

    def resetValues(self, infile):
        para = PylotParameter(infile)
        datstruct = para.get('datastructure')
        if datstruct == 'SeisComp':
            index = 0
        else:
            index = 2
        datapath = para.get('datapath') if not para.get('datapath') is None else ''
        values = {"data/dataRoot": self.dataDirEdit.setText("%s" % datapath),
                  "user/FullName": self.fullNameEdit.text(),
                  "data/Structure": self.structureSelect.setCurrentIndex(index),
                  "tstart": self.tstartBox.setValue(0),
                  "tstop": self.tstopBox.setValue(1e6),
                  "autosaveXML": self.autosaveXML_checkbox.setChecked(True), }
        return values


class OutputsTab(PropTab):
    def __init__(self, parent=None, infile=None):
        super(OutputsTab, self).__init__(parent)

        self.dirty = False
        settings = QSettings()
        curval = settings.value("output/Format", None)

        eventOutputLabel = QLabel("event/picks output format")
        self.eventOutputComboBox = QComboBox()
        eventoutputformats = OUTPUTFORMATS.keys()
        self.eventOutputComboBox.addItems(eventoutputformats)

        ind = findComboBoxIndex(self.eventOutputComboBox, curval)

        self.eventOutputComboBox.setCurrentIndex(ind)
        layout = QGridLayout()
        layout.addWidget(eventOutputLabel, 0, 0)
        layout.addWidget(self.eventOutputComboBox, 0, 1)

        self.setLayout(layout)

    def getValues(self):
        values = {"output/Format": self.eventOutputComboBox.currentText()}
        return values

    def resetValues(self, infile):
        values = {"output/Format": self.eventOutputComboBox.setCurrentIndex(1)}
        return values


class PhasesTab(PropTab):
    def __init__(self, parent=None, inputs=None):
        super(PhasesTab, self).__init__(parent)
        self.inputs = inputs
        self.dirty = False

        self.Pphases = 'P, Pg, Pn, PmP, P1, P2, P3'
        self.Sphases = 'S, Sg, Sn, SmS, S1, S2, S3'

        self.PphasesEdit = QLineEdit()
        self.SphasesEdit = QLineEdit()

        self.pickDefaultsButton = QtWidgets.QPushButton('Choose default phases...')
        PphasesLabel = QLabel("P Phases to pick")
        SphasesLabel = QLabel("S Phases to pick")
        notes_label = QLabel(
            'Note: Only for manual picking. '
            'Defaults depend on parameter "Array extent"')

        settings = QSettings()
        Pphases = settings.value('p_phases')
        Sphases = settings.value('s_phases')

        self.PphasesEdit.setText("%s" % Pphases)
        self.SphasesEdit.setText("%s" % Sphases)

        self.main_layout = QtWidgets.QHBoxLayout()
        layout = QGridLayout()
        layout.addWidget(PphasesLabel, 0, 0)
        layout.addWidget(SphasesLabel, 1, 0)

        layout.addWidget(self.PphasesEdit, 0, 1)
        layout.addWidget(self.SphasesEdit, 1, 1)

        layout.addWidget(self.pickDefaultsButton, 2, 0)
        layout.addWidget(notes_label, 2, 1)
        self.main_layout.addLayout(layout)
        self.setLayout(self.main_layout)

        self.connectSignals()

    def connectSignals(self):
        self.pickDefaultsButton.clicked.connect(self.get_defaults)

    def get_defaults(self):
        phases = [p.strip() for p in self.Pphases.split(',')] + [s.strip() for s in self.Sphases.split(',')]
        p_current = [p.strip() for p in self.PphasesEdit.text().split(',')]
        s_current = [s.strip() for s in self.SphasesEdit.text().split(',')]

        current_phases = p_current + s_current
        if self.inputs:
            parameter = self.inputs
            if parameter.get('extent') == 'global':
                # get all default phase names known to obspy.taup
                # in a list and remove duplicates
                phases = list(set(get_phase_names('ttall')))
                phases.sort()

        phaseDefaults = PhaseDefaults(self, phase_defaults=phases,
                                      current_phases=current_phases)
        if phaseDefaults.exec_():
            phase_dict = self.sortPhases(phaseDefaults.selected_phases)
            p_phases = ''
            s_phases = ''
            for index, p in enumerate(phase_dict['P']):
                p_phases += p
                if not index == len(phase_dict['P']) - 1:
                    p_phases += ', '
            for index, s in enumerate(phase_dict['S']):
                s_phases += s
                if not index == len(phase_dict['S']) - 1:
                    s_phases += ', '
            self.PphasesEdit.setText(p_phases)
            self.SphasesEdit.setText(s_phases)

    @staticmethod
    def sortPhases(phases):
        sorted_phases = {'P': [],
                         'S': []}
        for phase in phases:
            idf_phase = loopIdentifyPhase(phase)
            if idf_phase:
                sorted_phases[identifyPhase(idf_phase)].append(phase)
        return sorted_phases

    def getValues(self):
        p_phases = self.PphasesEdit.text()
        s_phases = self.SphasesEdit.text()
        values = {'p_phases': p_phases,
                  's_phases': s_phases}
        return values

    def resetValues(self, infile=None):
        values = {'p_phases': self.PphasesEdit.setText(self.Pphases),
                  's_phases': self.SphasesEdit.setText(self.Sphases)}
        return values


class GraphicsTab(PropTab):
    def __init__(self, parent=None):
        super(GraphicsTab, self).__init__(parent)
        self.dirty = False
        self.pylot_mainwindow = parent._pylot_mainwindow
        self.init_layout()
        # self.add_pg_cb()
        self.add_nth_sample()
        self.add_style_settings()
        self.setLayout(self.main_layout)

    def init_layout(self):
        self.main_layout = QGridLayout()

    def add_style_settings(self):
        styles = self.pylot_mainwindow._styles
        active_stylename = self.pylot_mainwindow._stylename
        label = QtWidgets.QLabel('Application style (might require Application restart):')
        self.style_cb = QComboBox()
        for stylename, style in styles.items():
            self.style_cb.addItem(stylename, style)
        index_current_style = self.style_cb.findText(active_stylename)
        self.style_cb.setCurrentIndex(index_current_style)
        self.main_layout.addWidget(label, 2, 0)
        self.main_layout.addWidget(self.style_cb, 2, 1)
        self.style_cb.activated.connect(self.set_current_style)

    def add_nth_sample(self):
        settings = QSettings()
        nth_sample = settings.value("nth_sample")
        if not nth_sample:
            nth_sample = 1

        self.spinbox_nth_sample = QSpinBox()
        label = QLabel('nth sample')
        label.setToolTip('Plot every nth sample (to speed up plotting)')
        self.spinbox_nth_sample.setMinimum(1)
        self.spinbox_nth_sample.setMaximum(10e3)
        self.spinbox_nth_sample.setValue(int(nth_sample))
        self.main_layout.addWidget(label, 1, 0)
        self.main_layout.addWidget(self.spinbox_nth_sample, 1, 1)
        self.spinbox_nth_sample.valueChanged.connect(self.setDirty)

    def setDirty(self):
        self.dirty = 1

    def set_current_style(self):
        selected_style = self.style_cb.currentText()
        self.pylot_mainwindow.set_style(selected_style)

    def getValues(self):
        values = {'nth_sample': self.spinbox_nth_sample.value()}
        return values

    def resetValues(self, infile=None):
        values = {'nth_sample': self.spinbox_nth_sample.setValue(1)}
        return values


class ChannelOrderTab(PropTab):
    def __init__(self, parent=None, infile=None):
        super(ChannelOrderTab, self).__init__(parent)

        self.dirty = False
        settings = QSettings()
        compclass = SetChannelComponents.from_qsettings(settings)

        ChannelOrderLabelZ = QLabel("Channel Z [up/down, default=3]")
        ChannelOrderLabelN = QLabel("Channel N [north/south, default=1]")
        ChannelOrderLabelE = QLabel("Channel E [east/west, default=2]")
        self.ChannelOrderZEdit = QLineEdit()
        self.ChannelOrderZEdit.setMaxLength(1)
        self.ChannelOrderZEdit.setFixedSize(20, 20)
        self.ChannelOrderNEdit = QLineEdit()
        self.ChannelOrderNEdit.setMaxLength(1)
        self.ChannelOrderNEdit.setFixedSize(20, 20)
        self.ChannelOrderEEdit = QLineEdit()
        self.ChannelOrderEEdit.setMaxLength(1)
        self.ChannelOrderEEdit.setFixedSize(20, 20)
        # get channel order settings
        zcomp = compclass.getCompPosition('Z')
        ncomp = compclass.getCompPosition('N')
        ecomp = compclass.getCompPosition('E')
        self.ChannelOrderZEdit.setText("%s" % zcomp)
        self.ChannelOrderNEdit.setText("%s" % ncomp)
        self.ChannelOrderEEdit.setText("%s" % ecomp)

        layout = QGridLayout()
        layout.addWidget(ChannelOrderLabelZ, 0, 0)
        layout.addWidget(ChannelOrderLabelN, 1, 0)
        layout.addWidget(ChannelOrderLabelE, 2, 0)
        layout.addWidget(self.ChannelOrderZEdit, 0, 1)
        layout.addWidget(self.ChannelOrderNEdit, 1, 1)
        layout.addWidget(self.ChannelOrderEEdit, 2, 1)

        self.setLayout(layout)
        self.connectSignals()

    def connectSignals(self):
        self.ChannelOrderZEdit.textEdited.connect(self.checkDoubleZ)
        self.ChannelOrderNEdit.textEdited.connect(self.checkDoubleN)
        self.ChannelOrderEEdit.textEdited.connect(self.checkDoubleE)
        self.ChannelOrderZEdit.textEdited.connect(self.setDirty)
        self.ChannelOrderNEdit.textEdited.connect(self.setDirty)
        self.ChannelOrderEEdit.textEdited.connect(self.setDirty)

    def setDirty(self):
        self.dirty = 1

    def checkDoubleZ(self, text):
        self.checkDouble(text, 'Z')

    def checkDoubleN(self, text):
        self.checkDouble(text, 'N')

    def checkDoubleE(self, text):
        self.checkDouble(text, 'E')

    def checkDouble(self, text, comp):
        channelOrderEdits = {
            'Z': self.ChannelOrderZEdit,
            'N': self.ChannelOrderNEdit,
            'E': self.ChannelOrderEEdit
        }
        for key in channelOrderEdits.keys():
            if key == comp:
                continue
            if str(channelOrderEdits[key].text()) == str(text):
                channelOrderEdits[key].setText('')

    def getValues(self):
        values = {"Z": int(self.ChannelOrderZEdit.text()),
                  "N": int(self.ChannelOrderNEdit.text()),
                  "E": int(self.ChannelOrderEEdit.text())}
        return values

    def resetValues(self, infile=None):
        Zdefault = 3
        Ndefault = 1
        Edefault = 2
        values = {"Z": self.ChannelOrderZEdit.setText("%d" % Zdefault),
                  "N": self.ChannelOrderNEdit.setText("%d" % Ndefault),
                  "E": self.ChannelOrderEEdit.setText("%d" % Edefault)}
        return values

        # MP MP: No idea why this function exists!?
        # def getComponents(self):
        #     self.CompName = dict(Z='10', N='11', E='12')


class LocalisationTab(PropTab):
    def __init__(self, parent=None, infile=None):
        super(LocalisationTab, self).__init__(parent)

        self.dirty = False
        settings = QSettings()
        curtool = settings.value("loc/tool", None)

        # loctoollabel = QLabel("location tool")
        self.locToolComboBox = QComboBox()
        # loctools = LOCTOOLS.keys()
        # self.locToolComboBox.addItems(loctools)

        # toolind = findComboBoxIndex(self.locToolComboBox, curtool)

        # self.locToolComboBox.setCurrentIndex(toolind)

        curroot = settings.value("{0}/rootPath".format(curtool), None)
        curbin = settings.value("{0}/binPath".format(curtool), None)

        self.rootlabel = QLabel("root directory")
        self.binlabel = QLabel("bin directory")

        self.rootedit = QLineEdit('')
        self.binedit = QLineEdit('')

        if curroot is not None:
            self.rootedit.setText(curroot)
        if curbin is not None:
            self.binedit.setText(curbin)

        rootBrowse = QPushButton('...', self)
        rootBrowse.clicked.connect(lambda: self.selectDirectory(self.rootedit))

        binBrowse = QPushButton('...', self)
        binBrowse.clicked.connect(lambda: self.selectDirectory(self.binedit))

        # self.locToolComboBox.currentIndexChanged.connect(self.updateUi)

        self.updateUi()

        layout = QGridLayout()
        # layout.addWidget(loctoollabel, 0, 0)
        # layout.addWidget(self.locToolComboBox, 0, 1)
        layout.addWidget(self.rootlabel, 1, 0)
        layout.addWidget(self.rootedit, 1, 1)
        layout.addWidget(rootBrowse, 1, 2)
        layout.addWidget(self.binlabel, 2, 0)
        layout.addWidget(self.binedit, 2, 1)
        layout.addWidget(binBrowse, 2, 2)

        self.setLayout(layout)

    def updateUi(self):
        curtool = self.locToolComboBox.currentText()
        # if curtool is not None:
        self.rootlabel.setText("{0} root directory".format(curtool))
        self.binlabel.setText("{0} bin directory".format(curtool))

    @staticmethod
    def selectDirectory(edit):
        selected_directory = QFileDialog.getExistingDirectory()
        # check if string is empty
        if selected_directory:
            edit.setText(selected_directory)

    def getValues(self):
        loctool = self.locToolComboBox.currentText()
        values = {"{0}/rootPath".format(loctool): self.rootedit.text(),
                  "{0}/binPath".format(loctool): self.binedit.text()}
        # "loc/tool": loctool}
        return values

    def resetValues(self, infile):
        para = PylotParameter(infile)
        nllocroot = para.get('nllocroot')
        nllocbin = para.get('nllocbin')
        loctool = self.locToolComboBox.setCurrentIndex(3)
        values = {"nll/rootPath": self.rootedit.setText("%s" % nllocroot),
                  "nll/binPath": self.binedit.setText("%s" % nllocbin)}


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
        return {'origintime': self.eventTimeEdit.dateTime().toPython(),
                'latitude': self.latEdit.text(),
                'longitude': self.lonEdit.text(),
                'depth': self.depEdit.text()}

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
        super(FilterOptionsDialog, self).__init__(parent)
        if parent is not None and parent.getFilters():
            self.filterOptions = parent.getFilters()
        # elif filterOptions is not None:
        #     self.filterOptions = filterOptions
        else:
            self.filterOptions = {'P': FilterOptions(),
                                  'S': FilterOptions()}

        self.setWindowTitle(titleString)
        self.filterOptionWidgets = {
            'P': FilterOptionsWidget(self.filterOptions['P'], self.parent().getAutoFilteroptions('P')),
            'S': FilterOptionsWidget(self.filterOptions['S'], self.parent().getAutoFilteroptions('S'))}
        self.setupUi()
        self.updateUi()
        self.connectButtons()

    def setupUi(self):
        self.main_layout = QtWidgets.QVBoxLayout()
        self.filter_layout = QtWidgets.QHBoxLayout()
        self.groupBoxes = {'P': QtWidgets.QGroupBox('P Filter'),
                           'S': QtWidgets.QGroupBox('S Filter')}

        settings = QSettings()
        overwriteFilter = get_bool(settings.value('useGuiFilter'))

        self.overwriteFilterCheckbox = QCheckBox('Overwrite filteroptions')
        self.overwriteFilterCheckbox.setToolTip('Overwrite filter settings for refined pick with GUI settings')
        self.overwriteFilterCheckbox.setChecked(overwriteFilter)
        self.overwriteFilterCheckbox.clicked.connect(self.toggleFilterOverwrite)

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok |
                                          QDialogButtonBox.Cancel)

        for phase in ['P', 'S']:
            groupbox = self.groupBoxes[phase]
            box_layout = QtWidgets.QVBoxLayout()
            groupbox.setLayout(box_layout)

            self.filter_layout.addWidget(groupbox)
            box_layout.addWidget(self.filterOptionWidgets[phase])

        self.main_layout.addLayout(self.filter_layout)
        self.main_layout.addWidget(self.overwriteFilterCheckbox)
        self.main_layout.addWidget(self.buttonBox)
        self.setLayout(self.main_layout)

    def toggleFilterOverwrite(self):
        if self.overwriteFilterCheckbox.isChecked():
            qmb = QMessageBox(self, icon=QMessageBox.Warning,
                              text='Warning: Overwriting automatic filter settings'
                                   ' for final pick will contaminate comparability'
                                   ' of automatic and manual picks! Continue?')
            qmb.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            qmb.setDefaultButton(QMessageBox.Yes)
            ret = qmb.exec_()
            if not ret == qmb.Yes:
                self.overwriteFilterCheckbox.setChecked(False)

        settings = QSettings()
        settings.setValue('useGuiFilter', self.overwriteFilterCheckbox.isChecked())

    def connectButtons(self):
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

    def accept(self):
        if not self.checkMinMax():
            QMessageBox.warning(self, "Value error",
                                "Maximum frequency must be at least the "
                                "same value as minimum frequency (notch)! "
                                "Adjusted maximum frequency automatically!")
            return
        self.updateUi()
        QDialog.accept(self)

    def checkMinMax(self):
        returnvals = []
        for foWidget in self.filterOptionWidgets.values():
            if foWidget.filterOptions._filtertype in ['highpass', 'lowpass']:
                returnvals.append(True)
                continue
            returnvals.append(foWidget.checkMin())
            returnvals.append(foWidget.checkMax())
        if all(returnvals):
            return True
        else:
            return False

    def updateUi(self):
        returnvals = []
        for foWidget in self.filterOptionWidgets.values():
            foWidget.updateUi()

    def getFilterOptions(self):
        filteroptions = {'P': self.filterOptionWidgets['P'].getFilterOptions(),
                         'S': self.filterOptionWidgets['S'].getFilterOptions()}
        return filteroptions

    def isComparable(self):
        return {phase: fltOptWid.comparable for (phase, fltOptWid) in self.filterOptionWidgets.items()}


class FilterOptionsWidget(QWidget):
    def __init__(self, filterOptions, filterOptionsAuto):
        super(FilterOptionsWidget, self).__init__()
        self.filterOptions = filterOptions
        self.filterOptionsAuto = filterOptionsAuto

        self.comparable = None

        _enable = True
        if self.getFilterOptions().getFilterType() is None:
            _enable = False

        self.freqminLabel = QLabel()
        self.freqminLabel.setText("minimum:")
        self.freqminSpinBox = QDoubleSpinBox()
        self.freqminSpinBox.setRange(5e-7, 1e6)
        self.freqminSpinBox.setDecimals(5)
        self.freqminSpinBox.setSingleStep(0.01)
        self.freqminSpinBox.setSuffix(' Hz')
        self.freqminSpinBox.setEnabled(_enable)

        self.freqmaxLabel = QLabel()
        self.freqmaxLabel.setText("maximum:")
        self.freqmaxSpinBox = QDoubleSpinBox()
        self.freqmaxSpinBox.setRange(5e-7, 1e6)
        self.freqmaxSpinBox.setDecimals(5)
        self.freqmaxSpinBox.setSingleStep(0.01)
        self.freqmaxSpinBox.setSuffix(' Hz')

        # if _enable:
        #     self.freqminSpinBox.setValue(self.getFilterOptions().getFreq()[0])
        #     if self.getFilterOptions().getFilterType() in ['bandpass',
        #                                                    'bandstop']:
        #         self.freqmaxSpinBox.setValue(
        #             self.getFilterOptions().getFreq()[1])
        # else:

        self.typeOptions = ["bandpass", "bandstop", "lowpass", "highpass"]

        self.resetButton = QPushButton('Reset')
        self.resetButton.setToolTip('Reset filter settings to settings for automatic picking.')

        self.orderLabel = QLabel()
        self.orderLabel.setText("Order:")
        self.orderSpinBox = QSpinBox()
        self.orderSpinBox.setRange(2, 10)
        self.orderSpinBox.setEnabled(_enable)
        self.selectTypeLabel = QLabel()
        self.selectTypeLabel.setText("Select filter type:")
        self.selectTypeCombo = QComboBox()
        self.selectTypeCombo.addItems(self.typeOptions)

        self.autoOrder = QLabel('{}')
        self.autoType = QLabel('{}')

        self.selectTypeLayout = QGridLayout()
        self.selectTypeLayout.addWidget(self.orderLabel, 0, 0)
        self.selectTypeLayout.addWidget(self.orderSpinBox, 1, 0)
        self.selectTypeLayout.addWidget(self.autoOrder, 1, 1)
        self.selectTypeLayout.addWidget(self.selectTypeLabel, 2, 0)
        self.selectTypeLayout.addWidget(self.selectTypeCombo, 3, 0)
        self.selectTypeLayout.addWidget(self.autoType, 3, 1)
        self.selectTypeLayout.setColumnStretch(0, 1)
        self.selectTypeLayout.setColumnStretch(1, 1)

        self.autoMinFreq = QLabel('{} Hz')
        self.autoMaxFreq = QLabel('{} Hz')

        self.manuLabel = QLabel('manual')
        self.autoLabel = QLabel('auto')

        self.freqGroupBox = QGroupBox("Frequency range")
        self.freqGroupLayout = QGridLayout()
        self.freqGroupLayout.addWidget(self.manuLabel, 0, 1)
        self.freqGroupLayout.addWidget(self.autoLabel, 0, 2)
        self.freqGroupLayout.addWidget(self.freqminLabel, 1, 0)
        self.freqGroupLayout.addWidget(self.freqminSpinBox, 1, 1)
        self.freqGroupLayout.addWidget(self.autoMinFreq, 1, 2)
        self.freqGroupLayout.addWidget(self.freqmaxLabel, 2, 0)
        self.freqGroupLayout.addWidget(self.freqmaxSpinBox, 2, 1)
        self.freqGroupLayout.addWidget(self.autoMaxFreq, 2, 2)
        self.freqGroupLayout.setColumnStretch(0, 1)
        self.freqGroupLayout.setColumnStretch(1, 1)
        self.freqGroupLayout.setColumnStretch(2, 1)
        self.freqGroupBox.setLayout(self.freqGroupLayout)

        self.freqmaxSpinBox.setEnabled(_enable)

        grid = QGridLayout()
        grid.addWidget(self.freqGroupBox, 0, 0)
        grid.addLayout(self.selectTypeLayout, 1, 0)
        grid.addWidget(self.resetButton, 2, 0)

        self.setLayout(grid)

        self.setMFtoWidget()

        self.resetButton.clicked.connect(self.setManuToAuto)
        self.orderSpinBox.valueChanged.connect(self.updateUi)
        self.selectTypeCombo.currentIndexChanged.connect(self.updateUi)
        self.orderSpinBox.valueChanged.connect(self.checkAutoManu)
        self.selectTypeCombo.currentIndexChanged.connect(self.checkAutoManu)
        self.freqminSpinBox.valueChanged.connect(self.checkAutoManu)
        self.freqmaxSpinBox.valueChanged.connect(self.checkAutoManu)
        self.checkAutoManu()

        self.setWindowModality(QtCore.Qt.WindowModality.ApplicationModal)

    def checkAutoManu(self):
        self.updateMFfromWidget()

        manuFilter = self.filterOptions
        autoFilter = self.filterOptionsAuto
        labels = {'order': self.autoOrder,
                  'type': self.autoType,
                  'freqmin': self.autoMinFreq,
                  'freqmax': self.autoMaxFreq}
        autoOptions = {'order': autoFilter.getOrder(),
                       'type': autoFilter.getFilterType(),
                       'freqmin': autoFilter.getFreq()[0],
                       'freqmax': autoFilter.getFreq()[1]}
        manuOptions = {'order': manuFilter.getOrder(),
                       'type': manuFilter.getFilterType(),
                       'freqmin': manuFilter.getFreq()[0],
                       'freqmax': manuFilter.getFreq()[1]}

        comparable = []

        for key, label in labels.items():
            text = label.text().format(autoOptions[key])
            label.setText(text)
            isEqual = autoOptions[key] == manuOptions[key]
            if isEqual:
                label.setStyleSheet('color: green')
            else:
                label.setStyleSheet('color: red')
            comparable.append(isEqual)

        self.comparable = all(comparable)

    def setMFtoWidget(self):
        try:
            freqmin, freqmax = self.filterOptions.getFreq()
        except TypeError as e:
            print(e)
            freqmin, freqmax = (0.1, 1.0)
        try:
            order = self.getFilterOptions().getOrder()
        except TypeError as e:
            print(e)
            order = 4
        ftype = self.getFilterOptions().getFilterType()

        self.selectTypeCombo.setCurrentIndex(self.typeOptions.index(ftype))
        self.freqminSpinBox.setValue(freqmin)
        self.freqmaxSpinBox.setValue(freqmax)
        self.orderSpinBox.setValue(order)

    def updateMFfromWidget(self):
        type = self.selectTypeCombo.currentText()
        freq = [self.freqminSpinBox.value(), self.freqmaxSpinBox.value()]
        self.filterOptions.setFilterType(type)
        self.filterOptions.setFreq(freq)
        self.filterOptions.setOrder(self.orderSpinBox.value())

    def setManuToAuto(self):
        self.filterOptions.setOrder(self.filterOptionsAuto.getOrder())
        self.filterOptions.setFilterType(self.filterOptionsAuto.getFilterType())
        self.filterOptions.setFreq(self.filterOptionsAuto.getFreq())
        self.setMFtoWidget()
        self.updateUi()

    def checkMin(self):
        if not self.freqminSpinBox.value() <= self.freqmaxSpinBox.value():
            self.freqmaxSpinBox.setValue(self.freqminSpinBox.value())
            return False
        return True

    def checkMax(self):
        if not self.freqminSpinBox.value() <= self.freqmaxSpinBox.value():
            self.freqminSpinBox.setValue(self.freqmaxSpinBox.value())
            return False
        return True

    def updateUi(self):
        type = self.selectTypeCombo.currentText()
        _enable = type in ['bandpass', 'bandstop']
        freq = [self.freqminSpinBox.value(), self.freqmaxSpinBox.value()]
        self.freqmaxLabel.setEnabled(True)
        self.freqmaxSpinBox.setEnabled(True)
        self.freqminLabel.setEnabled(True)
        self.freqminSpinBox.setEnabled(True)
        self.freqminLabel.setText("minimum:")
        self.freqmaxLabel.setText("maximum:")

        if not _enable:
            if type == 'highpass':
                self.freqminLabel.setText("cutoff:")
                self.freqmaxLabel.setEnabled(False)
                self.freqmaxSpinBox.setEnabled(False)
            elif type == 'lowpass':
                self.freqmaxLabel.setText("cutoff:")
                self.freqminLabel.setEnabled(False)
                self.freqminSpinBox.setEnabled(False)
        else:
            if not isSorted(freq):
                QMessageBox.warning(self, "Value error",
                                    "Maximum frequency must be at least the "
                                    "same value as minimum frequency (notch)! "
                                    "Adjusted maximum frequency automatically!")
                freq[1] = freq[0]
                self.freqmaxSpinBox.setValue(freq[1])
                self.freqmaxSpinBox.selectAll()
                self.freqmaxSpinBox.setFocus()

        self.updateMFfromWidget()

    def getFilterOptions(self):
        return self.filterOptions

    @staticmethod
    def getFilterObject():
        dlg = FilterOptionsDialog()
        if dlg.exec_():
            return dlg.getFilterOptions()
        return None


class LoadDataDlg(QDialog):
    def __init__(self, parent=None):
        super(LoadDataDlg, self).__init__(parent)

        pass


class HelpForm(QDialog):
    def __init__(self, parent=None,
                 page=QUrl('https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/')):
        super(HelpForm, self).__init__(parent)
        self.setAttribute(Qt.WA_DeleteOnClose)
        self.setAttribute(Qt.WA_GroupLeader)

        self.home_page = page

        back_icon = QIcon()
        back_icon.addPixmap(QPixmap(':/icons/back.png'))

        home_icon = QIcon()
        home_icon.addPixmap(QPixmap(':/icons/home.png'))

        backAction = QAction(back_icon, "&Back", self)
        backAction.setShortcut(QKeySequence.Back)
        homeAction = QAction(home_icon, "&Home", self)
        homeAction.setShortcut("Home")
        self.pageLabel = QLabel()

        toolBar = QToolBar()
        toolBar.addAction(backAction)
        toolBar.addAction(homeAction)
        toolBar.addWidget(self.pageLabel)
        self.webBrowser = QWebView()
        self.webBrowser.load(page)
        # self.webBrowser.load('C:/Shared/code/git/pylot/pylot/core/util/map_test.html')

        layout = QVBoxLayout()
        layout.addWidget(toolBar)
        layout.addWidget(self.webBrowser)
        self.setLayout(layout)

        backAction.triggered.connect(self.webBrowser.back)
        homeAction.triggered.connect(self.home)
        self.webBrowser.urlChanged.connect(self.updatePageTitle)

        self.resize(1280, 720)
        self.setWindowTitle("{0} Help".format(QApplication.applicationName()))

    def home(self):
        self.webBrowser.load(self.home_page)

    def updatePageTitle(self):
        self.pageLabel.setText(self.webBrowser.title())


class PickQualitiesFromXml(QWidget):
    """
            PyLoT widget PickQualitiesFromXml is a QWidget object. It is an UI that enables the user
            to create a plot showing the pick qualities in the event selected inside the QComboBox created
            by this Widget. The user can also choose to select all Events.
            The plot is being shown by a FigureCanvas from matlplotlib.
    """

    def __init__(self, parent=None, figure=Figure(), path="", inputVar=None):
        super(PickQualitiesFromXml, self).__init__(parent)

        self.fig = figure
        self.chooseBox = QComboBox()
        self.path = path
        self.inputs = inputVar
        self.setupUi()

    def setupUi(self):
        self.setWindowTitle("Get pick qualities from xml files")
        self.main_layout = QtWidgets.QVBoxLayout()
        self.figureC = FigureCanvas(self.fig)
        self.chooseBox = self.createComboBox()
        if self.chooseBox:
            self.main_layout.addWidget(self.chooseBox)
        self.main_layout.addWidget(self.figureC)
        self.setLayout(self.main_layout)

    def showUI(self):
        self.show()

    # Creates a QComboBox and adds all events in the current folder as options. Also gives the option to choose all events
    def createComboBox(self):
        self.chooseBox.addItems(glob.glob(os.path.join(os.path.dirname(self.path) + "/*/", '*.xml')))
        self.chooseBox.addItem("All")
        self.chooseBox.currentIndexChanged.connect(self.selectionChanged)
        self.chooseBox.setCurrentIndex(self.chooseBox.count() - 1)
        return self.chooseBox

    # Function that gets called when the user changes the current selection in the QComboBox. Redraws the plot with the new data
    def selectionChanged(self):
        self.figureC.setParent(None)
        if self.chooseBox.currentIndex() == self.chooseBox.count() - 1:
            (_, _, plot) = getQualitiesfromxml(self.path, self.inputs.get('timeerrorsP'),
                                               self.inputs.get('timeerrorsS'), plotflag=1)
            self.figureC = FigureCanvas(plot)
        else:
            (_, _, plot) = getQualitiesfromxml(self.path, self.inputs.get('timeerrorsP'),
                                               self.inputs.get('timeerrorsS'), plotflag=1,
                                               xmlnames=[self.chooseBox.currentText()])
            self.figureC = FigureCanvas(plot)
        self.figureC.draw()
        self.main_layout.addWidget(self.figureC)
        self.setLayout(self.main_layout)


class SourceSpecWindow(QWidget):
    def __init__(self, parent=None, figure=Figure()):
        super(SourceSpecWindow, self).__init__(parent)
        self.main_layout = QVBoxLayout()
        self.setWindowTitle("Display source spectrum from selected trace")
        self.fig = figure

    def setupUi(self):
        self.figureC = FigureCanvas(self.fig)
        self.main_layout.addWidget(self.figureC)
        self.setLayout(self.main_layout)


class ChooseWaveFormWindow(QWidget):
    def __init__(self, parent=None, WaveForms=[], traces=[], stream=None, chooseB=False):
        super(ChooseWaveFormWindow, self).__init__(parent)
        self.main_layout = QVBoxLayout()
        self.setWindowTitle("Choose trace to display source spectrum")
        self.wFs = WaveForms
        self.chooseBoxTraces = QComboBox()
        self.chooseBoxComponent = QComboBox()
        self.submitButton = QPushButton(text='test')

        self.chooseB = chooseB

        self.verticalButton = QPushButton(text='Z')
        self.northButton = QPushButton(text='N')
        self.eastButton = QPushButton(text='E')
        self.component = ''

        self.traces = traces
        self.stream = stream
        self.setupUI()
        self.currentSpectro = Figure()

    def setupUI(self):
        self.submitButton.clicked.connect(self.submit)
        self.createComboBoxTraces()
        self.main_layout.addWidget(self.chooseBoxTraces)

        if self.chooseB:
            self.createComboBoxComponent()
            self.main_layout.addWidget(self.submitButton)
            self.main_layout.addWidget(self.chooseBoxComponent)
        else:
            self.createButtonsComponent()

        self.setLayout(self.main_layout)

    def submit(self):
        matplotlib.pyplot.close(self.currentSpectro)
        t = self.chooseBoxTraces.currentText() + " " + self.chooseBoxComponent.currentText()
        #self.currentSpectro = self.traces[
        #    self.chooseBoxTraces.currentText()[3:]][self.chooseBoxComponent.currentText()].spectrogram(show=False, title=t)
        #self.currentSpectro.show()
        self.applyFFT()

    def applyFFT(self, trace):
        tra = self.traces[self.chooseBoxTraces.currentText()[3:]]['Z']
        transformed = abs(np.fft.rfft(tra.data))
        print ( transformed )
        matplotlib.pyplot.plot ( transformed )
        matplotlib.pyplot.show()

    def applyFFTs(self, tra):
        transformed = abs(np.fft.rfft(tra.data))
        print ( transformed )
        matplotlib.pyplot.plot ( transformed )
        matplotlib.pyplot.show()

    def submitN(self):
        matplotlib.pyplot.close(self.currentSpectro)
        t = self.chooseBoxTraces.currentText() + " " + self.chooseBoxComponent.currentText()
        self.currentSpectro = self.traces[
            self.chooseBoxTraces.currentText()[3:]]['N'].spectrogram(show=False, title=t)
        self.currentSpectro.show()

    def submitE(self):
        matplotlib.pyplot.close(self.currentSpectro)
        t = self.chooseBoxTraces.currentText() + " " + self.chooseBoxComponent.currentText()
        self.currentSpectro = self.traces[
            self.chooseBoxTraces.currentText()[3:]]['E'].spectrogram(show=False, title=t)
        self.currentSpectro.show()

    # Creates a QComboBox and adds all traces provided
    def createComboBoxTraces(self):
        if len(self.wFs) <= 0:
            raise 'No traces provided'
        self.chooseBoxTraces.addItems(self.wFs)
        self.chooseBoxTraces.currentIndexChanged.connect(self.selectionChanged)
        self.chooseBoxTraces.setCurrentIndex(0)
        return self.chooseBoxTraces

    def createButtonsComponent(self):
        self.northButton.clicked.connect(self.submitN)
        self.eastButton.clicked.connect(self.submitE)
        self.verticalButton.clicked.connect(self.submitZ)

        self.main_layout.addWidget(self.verticalButton)
        self.main_layout.addWidget(self.northButton)
        self.main_layout.addWidget(self.eastButton)

    def createComboBoxComponent(self):
        self.chooseBoxComponent.addItems(['Z', 'N', 'E'])

    # Function that gets called when the user changes the current selection in the QComboBox.
    def selectionChanged(self):
        pass


class SpectrogramTab(QWidget):
    def __init__(self, traces, wfdata, parent=None):
        super(SpectrogramTab, self).__init__(parent)
        self.setupUi()
        self.traces = traces
        self.wfdata = wfdata

    def setupUi(self):
        pass
    def makeSpecFig(self, direction = 'Z', height = 0, width = 0, parent = None):

        i = 0
        grams = []
        figure, axis = matplotlib.pyplot.subplots(len(self.traces), sharex=True)

        start, end = full_range(self.wfdata)

        if height != 0 and width != 0:
            figure.figsize = (width, height)
            figure.set_figwidth = width
            figure.set_figheight = height

        #figure.tight_layout()

        for t in self.traces:
            tra = self.traces[t][direction]
            #print(start, end)

            # Set Title
            if i == 0:
                if direction == 'Z':
                    figure.suptitle("section: vertical components")
                elif direction == 'E':
                    figure.suptitle("section: east-west components")
                elif direction == 'N':
                    figure.suptitle("section: north-south components")
                axis[i].vlines(0, axis[i].get_ylim()[0], axis[i].get_ylim()[1],
                          colors='m', linestyles='dashed',
                          linewidth=2)

            # Different axis settings for visual improvements
            # axis[i].set_xlim(left=0, right=end - start)
            # axis[i].spines['top'].set_visible(False)
            # axis[i].spines['right'].set_visible(False)
            # # axis[i].spines['left'].set_visible(False)

            # axis[i].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            # if not (len(self.traces) == i - 1):
                # axis[i].spines['bottom'].set_visible(False)
            # axis[i].set_yticks([])
            # axis[i].set_ylabel(t, loc='center', rotation='horizontal')

            #ax.axhline(n, color="0.5", lw=0.5)


            grams.append(tra.spectrogram(show=False, axes=axis[i]))
            i+=1

        #figure.setXLims([0, end - start])
        figure.set_tight_layout(True)
        fC = FigureCanvas(figure)
        return fC

    #for t in self.traces:
       # tra = self.traces[t]['Z']
       # transformed = abs(np.fft.rfft(tra.data))
       # axis[i].plot(transformed, label=t)
       # # axis[i].tick_params(labelbottom=False)
       # axis[i].spines['top'].set_visible(False)
       # axis[i].spines['right'].set_visible(False)
       # axis[i].spines['left'].set_visible(False)
       # if not (len(self.traces) == i - 1):
       #     axis[i].spines['bottom'].set_visible(False)
       # axis[i].set_yticks([])
       # axis[i].set_ylabel(t, loc='center', rotation='horizontal')
       # # axis[i].axis('off')
       # i += 1
       # # self.applyFFTs(t)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
