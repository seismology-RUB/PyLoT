# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:27:35 2014

@author: sebastianw
"""

import os
import sys
import getpass
import warnings
import copy
import datetime
import numpy as np
import matplotlib.pyplot as plt

try:
    import pyqtgraph as pg
except:
    pg = None

from matplotlib.figure import Figure
from pylot.core.util.utils import find_horizontals

try:
    from matplotlib.backends.backend_qt4agg import FigureCanvas
except ImportError:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.widgets import MultiCursor
from PySide import QtCore, QtGui
from PySide.QtGui import QAction, QApplication, QCheckBox, QComboBox, \
    QDateTimeEdit, QDialog, QDialogButtonBox, QDoubleSpinBox, QGroupBox, \
    QGridLayout, QIcon, QKeySequence, QLabel, QLineEdit, QMessageBox, \
    QPixmap, QSpinBox, QTabWidget, QToolBar, QVBoxLayout, QWidget, \
    QPushButton, QFileDialog, QInputDialog, QKeySequence
from PySide.QtCore import QSettings, Qt, QUrl, Signal, Slot
from PySide.QtWebKit import QWebView
from obspy import Stream, UTCDateTime
from obspy.core.util import AttribDict
from obspy.taup import TauPyModel
from obspy.taup.utils import get_phase_names
from pylot.core.io.data import Data
from pylot.core.io.inputs import FilterOptions, PylotParameter
from pylot.core.pick.utils import getSNR, earllatepicker, getnoisewin, \
    getResolutionWindow
from pylot.core.pick.compare import Comparison
from pylot.core.util.defaults import OUTPUTFORMATS, FILTERDEFAULTS, ALTSUFFIX,\
    LOCTOOLS, SetChannelComponents
from pylot.core.util.utils import prepTimeAxis, full_range, scaleWFData, \
    demeanTrace, isSorted, findComboBoxIndex, clims
from autoPyLoT import autoPyLoT
from pylot.core.util.thread import Thread

if sys.version_info.major == 3:
    import icons_rc_3 as icons_rc
elif sys.version_info.major == 2:
    import icons_rc_2 as icons_rc
else:
    raise ImportError ('Could not determine python version.')

if pg:
    pg.setConfigOption ('background', 'w')
    pg.setConfigOption ('foreground', 'k')
    pg.setConfigOptions (antialias=True)
    # pg.setConfigOption('leftButtonPan', False)


def getDataType(parent):
    type = QInputDialog ().getItem (parent, "Select phases type", "Type:",
                                    ["manual", "automatic"])

    if type[0].startswith ('auto'):
        type = 'auto'
    else:
        type = type[0]

    return type


def plot_pdf(_axes, x, y, annotation, bbox_props, xlabel=None, ylabel=None,
             title=None):
    # try method or data
    try:
        _axes.plot (x, y ())  # y provided as method
    except:
        _axes.plot (x, y)  # y provided as data

    if title:
        _axes.set_title (title)
    if xlabel:
        _axes.set_xlabel (xlabel)
    if ylabel:
        _axes.set_ylabel (ylabel)
    _anno = _axes.annotate (annotation, xy=(.05, .5), xycoords='axes fraction')
    _anno.set_bbox (bbox_props)

    return _axes


def createAction(parent, text, slot=None, shortcut=None, icon=None,
                 tip=None, checkable=False):
    """
    :rtype : ~PySide.QtGui.QAction
    """
    action = QAction (text, parent)
    if icon is not None:
        action.setIcon (icon)
    if shortcut is not None:
        action.setShortcut (shortcut)
    if tip is not None:
        action.setToolTip (tip)
    if slot is not None:
        action.triggered.connect (slot)
    if checkable:
        action.setCheckable (True)
    return action

def loopIdentifyPhase(phase):
    phase_copy = phase
    while not identifyPhase (phase_copy):
        identified = False
        for alt_suf in ALTSUFFIX:
            if phase_copy.endswith (alt_suf):
                phase_copy = phase_copy.split (alt_suf)[0]
                identified = True
        if not identified:
            phase_copy = phase_copy[:-1]
        if len (phase_copy) < 1:
            print ('Warning: Could not identify phase {}!'.format (phase))
            return
    return phase_copy

def identifyPhase(phase):
    # common phase suffix for P and S
    common_P = ['P', 'p']
    common_S = ['S', 's']
    if phase[-1] in common_P:
        return 'P'
    if phase[-1] in common_S:
        return 'S'
    else:
        return False


class ComparisonDialog (QDialog):
    def __init__(self, c, parent=None):
        self._data = c
        self._stats = c.stations
        self._canvas = PlotWidget (self)
        self._widgets = dict (stationsComboBox=None,
                              phasesComboBox=None,
                              histCheckBox=None)
        self._phases = 'PS'
        self._plotprops = dict (station=self.stations[0], phase=self.phases[0])
        super (ComparisonDialog, self).__init__ (parent)
        self.setupUI ()
        self.plotcomparison ()

    def setupUI(self):

        _outerlayout = QVBoxLayout (self)
        _innerlayout = QVBoxLayout ()

        _stats_combobox = QComboBox (self)
        _stats_combobox.setObjectName ('stationsComboBox')
        _stats_combobox.setEditable (True)
        _stats_combobox.setInsertPolicy (QComboBox.NoInsert)
        _stats_combobox.addItems (sorted (self.stations))
        _stats_combobox.editTextChanged.connect (self.prepareplot)
        self.widgets = _stats_combobox

        _phases_combobox = QComboBox (self)
        _phases_combobox.setObjectName ('phasesComboBox')
        _phases_combobox.addItems (['P', 'S'])
        _phases_combobox.currentIndexChanged.connect (self.prepareplot)
        self.widgets = _phases_combobox

        _hist_checkbox = QCheckBox ('Show histograms', self)
        _hist_checkbox.setObjectName ('histCheckBox')
        _hist_checkbox.stateChanged.connect (self.plothist)
        self.widgets = _hist_checkbox

        _toolbar = QToolBar (self)
        _toolbar.addWidget (_stats_combobox)
        _toolbar.addWidget (_phases_combobox)
        _toolbar.addWidget (_hist_checkbox)

        _buttonbox = QDialogButtonBox (QDialogButtonBox.Close)

        _innerlayout.addWidget (self.canvas)
        _innerlayout.addWidget (_buttonbox)

        _outerlayout.addWidget (_toolbar)
        _outerlayout.addLayout (_innerlayout)

        _buttonbox.rejected.connect (self.reject)

        # finally layout the entire dialog
        self.setLayout (_outerlayout)

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
            if key not in self.plotprops.keys ():
                raise KeyError ("'key' {0} not found in "
                                "ComparisonDialog.plotprops keys.".format (key))
        except ValueError:
            raise ValueError ("Pass an iterable with two items")
        else:
            self._plotprops[key] = value

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        assert not isinstance (data, Comparison)
        self.stations = data.stations
        self._data = data

    @property
    def widgets(self):
        return self._widgets

    @widgets.setter
    def widgets(self, widget):
        name = widget.objectName ()
        if name in self.widgets.keys ():
            self._widgets[name] = widget

    def clf(self):
        self.canvas.figure.clf ()

    def hasvalue(self, sender):
        text = sender.currentText ()
        index = sender.findText (text.upper ())
        return index

    def prepareplot(self):
        try:
            _widget = self.sender ()
            name = _widget.objectName ()
            text = _widget.currentText ().upper ()
            index = self.hasvalue (_widget)
            if name == 'stationsComboBox' and index is not -1:
                _widget.setCurrentIndex (index)
                self.plotprops = ('station', text)
            elif name == 'phasesComboBox':
                self.plotprops = ('phase', text)
        except ValueError:
            raise ValueError ('No sender widget given!')
        finally:
            self.plotcomparison ()

    def plotcomparison(self):
        from matplotlib import gridspec

        _gs = gridspec.GridSpec (3, 2)
        self.clf ()
        _axes = self.canvas.figure.add_subplot (_gs[0:2, :])
        _ax1 = self.canvas.figure.add_subplot (_gs[2, 0])
        _ax2 = self.canvas.figure.add_subplot (_gs[2, 1])

        # _axes.cla()
        station = self.plotprops['station']
        phase = self.plotprops['phase']
        pdf = self.data.comparison[station][phase]
        x, y, std, exp = pdf.axis, pdf.data, pdf.standard_deviation (), \
                         pdf.expectation ()

        annotation = "%s difference on %s\n" \
                     "expectation: %7.4f s\n" \
                     "std: %7.4f s" % (phase, station,
                                       exp, std)
        bbox_props = dict (boxstyle='round', facecolor='lightgrey', alpha=.7)

        plot_pdf (_axes, x, y, annotation, bbox_props, 'time difference [s]',
                  'propability density [-]', phase)

        pdf_a = copy.deepcopy (self.data.get ('auto')[station][phase])
        pdf_m = copy.deepcopy (self.data.get ('manu')[station][phase])

        xauto, yauto, stdauto, expauto, alim = pdf_a.axis, pdf_a.data (), \
                                               pdf_a.standard_deviation (), \
                                               pdf_a.expectation (), \
                                               pdf_a.limits ()
        xmanu, ymanu, stdmanu, expmanu, mlim = pdf_m.axis, pdf_m.data (), \
                                               pdf_m.standard_deviation (), \
                                               pdf_m.expectation (), \
                                               pdf_m.limits ()
        # find common limits
        lims = clims (alim, mlim)
        # relative x axis
        x0 = lims[0]
        xmanu -= x0
        xauto -= x0
        lims = [lim - x0 for lim in lims]
        x0 = UTCDateTime (x0)

        # set annotation text
        mannotation = "probability density of manual pick\n" \
                      "expectation: %7.4f s\n" \
                      "std: %7.4f s" % (expmanu - x0.timestamp, stdmanu)

        aannotation = "probability density of automatic pick\n" \
                      "expectation: %7.4f s\n" \
                      "std: %7.4f s" % (expauto - x0.timestamp, stdauto)

        _ax1 = plot_pdf (_ax1, xmanu, ymanu, mannotation,
                         bbox_props=bbox_props, xlabel='seconds since '
                                                       '{0}'.format (x0),
                         ylabel='probability density [-]')
        _ax1.set_xlim (lims)

        _ax2 = plot_pdf (_ax2, xauto, yauto, aannotation,
                         bbox_props=bbox_props, xlabel='seconds since '
                                                       '{0}'.format (x0))
        _ax2.set_xlim (lims)

        _gs.update (wspace=0.5, hspace=0.5)

        self.canvas.draw ()

    def plothist(self):
        name = self.sender ().objectName ()
        if self.widgets[name].isChecked ():
            for wname, widget in self.widgets.items ():
                if wname != name:
                    self.widgets[wname].setEnabled (False)
            self.canvas.figure.clf ()
            _axPstd, _axPexp = self.canvas.figure.add_subplot (221), self.canvas.figure.add_subplot (223)
            _axSstd, _axSexp = self.canvas.figure.add_subplot (222), self.canvas.figure.add_subplot (224)
            axes_dict = dict (P=dict (std=_axPstd, exp=_axPexp),
                              S=dict (std=_axSstd, exp=_axSexp))
            bbox_props = dict (boxstyle='round', facecolor='lightgrey', alpha=.7)
            for phase in self.phases:
                std = self.data.get_std_array (phase)
                std = std[np.isfinite (std)]
                stdxlims = [0., 1.2 * max (std)]
                exp = self.data.get_expectation_array (phase)
                exp = exp[np.isfinite (exp)]
                eps_exp = 0.05 * (max (exp) - min (exp))
                expxlims = [min (exp) - eps_exp, max (exp) + eps_exp]
                axes_dict[phase]['std'].hist (std, range=stdxlims, bins=20, normed=False)
                axes_dict[phase]['exp'].hist (exp, range=expxlims, bins=20,
                                              normed=False)
                std_annotation = "Distribution curve for {phase} differences'\n" \
                                 "standard deviations (all stations)\n" \
                                 "number of samples: {nsamples}".format (phase=phase, nsamples=len (std))
                _anno_std = axes_dict[phase]['std'].annotate (std_annotation, xy=(.05, .8), xycoords='axes fraction')
                _anno_std.set_bbox (bbox_props)
                exp_annotation = "Distribution curve for {phase} differences'\n" \
                                 "expectations (all stations)\n" \
                                 "number of samples: {nsamples}".format (phase=phase, nsamples=len (exp))
                _anno_exp = axes_dict[phase]['exp'].annotate (exp_annotation, xy=(.05, .8), xycoords='axes fraction')
                _anno_exp.set_bbox (bbox_props)
                axes_dict[phase]['exp'].set_xlabel ('expectation [s]')
                axes_dict[phase]['std'].set_xlabel ('standard deviation [s]')

            for ax in axes_dict['P'].values ():
                ax.set_ylabel ('frequency [-]')

            self.canvas.draw ()
        else:
            for wname, widget in self.widgets.items ():
                if wname != name:
                    self.widgets[wname].setEnabled (True)
            self.canvas.figure.clf ()
            self.plotcomparison ()


class PlotWidget (FigureCanvas):
    def __init__(self, parent=None, xlabel='x', ylabel='y', title='Title'):
        self._parent = parent
        self._fig = Figure ()
        self._xl = xlabel
        self._yl = ylabel
        self._title = title
        super (PlotWidget, self).__init__ (self.figure)

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


class WaveformWidgetPG (QtGui.QWidget):
    def __init__(self, parent=None, xlabel='x', ylabel='y', title='Title'):
        QtGui.QWidget.__init__ (self, parent)  # , 1)
        self.setParent (parent)
        self._parent = parent
        # attribute plotdict is a dictionary connecting position and a name
        self.plotdict = dict ()
        # create plot
        self.main_layout = QtGui.QVBoxLayout ()
        self.label = QtGui.QLabel ()
        self.setLayout (self.main_layout)
        self.plotWidget = pg.PlotWidget (title=title, autoDownsample=True)
        self.main_layout.addWidget (self.plotWidget)
        self.main_layout.addWidget (self.label)
        self.plotWidget.showGrid (x=False, y=True, alpha=0.2)
        self.plotWidget.hideAxis ('bottom')
        self.plotWidget.hideAxis ('left')
        self.reinitMoveProxy ()
        self._proxy = pg.SignalProxy (self.plotWidget.scene ().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)

    def reinitMoveProxy(self):
        self.vLine = pg.InfiniteLine (angle=90, movable=False)
        self.hLine = pg.InfiniteLine (angle=0, movable=False)
        self.plotWidget.addItem (self.vLine, ignoreBounds=True)
        self.plotWidget.addItem (self.hLine, ignoreBounds=True)

    def mouseMoved(self, evt):
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if self.plotWidget.sceneBoundingRect ().contains (pos):
            mousePoint = self.plotWidget.getPlotItem ().vb.mapSceneToView (pos)
            x, y, = (mousePoint.x (), mousePoint.y ())
            # if x > 0:# and index < len(data1):
            wfID = self._parent.getWFID (y)
            station = self._parent.getStationName (wfID)
            if self._parent.get_current_event ():
                self.label.setText ("station = {}, t = {} [s]".format (station, x))
            self.vLine.setPos (mousePoint.x ())
            self.hLine.setPos (mousePoint.y ())

    def getPlotDict(self):
        return self.plotdict

    def setPlotDict(self, key, value):
        self.plotdict[key] = value

    def clearPlotDict(self):
        self.plotdict = dict ()

    def getParent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    def plotWFData(self, wfdata, title=None, zoomx=None, zoomy=None,
                   noiselevel=None, scaleddata=False, mapping=True,
                   component='*', nth_sample=1, iniPick=None, verbosity=0):
        self.title = title
        self.clearPlotDict ()
        wfstart, wfend = full_range (wfdata)
        nmax = 0

        settings = QSettings ()
        compclass = settings.value ('compclass')
        if not compclass:
            print ('Warning: No settings for channel components found. Using default')
            compclass = SetChannelComponents ()

        if not component == '*':
            alter_comp = compclass.getCompPosition (component)
            # alter_comp = str(alter_comp[0])

            st_select = wfdata.select (component=component)
            st_select += wfdata.select (component=alter_comp)
        else:
            st_select = wfdata

        # list containing tuples of network, station, channel (for sorting)
        nsc = []
        for trace in st_select:
            nsc.append ((trace.stats.network, trace.stats.station, trace.stats.channel))
        nsc.sort ()
        nsc.reverse ()
        plots = []

        try:
            self.plotWidget.getPlotItem ().vb.setLimits (xMin=float (0),
                                                         xMax=float (wfend - wfstart),
                                                         yMin=-0.5,
                                                         yMax=len (nsc) + 0.5)
        except:
            print ('Warning: Could not set zoom limits')

        for n, (network, station, channel) in enumerate (nsc):
            st = st_select.select (network=network, station=station, channel=channel)
            trace = st[0]
            if mapping:
                comp = channel[-1]
                n = compclass.getPlotPosition (str (comp))
                # n = n[0]
            if n > nmax:
                nmax = n
            if verbosity:
                msg = 'plotting %s channel of station %s' % (channel, station)
                print (msg)
            stime = trace.stats.starttime - wfstart
            time_ax = prepTimeAxis (stime, trace)
            if time_ax is not None:
                if not scaleddata:
                    trace.detrend ('constant')
                    trace.normalize (np.max (np.abs (trace.data)) * 2)
                times = [time for index, time in enumerate (time_ax) if not index % nth_sample]
                data = [datum + n for index, datum in enumerate (trace.data) if not index % nth_sample]
                plots.append ((times, data))
                self.setPlotDict (n, (station, channel, network))
        self.xlabel = 'seconds since {0}'.format (wfstart)
        self.ylabel = ''
        self.setXLims ([0, wfend - wfstart])
        self.setYLims ([-0.5, nmax + 0.5])
        return plots

    # def getAxes(self):
    #     return self.axes

    # def getXLims(self):
    #     return self.getAxes().get_xlim()

    # def getYLims(self):
    #     return self.getAxes().get_ylim()

    def setXLims(self, lims):
        vb = self.plotWidget.getPlotItem ().getViewBox ()
        vb.setXRange (float (lims[0]), float (lims[1]), padding=0)

    def setYLims(self, lims):
        vb = self.plotWidget.getPlotItem ().getViewBox ()
        vb.setYRange (float (lims[0]), float (lims[1]), padding=0)

    def setYTickLabels(self, pos, labels):
        ticks = zip (pos, labels)
        minorTicks = [(0, 0) for item in labels]
        # leftAx.tickLength = 5
        # leftAx.orientation = 'right'
        self.getAxItem ('left').setTicks ([ticks, minorTicks])

    def updateXLabel(self, text):
        self.getAxItem ('bottom').setLabel (text)
        self.draw ()

    def updateYLabel(self, text):
        self.getAxItem ('left').setLabel (text)
        self.draw ()

    def getAxItem(self, position):
        return self.plotWidget.getPlotItem ().axes[position]['item']

    def updateTitle(self, text):
        self.plotWidget.getPlotItem ().setTitle (text)
        self.draw ()

    def updateWidget(self):  # , xlabel, ylabel, title):
        self.updateXLabel (self.xlabel)
        self.updateYLabel (self.ylabel)
        self.updateTitle (self.title)

    def draw(self):
        pass


class WaveformWidget (FigureCanvas):
    def __init__(self, parent=None, xlabel='x', ylabel='y', title='Title'):

        self._parent = parent
        self.figure = Figure ()
        self.figure.set_facecolor ((.92, .92, .92))
        # attribute plotdict is a dictionary connecting position and a name
        self.plotdict = dict ()
        # create axes
        self.axes = self.figure.add_subplot (111)
        # initialize super class
        super (WaveformWidget, self).__init__ (self.figure)
        # add an cursor for station selection
        self.multiCursor = MultiCursor (self.figure.canvas, (self.axes,),
                                        horizOn=True, useblit=True,
                                        color='m', lw=1)
        # update labels of the entire widget
        self.updateWidget (xlabel, ylabel, title)
        try:
            self.figure.tight_layout ()
        except:
            pass

    def getPlotDict(self):
        return self.plotdict

    def setPlotDict(self, key, value):
        self.plotdict[key] = value

    def clearPlotDict(self):
        self.plotdict = dict ()

    def getParent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    def plotWFData(self, wfdata, title=None, zoomx=None, zoomy=None,
                   noiselevel=None, scaleddata=False, mapping=True,
                   component='*', nth_sample=1, iniPick=None, verbosity=0):
        self.getAxes ().cla ()
        self.clearPlotDict ()
        wfstart, wfend = full_range (wfdata)
        nmax = 0

        settings = QSettings ()
        compclass = settings.value ('compclass')
        if not compclass:
            print ('Warning: No settings for channel components found. Using default')
            compclass = SetChannelComponents ()

        if not component == '*':
            alter_comp = compclass.getCompPosition (component)
            # alter_comp = str(alter_comp[0])

            st_select = wfdata.select (component=component)
            st_select += wfdata.select (component=alter_comp)
        else:
            st_select = wfdata

        # list containing tuples of network, station, channel (for sorting)
        nsc = []
        for trace in st_select:
            nsc.append ((trace.stats.network, trace.stats.station, trace.stats.channel))
        nsc.sort ()
        nsc.reverse ()

        for n, (network, station, channel) in enumerate (nsc):
            st = st_select.select (network=network, station=station, channel=channel)
            trace = st[0]
            if mapping:
                comp = channel[-1]
                n = compclass.getPlotPosition (str (comp))
                # n = n[0]
            if n > nmax:
                nmax = n
            if verbosity:
                msg = 'plotting %s channel of station %s' % (channel, station)
                print (msg)
            stime = trace.stats.starttime - wfstart
            time_ax = prepTimeAxis (stime, trace)
            if time_ax is not None:
                if not scaleddata:
                    trace.detrend ('constant')
                    trace.normalize (np.max (np.abs (trace.data)) * 2)
                times = [time for index, time in enumerate (time_ax) if not index % nth_sample]
                data = [datum + n for index, datum in enumerate (trace.data) if not index % nth_sample]
                self.getAxes ().plot (times, data, 'k', linewidth=0.7)
                if noiselevel is not None:
                    for level in noiselevel:
                        self.getAxes ().plot ([time_ax[0], time_ax[-1]],
                                              [level, level], '--k')
                if iniPick:
                    ax = self.getAxes ()
                    ax.vlines (iniPick, ax.get_ylim ()[0], ax.get_ylim ()[1],
                               colors='m', linestyles='dashed',
                               linewidth=2)
                self.setPlotDict (n, (station, channel, network))
        xlabel = 'seconds since {0}'.format (wfstart)
        ylabel = ''
        self.updateWidget (xlabel, ylabel, title)
        self.setXLims ([0, wfend - wfstart])
        self.setYLims ([-0.5, nmax + 0.5])
        if zoomx is not None:
            self.setXLims (zoomx)
        if zoomy is not None:
            self.setYLims (zoomy)
        self.draw ()

    def getAxes(self):
        return self.axes

    def getXLims(self):
        return self.getAxes ().get_xlim ()

    def getYLims(self):
        return self.getAxes ().get_ylim ()

    def setXLims(self, lims):
        self.getAxes ().set_xlim (lims)

    def setYLims(self, lims):
        self.getAxes ().set_ylim (lims)

    def setYTickLabels(self, pos, labels):
        self.getAxes ().set_yticks (list (pos))
        self.getAxes ().set_yticklabels (labels)
        self.draw ()

    def updateXLabel(self, text):
        self.getAxes ().set_xlabel (text)
        self.draw ()

    def updateYLabel(self, text):
        self.getAxes ().set_ylabel (text)
        self.draw ()

    def updateTitle(self, text):
        self.getAxes ().set_title (text, verticalalignment='bottom')
        self.draw ()

    def updateWidget(self, xlabel, ylabel, title):
        self.updateXLabel (xlabel)
        self.updateYLabel (ylabel)
        self.updateTitle (title)

    def insertLabel(self, pos, text):
        pos = pos / max (self.getAxes ().ylim)
        axann = self.getAxes ().annotate (text, xy=(.03, pos),
                                          xycoords='axes fraction')
        axann.set_bbox (dict (facecolor='lightgrey', alpha=.6))


class PhaseDefaults (QtGui.QDialog):
    def __init__(self, parent=None, nrow=10,
                 phase_defaults=['P', 'S'],
                 current_phases=[]):
        super (PhaseDefaults, self).__init__ (parent)
        self.nrow = nrow
        self.checktoggle = True
        self.main_layout = QtGui.QVBoxLayout ()
        self.sub_layout = QtGui.QGridLayout ()
        self.phase_names = phase_defaults
        self.current_phases = current_phases
        self.setButtons ()
        self.setupUi ()
        self.connectSignals ()
        self.setWindowTitle ('Default Phases')
        self.selected_phases = []

    def setButtons(self):
        self._check_all_button = QtGui.QPushButton('Check/Uncheck all')
        self._buttonbox = QDialogButtonBox (QDialogButtonBox.Ok |
                                            QDialogButtonBox.Cancel)

    def setupUi(self):
        self.setLayout (self.main_layout)
        self.main_layout.addWidget(self._check_all_button)
        self.main_layout.addLayout (self.sub_layout)
        self.init_phase_layout ()
        self.main_layout.addWidget (self._buttonbox)

    def connectSignals(self):
        self._check_all_button.clicked.connect(self.toggleAllChecked)
        self._buttonbox.accepted.connect (self.accept)
        self._buttonbox.rejected.connect (self.reject)

    def toggleAllChecked(self):
        for box in self._checkboxes.values():
            box.setChecked(self.checktoggle)
        self.checktoggle = not self.checktoggle

    def init_phase_layout(self):
        self._checkboxes = {}
        row = 0
        column = 0
        for index, phase in enumerate (self.phase_names):
            if row > self.nrow:
                column += 1
                row = 0
            checkbox = self.create_phase_box (phase)
            self.sub_layout.addWidget (checkbox,
                                       row, column)
            self._checkboxes[phase] = checkbox
            checkbox.setChecked (bool (phase in self.current_phases))
            row += 1

    def create_phase_box(self, phase_name):
        checkbox = QtGui.QCheckBox (phase_name)
        return checkbox

    def update_selected_phases(self):
        self.selected_phases = []
        for phase in self.phase_names:
            if self._checkboxes[phase].isChecked ():
                self.selected_phases.append (phase)

    def accept(self):
        self.update_selected_phases ()
        QtGui.QDialog.accept (self)


class PickDlg (QDialog):
    update_picks = QtCore.Signal (dict)

    def __init__(self, parent=None, data=None, station=None, network=None, picks=None,
                 autopicks=None, rotate=False, parameter=None, embedded=False, model='iasp91'):
        super (PickDlg, self).__init__ (parent)

        # initialize attributes
        self.parameter = parameter
        self._embedded = embedded
        self.station = station
        self.network = network
        self.rotate = rotate
        self.components = 'ZNE'
        self.currentPhase = None
        self.phaseText = []
        self.arrivals = []
        self.arrivalsText = []
        self.cidpick = []
        settings = QSettings ()
        pylot_user = getpass.getuser ()
        self._user = settings.value ('user/Login', pylot_user)
        self._dirty = False
        if picks:
            self.picks = picks
            self._init_picks = copy.deepcopy (picks)
        else:
            self.picks = {}
            self._init_picks = {}
        if autopicks:
            self.autopicks = autopicks
            self._init_autopicks = copy.deepcopy (autopicks)
        else:
            self.autopicks = {}
            self._init_autopicks = {}
        if hasattr (self.parent (), 'filteroptions'):
            self.filteroptions = self.parent ().filteroptions
        else:
            self.filteroptions = FILTERDEFAULTS
        self.pick_block = False
        self.nextStation = QtGui.QCheckBox ('Continue with next station.')

        # initialize panning attributes
        self.press = None
        self.xpress = None
        self.ypress = None
        self.cur_xlim = None
        self.cur_ylim = None

        # set attribute holding data
        if data is None:
            try:
                data = parent.get_data ().getWFData ().copy ()
                self.data = data.select (station=station)
            except AttributeError as e:
                errmsg = 'You either have to put in a data or an appropriate ' \
                         'parent (PyLoT MainWindow) object: {0}'.format (e)
                raise Exception (errmsg)
        else:
            self.data = data

        self.stime, self.etime = full_range (self.getWFData ())

        # initialize plotting widget
        self.multicompfig = WaveformWidget (self)
        self.phaseplot = PhasePlotWidget(self)
        self.phaseplot.hide()

        # setup ui
        self.setupUi ()

        # plot data
        self.getPlotWidget ().plotWFData (wfdata=self.getWFData (),
                                          title=self.getStation ())

        xlims = self.getPlotWidget ().getXLims ()
        ylims = self.getPlotWidget ().getYLims ()

        self.limits = {'x': xlims,
                       'y': ylims}

        self.updateCurrentLimits ()

        # set plot labels
        self.setPlotLabels ()

        # draw picks if present
        self.drawAllPicks ()

        # connect button press event to an action
        self.cidpress = self.connectPressEvent (self.panPress)
        self.cidmotion = self.connectMotionEvent (self.panMotion)
        self.cidrelease = self.connectReleaseEvent (self.panRelease)
        self.cidscroll = self.connectScrollEvent (self.scrollZoom)

        # init expected picks using obspy Taup
        try:
            if self.parent ().metadata:
                self.model = TauPyModel (model)
                self.get_arrivals ()
                self.drawArrivals ()
                self.activateArrivalsButton (True)
            else:
                self.activateArrivalsButton (False)
        except Exception as e:
            print ('Warning: Could not init expected picks from taup: {}'.format (e))
            self.activateArrivalsButton (False)

        # init pick delete (with right click)
        self.connect_pick_delete ()

    def setupUi(self):
        menuBar = QtGui.QMenuBar (self)
        if not self._embedded:
            exitMenu = menuBar.addMenu ('File')
            exitAction = QtGui.QAction ('Close', self)
            exitAction.triggered.connect (self.close)
            exitMenu.addAction (exitAction)

        self.addPickPhases (menuBar)

        # create matplotlib toolbar to inherit functionality
        self.figToolBar = NavigationToolbar2QT (self.getPlotWidget (), self)
        self.figToolBar.hide ()

        # create icons
        filter_icon = QIcon ()
        filter_icon.addPixmap (QPixmap (':/icons/filter.png'))
        zoom_icon = QIcon ()
        zoom_icon.addPixmap (QPixmap (':/icons/zoom_in.png'))
        home_icon = QIcon ()
        home_icon.addPixmap (QPixmap (':/icons/zoom_0.png'))
        del_icon = QIcon ()
        del_icon.addPixmap (QPixmap (':/icons/delete.png'))

        # create actions
        self.filterAction = createAction (parent=self, text='Filter',
                                          slot=self.filterWFData,
                                          icon=filter_icon,
                                          tip='Toggle filtered/original'
                                              ' waveforms')
        self.zoomAction = createAction (parent=self, text='Zoom',
                                        slot=self.zoom, icon=zoom_icon,
                                        tip='Zoom into waveform',
                                        checkable=True)
        self.resetZoomAction = createAction (parent=self, text='Home',
                                             slot=self.resetZoom, icon=home_icon,
                                             tip='Reset zoom to original limits')
        self.resetPicksAction = createAction (parent=self, text='Delete Picks',
                                              slot=self.delPicks, icon=del_icon,
                                              tip='Delete current picks.')

        # create other widget elements
        phaseitems = [None] + list (FILTERDEFAULTS.keys ())

        # create buttons for P and S filter and picking
        self.p_button = QPushButton ('P', self)
        self.s_button = QPushButton ('S', self)
        self.p_button.setCheckable (True)
        self.s_button.setCheckable (True)
        # set button tooltips
        # self.p_button.setToolTip('Hotkey: "1"')
        # self.s_button.setToolTip('Hotkey: "2"')

        self.plot_arrivals_button = QPushButton ('Plot phases')
        self.plot_arrivals_button.setCheckable(True)

        # create accept/reject button
        self.accept_button = QPushButton ('&Accept Picks')
        self.reject_button = QPushButton ('&Reject Picks')
        self.disable_ar_buttons ()

        # add hotkeys
        self._shortcut_space = QtGui.QShortcut (QtGui.QKeySequence (' '), self)
        self._shortcut_space.activated.connect (self.accept_button.clicked)
        # button shortcuts (1 for P-button, 2 for S-button)
        # self.p_button.setShortcut(QKeySequence('1'))
        # self.s_button.setShortcut(QKeySequence('2'))

        # layout the outermost appearance of the Pick Dialog
        _outerlayout = QVBoxLayout ()
        _dialtoolbar = QToolBar ()

        # fill toolbar with content
        _dialtoolbar.addAction (self.filterAction)
        _dialtoolbar.addWidget (self.p_button)
        _dialtoolbar.addWidget (self.s_button)
        _dialtoolbar.addAction (self.zoomAction)
        _dialtoolbar.addSeparator ()
        _dialtoolbar.addAction (self.resetZoomAction)
        _dialtoolbar.addSeparator ()
        _dialtoolbar.addAction (self.resetPicksAction)
        if self._embedded:
            _dialtoolbar.addWidget (self.accept_button)
            _dialtoolbar.addWidget (self.reject_button)
        else:
            _dialtoolbar.addWidget (self.nextStation)
        _dialtoolbar.addWidget (self.plot_arrivals_button)

        # layout the innermost widget
        _innerlayout = QVBoxLayout ()
        _innerinnerlayout = QtGui.QHBoxLayout()
        _innerinnerlayout.addWidget (self.multicompfig)
        _innerinnerlayout.addWidget (self.phaseplot)
        _innerlayout.addLayout(_innerinnerlayout)

        # add button box to the dialog
        _buttonbox = QDialogButtonBox (QDialogButtonBox.Ok |
                                       QDialogButtonBox.Cancel)

        # merge widgets and layouts to establish the dialog
        if not self._embedded:
            _innerlayout.addWidget (_buttonbox)
        _outerlayout.addWidget (menuBar)
        _outerlayout.addWidget (_dialtoolbar)
        _outerlayout.addLayout (_innerlayout)
        _outerlayout.setStretch (0, 0)
        _outerlayout.setStretch (1, 0)
        _outerlayout.setStretch (2, 1)

        # connect widget element signals with slots (methods to the dialog
        # object
        self.p_button.clicked.connect (self.p_clicked)
        self.s_button.clicked.connect (self.s_clicked)
        self.accept_button.clicked.connect (self.accept)
        self.reject_button.clicked.connect (self.reject)
        self.accept_button.clicked.connect (self.disable_ar_buttons)
        self.reject_button.clicked.connect (self.disable_ar_buttons)
        self.plot_arrivals_button.clicked.connect (self.toggle_arrivals_plot)
        _buttonbox.accepted.connect (self.accept)
        _buttonbox.rejected.connect (self.reject)

        # finally layout the entire dialog
        self.setLayout (_outerlayout)
        self.resize (1280, 720)

    def activateArrivalsButton(self, val=True):
        self.plot_arrivals_button.setEnabled (val)

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
            ax.legend()
            self.phaseplot.new = False
            self.phaseplot.draw()
        self.phaseplot.show()

    def setDirty(self, bool):
        self._dirty = bool

    def get_arrivals(self, plot=False):
        func = {True: self.model.get_ray_paths_geo,
                False: self.model.get_travel_times_geo}
        phases = self.prepare_phases ()
        station_id = self.data.traces[0].get_id ()
        parser = self.parent ().metadata[1]
        station_coords = parser.get_coordinates (station_id)
        origins = self.parent ().get_current_event ().origins
        if origins:
            source_origin = origins[0]
        else:
            raise ValueError ('No source origin given.')
        arrivals = func[plot] (source_origin.depth,
                               source_origin.latitude,
                               source_origin.longitude,
                               station_coords['latitude'],
                               station_coords['longitude'],
                               phases)
        self.arrivals = arrivals

    def prepare_phases(self):
        settings = QtCore.QSettings ()
        p_phases = settings.value ('p_phases')
        s_phases = settings.value ('s_phases')
        phases = p_phases + ',' + s_phases
        phases = phases.split (',')
        phases = [phase.strip () for phase in phases]
        return phases

    def drawArrivals(self, textOnly=False):
        if not self.arrivals:
            return
        ax = self.getPlotWidget ().axes
        if not textOnly:
            ylims = self.getGlobalLimits ('y')
        else:
            ylims = self.getPlotWidget ().getYLims ()
        stime = self.getStartTime ()
        source_origin = self.parent ().get_current_event ().origins[0]
        source_time = source_origin.time
        for arrival in self.arrivals:
            arrival_time_abs = source_time + arrival.time
            time_rel = arrival_time_abs - stime
            if not textOnly:
                ax.plot ([time_rel, time_rel], ylims, '0.3', linestyle='dashed')
            self.arrivalsText.append (ax.text (time_rel, ylims[0], arrival.name, color='0.5'))

    def drawArrivalsText(self):
        return self.drawArrivals (True)

    def refreshArrivalsText(self):
        self.removeArrivalsText ()
        self.drawArrivalsText ()

    def removeArrivalsText(self):
        for textItem in self.arrivalsText:
            try:
                textItem.remove ()
            except:
                pass
        self.arrivalsText = []

    def addPickPhases(self, menuBar):
        settings = QtCore.QSettings ()
        p_phases = settings.value ('p_phases')
        s_phases = settings.value ('s_phases')

        if p_phases:
            p_phases = p_phases.split (',')
        else:
            p_phases = []
        if s_phases:
            s_phases = s_phases.split (',')
        else:
            s_phases = []

        phases = {'P': p_phases,
                  'S': s_phases}
        if not 'P' in phases['P'] and not 'p' in phases['P']:
            phases['P'] = ['P'] + phases['P']
        if not 'S' in phases['S'] and not 's' in phases['S']:
            phases['S'] = ['S'] + phases['S']

        picksMenu = menuBar.addMenu ('Picks')
        self.picksActions = {}

        # dictionary points on corresponding phase_select function
        phaseSelect = {'P': self.p_phase_select,
                       'S': self.s_phase_select}

        nHotkey = 4  # max hotkeys per phase
        hotkey = 1  # start hotkey

        # loop over P and S (use explicit list instead of iter over dict.keys to keep order)
        for phaseIndex, phaseID in enumerate (['P', 'S']):
            # loop through phases in list
            for index, phase in enumerate (phases[phaseID]):
                # remove zeros
                phase = phase.strip ()
                # add hotkeys
                if not index >= nHotkey:
                    shortcut = str (hotkey)
                    hotkey += 1
                else:
                    shortcut = None
                # create action and add to menu
                # phase name transferred using lambda function
                slot = lambda phase=phase, phaseID=phaseID: phaseSelect[phaseID] (phase)
                picksAction = createAction (parent=self, text=phase,
                                            slot=slot,
                                            shortcut=shortcut)
                picksMenu.addAction (picksAction)
                self.picksActions[str (phase)] = picksAction  # save action in dictionary
            if phaseIndex == 0:
                picksMenu.addSeparator ()

    def disconnectPressEvent(self):
        widget = self.getPlotWidget ()
        widget.mpl_disconnect (self.cidpress)
        self.cidpress = None

    def connectPressEvent(self, slot):
        widget = self.getPlotWidget ()
        return widget.mpl_connect ('button_press_event', slot)

    def disconnectScrollEvent(self):
        widget = self.getPlotWidget ()
        widget.mpl_disconnect (self.cidscroll)
        self.cidscroll = None

    def connectScrollEvent(self, slot):
        widget = self.getPlotWidget ()
        return widget.mpl_connect ('scroll_event', slot)

    def disconnectMotionEvent(self):
        widget = self.getPlotWidget ()
        widget.mpl_disconnect (self.cidmotion)
        self.cidmotion = None

    def connectMotionEvent(self, slot):
        widget = self.getPlotWidget ()
        return widget.mpl_connect ('motion_notify_event', slot)

    def disconnectReleaseEvent(self):
        widget = self.getPlotWidget ()
        widget.mpl_disconnect (self.cidrelease)
        self.cidrelease = None

    def connectReleaseEvent(self, slot):
        widget = self.getPlotWidget ()
        return widget.mpl_connect ('button_release_event', slot)

    def disable_ar_buttons(self):
        self.enable_ar_buttons (False)

    def enable_ar_buttons(self, bool=True):
        self.accept_button.setEnabled (bool)
        self.reject_button.setEnabled (bool)

    def p_phase_select(self, phase):
        if not self.p_button.isChecked ():
            self.p_button.setChecked (True)
            self.p_button.setText (phase)
        else:
            if str (phase) == str (self.p_button.text ()):
                self.reset_p_button ()
            else:
                self.p_button.setText (phase)
        self.p_clicked ()

    def s_phase_select(self, phase):
        if not self.s_button.isChecked ():
            self.s_button.setChecked (True)
            self.s_button.setText (phase)
        else:
            if str (phase) == str (self.s_button.text ()):
                self.reset_s_button ()
            else:
                self.s_button.setText (phase)
        self.s_clicked ()

    def p_clicked(self):
        if self.p_button.isChecked ():
            self.reset_s_button ()
            self.s_button.setEnabled (False)
            self.init_p_pick ()
        else:
            self.leave_picking_mode ()

    def s_clicked(self):
        if self.s_button.isChecked ():
            self.reset_p_button ()
            self.p_button.setEnabled (False)
            self.init_s_pick ()
        else:
            self.leave_picking_mode ()

    def init_p_pick(self):
        self.set_button_color (self.p_button, 'yellow')
        self.updateCurrentLimits ()
        self.activatePicking ()
        self.currentPhase = str(self.p_button.text ())

    def init_s_pick(self):
        self.set_button_color (self.s_button, 'yellow')
        self.updateCurrentLimits ()
        self.activatePicking ()
        self.currentPhase = str(self.s_button.text ())

    def getPhaseID(self, phase):
        return identifyPhase (loopIdentifyPhase (phase))

    def set_button_color(self, button, color=None):
        if type (color) == QtGui.QColor:
            palette = button.palette ()
            role = button.backgroundRole ()
            palette.setColor (role, color)
            button.setPalette (palette)
            button.setAutoFillBackground (True)
        elif type (color) == str or not color:
            button.setStyleSheet ("background-color: {}".format (color))

    def reset_p_button(self):
        self.set_button_color (self.p_button)
        self.p_button.setEnabled (True)
        self.p_button.setChecked (False)
        self.p_button.setText ('P')

    def reset_s_button(self):
        self.set_button_color (self.s_button)
        self.s_button.setEnabled (True)
        self.s_button.setChecked (False)
        self.s_button.setText ('S')

    def leave_picking_mode(self):
        self.currentPhase = None
        self.reset_p_button ()
        self.reset_s_button ()
        self.getPlotWidget ().plotWFData (wfdata=self.getWFData (),
                                          title=self.getStation ())
        self.drawAllPicks ()
        self.drawArrivals ()
        self.setPlotLabels ()
        self.resetZoomAction.trigger ()
        self.deactivatePicking ()

    def activatePicking(self):
        if self.zoomAction.isChecked ():
            self.zoomAction.trigger ()
        self.disconnectReleaseEvent ()
        self.disconnectScrollEvent ()
        self.disconnectMotionEvent ()
        self.disconnectPressEvent ()
        self.cidpress = self.connectPressEvent (self.setIniPick)
        self.filterWFData ()
        self.pick_block = self.togglePickBlocker ()
        self.disconnect_pick_delete ()

    def deactivatePicking(self):
        self.disconnectPressEvent ()
        self.cidpress = self.connectPressEvent (self.panPress)
        self.cidmotion = self.connectMotionEvent (self.panMotion)
        self.cidrelease = self.connectReleaseEvent (self.panRelease)
        self.cidscroll = self.connectScrollEvent (self.scrollZoom)
        self.connect_pick_delete ()

    def getParameter(self):
        return self.parameter

    def getStartTime(self):
        return self.stime

    def getEndTime(self):
        return self.etime

    def getComponents(self):
        return self.components

    def getStation(self):
        if self.network and self.station:
            return self.network + '.' + self.station
        return self.station

    def getPlotWidget(self):
        return self.multicompfig

    def getChannelID(self, key):
        return self.getPlotWidget ().getPlotDict ()[int (key)][1]

    def getTraceID(self, channels):
        plotDict = self.getPlotWidget ().getPlotDict ()
        traceIDs = []
        for channel in channels:
            channel = channel.upper ()
            for traceID, channelID in plotDict.items ():
                if channelID[1].upper ().endswith (channel):
                    traceIDs.append (traceID)
        return traceIDs

    def getUser(self):
        return self._user

    def getFilterOptions(self, phase):
        options = self.filteroptions[self.getPhaseID(phase)]
        if type (options) == dict:
            return FilterOptions (**options)
        else:
            return options

    def getXLims(self):
        return self.cur_xlim

    def getYLims(self):
        return self.cur_ylim

    def setXLims(self, limits):
        self.cur_xlim = limits

    def setYLims(self, limits):
        self.cur_ylim = limits

    def getGlobalLimits(self, axis):
        return self.limits[axis]

    def updateCurrentLimits(self):
        self.setXLims (self.getPlotWidget ().getXLims ())
        self.setYLims (self.getPlotWidget ().getYLims ())

    def getWFData(self):
        return self.data

    def selectWFData(self, channel):
        component = channel[-1].upper ()
        wfdata = Stream ()

        def selectTrace(tr, components):
            if tr.stats.channel[-1].upper () in components:
                return tr

        if component == 'E' or component == 'N':
            for trace in self.getWFData ():
                trace = selectTrace (trace, 'NE')
                if trace:
                    wfdata.append (trace)
        elif component == '1' or component == '2':
            for trace in self.getWFData ():
                trace = selectTrace (trace, '12')
                if trace:
                    wfdata.append (trace)
        else:
            wfdata = self.getWFData ().select (component=component)
        return wfdata

    def getPicks(self, picktype='manual'):
        if picktype == 'manual':
            return self.picks
        elif picktype == 'auto':
            return self.autopicks
        else:
            raise TypeError ('Unknown picktype {0}'.format (picktype))

    def resetPicks(self):
        self.picks = {}

    def delPicks(self):
        self.resetPicks ()
        self.resetPlot ()

    def setIniPick(self, gui_event):

        trace_number = round (gui_event.ydata)

        channel = self.getChannelID (trace_number)
        wfdata = self.selectWFData (channel)

        self.disconnectScrollEvent ()
        self.disconnectPressEvent ()
        self.disconnectReleaseEvent ()
        self.disconnectMotionEvent ()
        self.cidpress = self.connectPressEvent (self.setPick)

        if self.getPhaseID(self.currentPhase) == 'P':
            self.set_button_color (self.p_button, 'green')
            self.setIniPickP (gui_event, wfdata, trace_number)
        elif self.getPhaseID(self.currentPhase) == 'S':
            self.set_button_color (self.s_button, 'green')
            self.setIniPickS (gui_event, wfdata)

        self.zoomAction.setEnabled (False)

        # reset labels
        self.setPlotLabels ()
        self.draw ()

    def setIniPickP(self, gui_event, wfdata, trace_number):

        parameter = self.parameter
        ini_pick = gui_event.xdata

        nfac = parameter.get ('nfacP')
        twins = parameter.get ('tsnrz')
        noise_win = twins[0]
        gap_win = twins[1]
        signal_win = twins[2]
        itrace = int (trace_number)

        while itrace > len (wfdata) - 1:
            itrace -= 1

        stime = self.getStartTime ()
        stime_diff = wfdata[itrace].stats.starttime - stime

        # copy data for plotting
        data = self.getWFData ().copy ()

        # filter data and trace on which is picked prior to determination of SNR
        phase = self.currentPhase
        filteroptions = self.getFilterOptions (self.getPhaseID(phase)).parseFilterOptions ()
        if filteroptions:
            try:
                data.filter (**filteroptions)
                wfdata.filter (**filteroptions)
            except ValueError as e:
                self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Information,
                                              'Denied', 'setIniPickP: Could not filter waveform: {}'.format (e))
                self.qmb.show ()

        result = getSNR (wfdata, (noise_win, gap_win, signal_win), ini_pick - stime_diff, itrace)

        snr = result[0]
        noiselevel = result[2]
        if noiselevel:
            noiselevel *= nfac
        else:
            noiselevel = nfac

        x_res = getResolutionWindow (snr, parameter.get ('extent'))

        # remove mean noise level from waveforms
        for trace in data:
            t = prepTimeAxis (trace.stats.starttime - stime, trace)
            inoise = getnoisewin (t, ini_pick, noise_win, gap_win)
            trace = demeanTrace (trace=trace, window=inoise)

        self.setXLims ([ini_pick - x_res, ini_pick + x_res])
        self.setYLims (np.array ([-noiselevel * 3.5, noiselevel * 3.5]) +
                       trace_number)
        self.getPlotWidget ().plotWFData (wfdata=data,
                                          title=self.getStation () +
                                                ' picking mode',
                                          zoomx=self.getXLims (),
                                          zoomy=self.getYLims (),
                                          noiselevel=(trace_number + noiselevel,
                                                      trace_number - noiselevel),
                                          iniPick=ini_pick)

    def setIniPickS(self, gui_event, wfdata):

        parameter = self.parameter
        ini_pick = gui_event.xdata

        nfac = parameter.get ('nfacS')
        twins = parameter.get ('tsnrh')
        noise_win = twins[0]
        gap_win = twins[1]
        signal_win = twins[2]

        stime = self.getStartTime ()
        stime_diff = wfdata[0].stats.starttime - stime

        # copy data for plotting
        data = self.getWFData ().copy ()

        # filter data and trace on which is picked prior to determination of SNR
        phase = self.currentPhase
        filteroptions = self.getFilterOptions (self.getPhaseID(phase)).parseFilterOptions ()
        if filteroptions:
            try:
                data.filter (**filteroptions)
                wfdata.filter (**filteroptions)
            except ValueError as e:
                self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Information,
                                              'Denied', 'setIniPickS: Could not filter waveform: {}'.format (e))
                self.qmb.show ()

        # determine SNR and noiselevel
        result = getSNR (wfdata, (noise_win, gap_win, signal_win), ini_pick - stime_diff)
        snr = result[0]
        noiselevel = result[2]

        if noiselevel:
            noiselevel *= nfac
        else:
            noiselevel = nfac

        # prepare plotting of data
        for trace in data:
            t = prepTimeAxis (trace.stats.starttime - stime, trace)
            inoise = getnoisewin (t, ini_pick, noise_win, gap_win)
            trace = demeanTrace (trace, inoise)

        # scale waveform for plotting
        horiz_comp = find_horizontals (data)
        data = scaleWFData (data, noiselevel * 2.5, horiz_comp)

        x_res = getResolutionWindow (snr, parameter.get ('extent'))

        self.setXLims (tuple ([ini_pick - x_res, ini_pick + x_res]))
        traces = self.getTraceID (horiz_comp)
        traces.sort ()
        self.setYLims (tuple (np.array ([-1.0, +1.0]) +
                              np.array (traces)))
        noiselevels = [trace + 1 / (2.5 * 2) for trace in traces] + \
                      [trace - 1 / (2.5 * 2) for trace in traces]

        self.getPlotWidget ().plotWFData (wfdata=data,
                                          title=self.getStation () +
                                                ' picking mode',
                                          zoomx=self.getXLims (),
                                          zoomy=self.getYLims (),
                                          noiselevel=noiselevels,
                                          scaleddata=True,
                                          iniPick=ini_pick)

    def setPick(self, gui_event):

        parameter = self.parameter

        # get axes limits
        self.updateCurrentLimits ()

        # setting pick
        pick = gui_event.xdata  # get pick time relative to the traces timeaxis not to the global
        channel = self.getChannelID (round (gui_event.ydata))

        # get name of phase actually picked
        phase = self.currentPhase

        # get filter parameter for the phase to be picked
        filteroptions = self.getFilterOptions (self.getPhaseID(phase)).parseFilterOptions ()

        # copy and filter data for earliest and latest possible picks
        wfdata = self.getWFData ().copy ().select (channel=channel)
        if filteroptions:
            try:
                wfdata.filter (**filteroptions)
            except ValueError as e:
                self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Information,
                                              'Denied', 'setPick: Could not filter waveform: {}'.format (e))
                self.qmb.show ()

        # get earliest and latest possible pick and symmetric pick error
        if wfdata[0].stats.channel[2] == 'Z' or wfdata[0].stats.channel[2] == '3':
            nfac = parameter.get ('nfacP')
            TSNR = parameter.get ('tsnrz')
        else:
            nfac = parameter.get ('nfacS')
            TSNR = parameter.get ('tsnrh')

        # return absolute time values for phases
        stime = self.getStartTime ()
        stime_diff = wfdata[0].stats.starttime - stime

        [epp, lpp, spe] = earllatepicker (wfdata, nfac, (TSNR[0], TSNR[1], TSNR[2]),
                                          pick - stime_diff, verbosity=1)

        mpp = stime + pick
        if epp:
            epp = stime + epp + stime_diff
        if lpp:
            lpp = stime + lpp + stime_diff

        # save pick times for actual phase
        phasepicks = dict (epp=epp, lpp=lpp, mpp=mpp, spe=spe,
                           picker='manual', channel=channel,
                           network=wfdata[0].stats.network)

        try:
            oldphasepick = self.picks[phase]
        except KeyError:
            self.picks[phase] = phasepicks
        else:
            self.picks[phase] = phasepicks
            oepp = oldphasepick['epp']
            ompp = oldphasepick['mpp']
            olpp = oldphasepick['lpp']
            msg = """Warning old phase information for phase {phase} has been
                     altered.\n
                     New phase times:\n
                     earliest possible pick: {epp}\n
                     most probable pick: {mpp}\n
                     latest possible pick: {lpp}\n
                     \n
                     Old phase times (overwritten):\n
                     earliest possible pick: {oepp}\n
                     most probable pick: {ompp}\n
                     latest possible pick: {olpp}\n""".format (phase=phase,
                                                               epp=epp,
                                                               mpp=pick,
                                                               lpp=lpp,
                                                               oepp=oepp,
                                                               ompp=ompp,
                                                               olpp=olpp)

        self.disconnectPressEvent ()
        self.enable_ar_buttons ()
        self.zoomAction.setEnabled (True)
        self.pick_block = self.togglePickBlocker ()
        self.leave_picking_mode ()
        self.setDirty (True)

    def drawAllPicks(self):
        self.removePhaseText ()
        self.drawPicks (picktype='manual')
        self.drawPicks (picktype='auto')

    def drawPicks(self, phase=None, picktype='manual', textOnly=False):
        # plotting picks
        ax = self.getPlotWidget ().axes
        if not textOnly:
            ylims = self.getGlobalLimits ('y')
        else:
            ylims = self.getPlotWidget ().getYLims ()
        phase_col = {
            'P': ('c', 'c--', 'b-', 'bv', 'b^', 'b', 'c:'),
            'S': ('m', 'm--', 'r-', 'rv', 'r^', 'r', 'm:')
        }
        if self.getPicks (picktype):
            if phase is not None:
                if (type (self.getPicks (picktype)[phase]) is dict
                    or type (self.getPicks (picktype)[phase]) is AttribDict):
                    picks = self.getPicks (picktype)[phase]
                    colors = phase_col[self.getPhaseID(phase)]
            elif phase is None:
                for phase in self.getPicks (picktype):
                    self.drawPicks (phase, picktype, textOnly)
                return
            else:
                return
        else:
            return

        mpp = picks['mpp'] - self.getStartTime ()
        if picks['epp'] and picks['lpp'] and not textOnly:
            epp = picks['epp'] - self.getStartTime ()
            lpp = picks['lpp'] - self.getStartTime ()
        spe = picks['spe']

        if picktype == 'manual':
            if not textOnly:
                if picks['epp'] and picks['lpp']:
                    ax.fill_between ([epp, lpp], ylims[0], ylims[1],
                                     alpha=.25, color=colors[0], label='EPP, LPP')
                if spe:
                    ax.plot ([mpp - spe, mpp - spe], ylims, colors[1], label='{}-SPE'.format (phase))
                    ax.plot ([mpp + spe, mpp + spe], ylims, colors[1])
                    ax.plot ([mpp, mpp], ylims, colors[2], label='{}-Pick'.format (phase), picker=5)
                else:
                    ax.plot ([mpp, mpp], ylims, colors[6], label='{}-Pick (NO PICKERROR)'.format (phase), picker=5)
            # append phase text (if textOnly: draw with current ylims)
            self.phaseText.append (ax.text (mpp, ylims[1], phase))
        elif picktype == 'auto':
            ax.plot (mpp, ylims[1], colors[3],
                     mpp, ylims[0], colors[4])
            ax.vlines (mpp, ylims[0], ylims[1], colors[5], linestyles='dotted')
        else:
            raise TypeError ('Unknown picktype {0}'.format (picktype))

        ax.legend ()

    def connect_pick_delete(self):
        self.cidpick = self.getPlotWidget ().mpl_connect ('pick_event', self.onpick_delete)

    def disconnect_pick_delete(self):
        if hasattr (self, 'cidpick'):
            self.getPlotWidget ().mpl_disconnect (self.cidpick)

    def onpick_delete(self, event):
        if not event.mouseevent.button == 3:
            return
        x = event.mouseevent.xdata
        self.remove_pick_by_x (x)
        self.resetPlot ()

    def remove_pick_by_x(self, x):
        if not self.picks and not self.autopicks:
            return
        # init empty list and get station starttime
        X = []
        starttime = self.getStartTime ()
        # init dictionarys to iterate through and iterate over them
        allpicks = [self.picks, self.autopicks]
        for index_ptype, picks in enumerate (allpicks):
            for phase in picks:
                pick_rel = picks[phase]['mpp'] - starttime
                # add relative pick time, phaseID and picktype index
                X.append ((pick_rel, phase, index_ptype))
        # find index and value closest to x
        index, value = min (enumerate ([val[0] for val in X]), key=lambda y: abs (y[1] - x))
        # unpack the found value
        pick_rel, phase, index_ptype = X[index]
        # delete the value from corresponding dictionary
        allpicks[index_ptype].pop (phase)
        # information output
        msg = 'Deleted pick for phase {}, {}[s] from starttime {}'
        print (msg.format (phase, pick_rel, starttime))
        self.setDirty (True)

    def drawPhaseText(self):
        return self.drawPicks (picktype='manual', textOnly=True)

    def removePhaseText(self):
        for textItem in self.phaseText:
            try:
                textItem.remove ()
            except:
                pass
        self.phaseText = []

    def refreshPhaseText(self):
        self.removePhaseText ()
        self.drawPhaseText ()

    def panPress(self, gui_event):
        ax = self.getPlotWidget ().axes
        if gui_event.inaxes != ax: return
        self.cur_xlim = ax.get_xlim ()
        self.cur_ylim = ax.get_ylim ()
        self.press = gui_event.xdata, gui_event.ydata
        self.xpress, self.ypress = self.press

    def panRelease(self, gui_event):
        ax = self.getPlotWidget ().axes
        self.press = None
        self.refreshPhaseText ()
        self.refreshArrivalsText ()
        ax.figure.canvas.draw ()

    def panMotion(self, gui_event):
        if self.press is None: return
        ax = self.getPlotWidget ().axes
        if gui_event.inaxes != ax: return
        dx = gui_event.xdata - self.xpress
        dy = gui_event.ydata - self.ypress
        self.cur_xlim -= dx
        self.cur_ylim -= dy
        ax.set_xlim (self.cur_xlim)
        ax.set_ylim (self.cur_ylim)

        ax.figure.canvas.draw ()

    def togglePickBlocker(self):
        return not self.pick_block

    def filterWFData(self):
        if self.pick_block:
            return
        self.updateCurrentLimits ()
        data = self.getWFData ().copy ()
        old_title = self.getPlotWidget ().getAxes ().get_title ()
        title = None
        phase = self.currentPhase
        filtoptions = None
        if phase:
            filtoptions = self.getFilterOptions (self.getPhaseID(phase)).parseFilterOptions ()
        if self.filterAction.isChecked ():
            if not phase:
                filtoptions = FilterOptionsDialog.getFilterObject ()
                filtoptions = filtoptions.parseFilterOptions ()
        if filtoptions is not None:
            data.filter (**filtoptions)
            if not old_title.endswith (')'):
                title = old_title + ' (filtered)'
            elif not old_title.endswith (' (filtered)') and not old_title.endswith (', filtered)'):
                title = old_title[:-1] + ', filtered)'
        else:
            if old_title.endswith (' (filtered)'):
                title = old_title.replace (' (filtered)', '')
            elif old_title.endswith (', filtered)'):
                title = old_title.replace (', filtered)', ')')
        if title is None:
            title = old_title
        self.getPlotWidget ().plotWFData (wfdata=data, title=title,
                                          zoomx=self.getXLims (),
                                          zoomy=self.getYLims ())
        self.setPlotLabels ()
        self.drawAllPicks ()
        self.draw ()

    def resetPlot(self):
        self.updateCurrentLimits ()
        data = self.getWFData ().copy ()
        title = self.getPlotWidget ().getAxes ().get_title ()
        self.getPlotWidget ().plotWFData (wfdata=data, title=title,
                                          zoomx=self.getXLims (),
                                          zoomy=self.getYLims ())
        self.setPlotLabels ()
        self.drawAllPicks ()
        self.draw ()

    def setPlotLabels(self):

        # get channel labels
        pos = self.getPlotWidget ().getPlotDict ().keys ()
        labels = [self.getPlotWidget ().getPlotDict ()[key][1] for key in pos]

        # set channel labels
        self.getPlotWidget ().setYTickLabels (pos, labels)
        self.getPlotWidget ().setXLims (self.getXLims ())
        self.getPlotWidget ().setYLims (self.getYLims ())

    def zoom(self):
        if self.zoomAction.isChecked () and self.pick_block:
            self.zoomAction.setChecked (False)
        elif self.zoomAction.isChecked ():
            self.disconnectPressEvent ()
            self.disconnectMotionEvent ()
            self.disconnectReleaseEvent ()
            self.disconnectScrollEvent ()
            self.figToolBar.zoom ()
        else:
            self.figToolBar.zoom ()
            self.cidpress = self.connectPressEvent (self.panPress)
            self.cidmotion = self.connectMotionEvent (self.panMotion)
            self.cidrelease = self.connectReleaseEvent (self.panRelease)
            self.cidscroll = self.connectScrollEvent (self.scrollZoom)

    def scrollZoom(self, gui_event, factor=2.):

        self.updateCurrentLimits ()

        if gui_event.button == 'up':
            scale_factor = 1 / factor
        elif gui_event.button == 'down':
            # deal with zoom out
            scale_factor = factor
        else:
            # deal with something that should never happen
            scale_factor = 1
            print (gui_event.button)

        new_xlim = gui_event.xdata - \
                   scale_factor * (gui_event.xdata - self.getXLims ())
        new_ylim = gui_event.ydata - \
                   scale_factor * (gui_event.ydata - self.getYLims ())

        new_xlim.sort ()
        global_x = self.getGlobalLimits ('x')
        global_y = self.getGlobalLimits ('y')
        new_xlim[0] = max (new_xlim[0], global_x[0])
        new_xlim[1] = min (new_xlim[1], global_x[1])
        new_ylim.sort ()
        new_ylim[0] = max (new_ylim[0], global_y[0])
        new_ylim[1] = min (new_ylim[1], global_y[1])

        self.getPlotWidget ().setXLims (new_xlim)
        self.getPlotWidget ().setYLims (new_ylim)
        self.refreshArrivalsText ()
        self.refreshPhaseText ()
        self.draw ()

    def resetZoom(self):
        self.getPlotWidget ().setXLims (self.getGlobalLimits ('x'))
        self.getPlotWidget ().setYLims (self.getGlobalLimits ('y'))
        self.draw ()

    def draw(self):
        self.getPlotWidget ().draw ()

    def apply(self):
        picks = self.getPicks ()
        self.update_picks.emit (picks)
        # for pick in picks:
        #     print(pick, picks[pick])

    def discard(self):
        picks = self._init_picks
        self.picks = picks
        self.update_picks.emit (picks)
        # for pick in picks:
        #     print(pick, picks[pick])

    def reject(self):
        self.discard ()
        if not self._embedded:
            QDialog.reject (self)
        else:
            self.resetPlot ()
            self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Information,
                                          'Denied', 'New picks rejected!')
            self.qmb.show ()

    def accept(self):
        self.apply ()
        if not self._embedded:
            QDialog.accept (self)
        else:
            self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Information,
                                          'Accepted', 'New picks applied!')
            self.qmb.show ()


class PhasePlotWidget(FigureCanvas):
    def __init__(self, parent=None):
        self._parent = parent
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111, projection='polar')
        self.new = True
        super (PhasePlotWidget, self).__init__ (self.fig)


class TuneAutopicker (QWidget):
    update = QtCore.Signal (str)
    '''
    QWidget used to modifiy and test picking parameters for autopicking algorithm.

    :param: parent
    :type: QtPyLoT Mainwindow
    '''

    def __init__(self, parent):
        QtGui.QWidget.__init__ (self, parent, 1)
        self.parent = parent
        self.setParent (parent)
        self.setWindowTitle ('PyLoT - Tune Autopicker')
        self.parameter = parent._inputs
        self.fig_dict = parent.fig_dict
        self.data = Data ()
        self.init_main_layouts ()
        self.init_eventlist ()
        self.init_figure_tabs ()
        self.init_stationlist ()
        self.init_pbwidget ()
        self.connect_signals ()
        self.add_parameters ()
        self.add_buttons ()
        self.add_log ()
        self.set_stretch ()
        self.resize (1280, 720)
        if hasattr (parent, 'metadata'):
            self.metadata = self.parent.metadata
        else:
            self.metadata = None
            # self.setWindowModality(QtCore.Qt.WindowModality.ApplicationModal)
            # self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)

    def init_main_layouts(self):
        self.main_layout = QtGui.QVBoxLayout ()
        self.tune_layout = QtGui.QHBoxLayout ()
        self.trace_layout = QtGui.QHBoxLayout ()
        self.parameter_layout = QtGui.QVBoxLayout ()

        self.main_layout.addLayout (self.trace_layout)
        self.main_layout.addLayout (self.tune_layout)
        self.setLayout (self.main_layout)

    def init_eventlist(self):
        self.eventBox = self.parent.createEventBox ()
        self.eventBox.setMaxVisibleItems (20)
        self.fill_eventbox ()
        self.trace_layout.addWidget (self.eventBox)

    def init_stationlist(self):
        self.stationBox = QtGui.QComboBox ()
        self.trace_layout.addWidget (self.stationBox)
        self.fill_stationbox ()
        self.figure_tabs.setCurrentIndex (0)

    def connect_signals(self):
        self.eventBox.activated.connect (self.fill_stationbox)
        self.eventBox.activated.connect (self.update_eventID)
        self.eventBox.activated.connect (self.fill_tabs)
        self.stationBox.activated.connect (self.fill_tabs)

    def fill_stationbox(self):
        fnames = self.parent.getWFFnames_from_eventbox (eventbox=self.eventBox)
        self.data.setWFData (fnames)
        self.stationBox.clear ()
        stations = []
        for trace in self.data.getWFData ():
            station = trace.stats.station
            network = trace.stats.network
            ns_tup = (str (network), str (station))
            if not ns_tup in stations:
                stations.append (ns_tup)
        stations.sort ()
        model = self.stationBox.model ()
        for network, station in stations:
            item = QtGui.QStandardItem (network + '.' + station)
            if station in self.get_current_event ().pylot_picks:
                item.setBackground (self.parent._colors['ref'])
            model.appendRow (item)

    def init_figure_tabs(self):
        self.figure_tabs = QtGui.QTabWidget ()
        self.fill_figure_tabs ()

    def init_pbwidget(self):
        self.pb_widget = QtGui.QWidget ()

    def init_tab_names(self):
        self.ptb_names = ['aicFig', 'slength', 'checkZ4s', 'refPpick', 'el_Ppick', 'fm_picker']
        self.stb_names = ['aicARHfig', 'refSpick', 'el_S1pick', 'el_S2pick']

    def add_parameters(self):
        self.paraBox = PylotParaBox (self.parameter)
        self.paraBox.set_tune_mode (True)
        self.update_eventID ()
        self.parameter_layout.addWidget (self.paraBox)
        self.parameter_layout.addWidget (self.pb_widget)
        self.tune_layout.insertLayout (1, self.parameter_layout)

    def add_buttons(self):
        self.pick_button = QtGui.QPushButton ('Pick Trace')
        self.pick_button.clicked.connect (self.call_picker)
        self.close_button = QtGui.QPushButton ('Close')
        self.close_button.clicked.connect (self.hide)
        self.trace_layout.addWidget (self.pick_button)
        self.trace_layout.setStretch (0, 1)
        self.parameter_layout.addWidget (self.close_button)

    def add_log(self):
        self.listWidget = QtGui.QListWidget ()
        self.figure_tabs.insertTab (4, self.listWidget, 'log')

    def add_log_item(self, text):
        self.listWidget.addItem (text)
        self.listWidget.scrollToBottom ()

    def get_current_event(self):
        path = self.eventBox.currentText ()
        return self.parent.project.getEventFromPath (path)

    def get_current_event_name(self):
        return self.eventBox.currentText ().split ('/')[-1]

    def get_current_event_fp(self):
        return self.eventBox.currentText ()

    def get_current_event_picks(self, station):
        event = self.get_current_event ()
        if station in event.pylot_picks.keys ():
            return event.pylot_picks[station]

    def get_current_event_autopicks(self, station):
        event = self.get_current_event ()
        if event.pylot_autopicks:
            return event.pylot_autopicks[station]

    def get_current_station(self):
        return str (self.stationBox.currentText ()).split ('.')[-1]

    def gen_tab_widget(self, name, canvas):
        widget = QtGui.QWidget ()
        v_layout = QtGui.QVBoxLayout ()
        v_layout.addWidget (canvas)
        v_layout.addWidget (NavigationToolbar2QT (canvas, self))
        widget.setLayout (v_layout)
        return widget

    def gen_pick_dlg(self):
        if not self.get_current_event ():
            self.pickDlg = None
            return
        station = self.get_current_station ()
        data = self.data.getWFData ()
        pickDlg = PickDlg (self, data=data.select (station=station),
                           station=station, parameter=self.parameter,
                           picks=self.get_current_event_picks (station),
                           autopicks=self.get_current_event_autopicks (station),
                           embedded=True)
        pickDlg.update_picks.connect (self.picks_from_pickdlg)
        pickDlg.update_picks.connect (self.fill_eventbox)
        pickDlg.update_picks.connect (self.fill_stationbox)
        pickDlg.update_picks.connect (lambda: self.parent.setDirty (True))
        pickDlg.update_picks.connect (self.parent.enableSaveManualPicksAction)
        self.pickDlg = QtGui.QWidget ()
        hl = QtGui.QHBoxLayout ()
        self.pickDlg.setLayout (hl)
        hl.addWidget (pickDlg)

    def picks_from_pickdlg(self, picks=None):
        station = self.get_current_station ()
        replot = self.parent.addPicks (station, picks)
        self.get_current_event ().setPick (station, picks)
        if self.get_current_event () == self.parent.get_current_event ():
            if replot:
                self.parent.plotWaveformDataThread ()
                self.parent.drawPicks ()
            else:
                self.parent.drawPicks (station)
            self.parent.draw ()

    def plot_manual_picks_to_figs(self):
        picks = self.get_current_event_picks (self.get_current_station ())
        if not picks:
            return
        st = self.data.getWFData ()
        tr = st.select (station=self.get_current_station ())[0]
        starttime = tr.stats.starttime
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
        for p_ax in p_axes:
            axes = self.parent.fig_dict[p_ax[0]].axes
            if not axes:
                continue
            ax = axes[p_ax[1]]
            self.plot_manual_Ppick_to_ax (ax, (picks['P']['mpp'] - starttime))
        for s_ax in s_axes:
            axes = self.parent.fig_dict[s_ax[0]].axes
            if not axes:
                continue
            ax = axes[s_ax[1]]
            self.plot_manual_Spick_to_ax (ax, (picks['S']['mpp'] - starttime))

    def plot_manual_Ppick_to_ax(self, ax, pick):
        y_top = 0.9 * ax.get_ylim ()[1]
        y_bot = 0.9 * ax.get_ylim ()[0]
        ax.vlines (pick, y_bot, y_top,
                   color='teal', linewidth=2, label='manual P Onset')
        ax.plot ([pick - 0.5, pick + 0.5],
                 [y_bot, y_bot], linewidth=2, color='teal')
        ax.plot ([pick - 0.5, pick + 0.5],
                 [y_top, y_top], linewidth=2, color='teal')
        ax.legend ()

    def plot_manual_Spick_to_ax(self, ax, pick):
        y_top = 0.9 * ax.get_ylim ()[1]
        y_bot = 0.9 * ax.get_ylim ()[0]
        ax.vlines (pick, y_bot, y_top,
                   color='magenta', linewidth=2, label='manual S Onset')
        ax.plot ([pick - 0.5, pick + 0.5],
                 [y_bot, y_bot], linewidth=2, color='magenta')
        ax.plot ([pick - 0.5, pick + 0.5],
                 [y_top, y_top], linewidth=2, color='magenta')
        ax.legend ()

    def fill_tabs(self, event=None, picked=False):
        self.clear_all ()
        canvas_dict = self.parent.canvas_dict
        self.gen_pick_dlg ()
        self.overview = self.gen_tab_widget ('Overview', canvas_dict['mainFig'])
        id0 = self.figure_tabs.insertTab (0, self.pickDlg, 'Traces Plot')
        id1 = self.figure_tabs.insertTab (1, self.overview, 'Overview')
        id2 = self.figure_tabs.insertTab (2, self.p_tabs, 'P')
        id3 = self.figure_tabs.insertTab (3, self.s_tabs, 'S')
        if picked and self.get_current_event ():
            self.fill_p_tabs (canvas_dict)
            self.fill_s_tabs (canvas_dict)
            self.toggle_autopickTabs (bool (self.fig_dict['mainFig'].axes))
            self.plot_manual_picks_to_figs ()
        else:
            self.disable_autopickTabs ()
        try:
            main_fig.tight_layout ()
        except:
            pass
        self.figure_tabs.setCurrentIndex (0)

    def fill_p_tabs(self, canvas_dict):
        for name in self.ptb_names:
            id = self.p_tabs.addTab (self.gen_tab_widget (name, canvas_dict[name]), name)
            self.p_tabs.setTabEnabled (id, bool (self.fig_dict[name].axes))
            try:
                self.fig_dict[name].tight_layout ()
            except:
                pass

    def fill_s_tabs(self, canvas_dict):
        for name in self.stb_names:
            id = self.s_tabs.addTab (self.gen_tab_widget (name, canvas_dict[name]), name)
            self.s_tabs.setTabEnabled (id, bool (self.fig_dict[name].axes))
            try:
                self.fig_dict[name].tight_layout ()
            except:
                pass

    def fill_figure_tabs(self):
        self.clear_all ()
        self.p_tabs = QtGui.QTabWidget ()
        self.s_tabs = QtGui.QTabWidget ()
        self.tune_layout.insertWidget (0, self.figure_tabs)
        self.init_tab_names ()

    def fill_eventbox(self):
        project = self.parent.project
        if not project:
            return
        # update own list
        self.parent.fill_eventbox (eventBox=self.eventBox, select_events='ref')
        index_start = self.parent.eventBox.currentIndex ()
        index = index_start
        if index == -1:
            index += 1
        nevents = self.eventBox.model ().rowCount ()
        path = self.eventBox.itemText (index)
        if project.getEventFromPath (path).isTestEvent ():
            for index in range (nevents):
                path = self.eventBox.itemText (index)
                if project.getEventFromPath (index):
                    if not project.getEventFromPath (index).isTestEvent ():
                        break
                # in case all events are marked as test events and last event is reached
                if index == nevents - 1:
                    index = -1
        self.eventBox.setCurrentIndex (index)
        if not index == index_start:
            self.eventBox.activated.emit (index)
        # update parent
        self.parent.fill_eventbox ()

    def update_eventID(self):
        self.paraBox.boxes['eventID'].setText (
            self.get_current_event_name ())
        self.figure_tabs.setCurrentIndex (0)

    def call_picker(self):
        self.parameter = self.params_from_gui ()
        station = self.get_current_station ()
        if not station:
            self._warn ('No station selected')
            return
        args = {'parameter': self.parameter,
                'station': station,
                'fnames': 'None',
                'eventid': self.get_current_event_fp (),
                'iplot': 2,
                'fig_dict': self.fig_dict,
                'locflag': 0}
        for key in self.fig_dict.keys ():
            self.fig_dict[key].clear ()
        self.ap_thread = Thread (self, autoPyLoT, arg=args,
                                 progressText='Picking trace...',
                                 pb_widget=self.pb_widget,
                                 redirect_stdout=True)
        self.enable (False)
        self.ap_thread.message.connect (self.add_log_item)
        self.ap_thread.finished.connect (self.finish_picker)
        self.figure_tabs.setCurrentIndex (4)
        self.ap_thread.start ()
        # picks = autoPyLoT(self.parameter, fnames='None', iplot=2, fig_dict=self.fig_dict)

    def finish_picker(self):
        self.enable (True)
        if not self.ap_thread._executed:
            self._warn ('Could not execute picker:\n{}'.format (
                self.ap_thread._executedError))
            return
        self.pylot_picks = self.ap_thread.data
        if not self.pylot_picks:
            self._warn ('No picks found. See terminal output.')
            return
        # renew tabs
        # self.fill_figure_tabs()
        self.set_stretch ()
        self.update.emit ('Update')
        self.figure_tabs.setCurrentIndex (1)

    def enable(self, bool):
        self.pick_button.setEnabled (bool)
        self.paraBox.setEnabled (bool)
        self.eventBox.setEnabled (bool)
        self.stationBox.setEnabled (bool)
        self.overview.setEnabled (bool)
        self.p_tabs.setEnabled (bool)
        self.s_tabs.setEnabled (bool)

    def params_from_gui(self):
        parameters = self.paraBox.params_from_gui ()
        if self.parent:
            self.parent._inputs = parameters
        return parameters

    def set_stretch(self):
        self.tune_layout.setStretch (0, 3)
        self.tune_layout.setStretch (1, 1)

    def clear_all(self):
        if hasattr (self, 'pickDlg'):
            if self.pickDlg:
                self.pickDlg.setParent (None)
            del (self.pickDlg)
        if hasattr (self, 'overview'):
            self.overview.setParent (None)
        if hasattr (self, 'p_tabs'):
            self.p_tabs.clear ()
            self.p_tabs.setParent (None)
        if hasattr (self, 's_tabs'):
            self.s_tabs.clear ()
            self.s_tabs.setParent (None)

    def disable_autopickTabs(self):
        self.toggle_autopickTabs (False)

    def toggle_autopickTabs(self, bool):
        self.figure_tabs.setTabEnabled (1, bool)
        self.figure_tabs.setTabEnabled (2, bool)
        self.figure_tabs.setTabEnabled (3, bool)

    def _warn(self, message):
        self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Warning,
                                      'Warning', message)
        self.qmb.show ()


class PylotParaBox (QtGui.QWidget):
    accepted = QtCore.Signal (str)
    rejected = QtCore.Signal (str)

    def __init__(self, parameter, parent=None):
        '''
        Generate Widget containing parameters for PyLoT.

        :param: parameter
        :type: PylotParameter (object)

        '''
        QtGui.QWidget.__init__ (self, parent)
        self.parameter = parameter
        self.tabs = QtGui.QTabWidget ()
        self.layout = QtGui.QVBoxLayout ()
        self._init_save_buttons ()
        self._init_tabs ()
        self._init_dialog_buttons ()
        self.labels = {}
        self.boxes = {}
        self.groupboxes = {}
        self._exclusive_widgets = []
        self._init_sublayouts ()
        self.setLayout (self.layout)
        self.add_main_parameters_tab ()
        self.add_special_pick_parameters_tab ()
        self.params_to_gui ()
        self._toggle_advanced_settings ()
        self.resize (720, 1280)
        self.setWindowModality (QtCore.Qt.WindowModality.ApplicationModal)
        self.accepted.connect (self.params_from_gui)
        self.rejected.connect (self.params_to_gui)

    def _init_sublayouts(self):
        self._main_layout = QtGui.QVBoxLayout ()
        self._advanced_layout = QtGui.QVBoxLayout ()
        self._create_advanced_cb ()

    def _init_save_buttons(self):
        self._buttons_layout = QtGui.QHBoxLayout ()
        self.loadButton = QtGui.QPushButton ('&Load settings')
        self.saveButton = QtGui.QPushButton ('&Save settings')
        self.defaultsButton = QtGui.QPushButton ('&Defaults')
        self._buttons_layout.addWidget (self.loadButton)
        self._buttons_layout.addWidget (self.saveButton)
        self._buttons_layout.addWidget (self.defaultsButton)
        self.layout.addLayout (self._buttons_layout)
        self.loadButton.clicked.connect (self.openFile)
        self.saveButton.clicked.connect (self.saveFile)
        self.defaultsButton.clicked.connect (self.restoreDefaults)

    def _init_tabs(self):
        self.layout.addWidget (self.tabs)

    def _init_dialog_buttons(self):
        self._dialog_buttons = QtGui.QHBoxLayout ()
        self._okay = QtGui.QPushButton ('Ok')
        self._close = QtGui.QPushButton ('Close')
        self._apply = QtGui.QPushButton ('Apply')
        self._dialog_buttons.addWidget (self._okay)
        self._dialog_buttons.addWidget (self._close)
        self._dialog_buttons.addWidget (self._apply)
        self._okay.clicked.connect (self.accept)
        self._okay.clicked.connect (self.close)
        self._apply.clicked.connect (self.accept)
        self._close.clicked.connect (self.close)
        self.layout.addLayout (self._dialog_buttons)

    def _create_advanced_cb(self):
        self._advanced_cb = QtGui.QCheckBox ('Enable Advanced Settings')
        self._advanced_layout.insertWidget (0, self._advanced_cb)
        self._advanced_cb.toggled.connect (self._toggle_advanced_settings)

    def _toggle_advanced_settings(self):
        if self._advanced_cb.isChecked ():
            self._enable_advanced (True)
        else:
            self._enable_advanced (False)

    def _enable_advanced(self, enable):
        for lst in self.parameter.get_special_para_names ().values ():
            for param in lst:
                box = self.boxes[param]
                if type (box) is not list:
                    box.setEnabled (enable)
                else:
                    for b in box:
                        b.setEnabled (enable)

    def set_tune_mode(self, bool):
        names = ['Directories', 'NLLoc',
                 'Seismic Moment']
        for name in names:
            self.hide_groupbox (name)
        if bool:
            self._apply.hide ()
            self._okay.hide ()
            self._close.hide ()
        else:
            self._apply.show ()
            self._okay.show ()
            self._close.show ()

    def init_boxes(self, parameter_names):
        grid = QtGui.QGridLayout ()

        for index1, name in enumerate (parameter_names):
            default_item = self.parameter.get_defaults ()[name]
            tooltip = default_item['tooltip']
            tooltip += ' | type: {}'.format (default_item['type'])
            if not type (default_item['type']) == tuple:
                typ = default_item['type']
                box = self.create_box (typ, tooltip)
                self.boxes[name] = box
                namestring = default_item['namestring']
            elif type (default_item['type']) == tuple:
                boxes = []
                values = self.parameter[name]
                for index2, val in enumerate (values):
                    typ = default_item['type'][index2]
                    boxes.append (self.create_box (typ, tooltip))
                headline = default_item['namestring'][1:]
                box, lower = self.create_multi_box (boxes, headline)
                self.boxes[name] = boxes
                namestring = default_item['namestring'][0]
            text = namestring + ' [?]'
            label = QtGui.QLabel (text)
            self.labels[name] = label
            label.setToolTip (tooltip)
            grid.addWidget (label, index1, 1)
            grid.addWidget (box, index1, 2)
        return grid

    def create_box(self, typ, tooltip):
        if typ == str:
            box = QtGui.QLineEdit ()
        elif typ == float:
            box = QtGui.QDoubleSpinBox ()
            box.setDecimals (4)
            box.setRange (-10e4, 10e4)
        elif typ == int:
            box = QtGui.QSpinBox ()
        elif typ == bool:
            box = QtGui.QCheckBox ()
        else:
            raise TypeError ('Unrecognized type {}'.format (typ))
        return box

    def create_multi_box(self, boxes, headline=None):
        box = QtGui.QWidget ()
        gl = QtGui.QGridLayout ()
        column = 0
        if headline:
            for index, item in enumerate (headline):
                if not item:
                    continue
                gl.addWidget (QtGui.QLabel (item), index, 0, 2)
                column = 1
        for index, b in enumerate (boxes):
            gl.addWidget (b, index, column)
        box.setLayout (gl)
        return box, column

    def add_tab(self, layout, name):
        widget = QtGui.QWidget ()
        scrollA = QtGui.QScrollArea ()
        scrollA.setWidgetResizable (True)
        scrollA.setWidget (widget)
        widget.setLayout (layout)
        self.tabs.addTab (scrollA, name)

    def add_main_parameters_tab(self):
        self.add_to_layout (self._main_layout, 'Directories',
                            self.parameter.get_main_para_names ()['dirs'], 0)
        self.add_to_layout (self._main_layout, 'NLLoc',
                            self.parameter.get_main_para_names ()['nlloc'], 1)
        self.add_to_layout (self._main_layout, 'Seismic Moment',
                            self.parameter.get_main_para_names ()['smoment'], 2)
        self.add_to_layout (self._main_layout, 'Local Magnitude',
                            self.parameter.get_main_para_names ()['localmag'], 3)
        self.add_to_layout (self._main_layout, 'Filter Settings',
                            self.parameter.get_main_para_names ()['filter'], 4)
        self.add_to_layout (self._main_layout, 'Common Settings Characteristic Function',
                            self.parameter.get_main_para_names ()['pick'], 5)
        self.add_tab (self._main_layout, 'Main Settings')

    def add_special_pick_parameters_tab(self):
        self.add_to_layout (self._advanced_layout, 'Z-component',
                            self.parameter.get_special_para_names ()['z'], 1)
        self.add_to_layout (self._advanced_layout, 'H-components',
                            self.parameter.get_special_para_names ()['h'], 2)
        self.add_to_layout (self._advanced_layout, 'First-motion picker',
                            self.parameter.get_special_para_names ()['fm'], 3)
        self.add_to_layout (self._advanced_layout, 'Quality assessment',
                            self.parameter.get_special_para_names ()['quality'], 4)
        self.add_tab (self._advanced_layout, 'Advanced Settings')

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
        for groupbox in self.groupboxes.values ():
            layout = groupbox._parentLayout
            position = groupbox._position
            layout.insertWidget (position, groupbox)

    def get_groupbox_exclusive(self, name):
        widget = QtGui.QWidget (self, 1)
        layout = QtGui.QVBoxLayout ()
        widget.setLayout (layout)
        layout.addWidget (self.groupboxes[name])
        self._exclusive_widgets.append (widget)
        return widget

    def get_groupbox_dialog(self, name):
        widget = self.get_groupbox_exclusive (name)
        dialog = QtGui.QDialog (self.parent ())
        layout = QtGui.QVBoxLayout ()
        dialog.setLayout (layout)
        buttonbox = QtGui.QDialogButtonBox (QDialogButtonBox.Ok |
                                            QDialogButtonBox.Cancel)
        buttonbox.accepted.connect (dialog.accept)
        buttonbox.accepted.connect (self.refresh)
        buttonbox.accepted.connect (self.params_from_gui)
        buttonbox.rejected.connect (dialog.reject)
        buttonbox.rejected.connect (self.refresh)
        buttonbox.rejected.connect (self.params_to_gui)
        layout.addWidget (widget)
        layout.addWidget (buttonbox)
        self._exclusive_dialog = dialog
        return dialog

    def add_to_layout(self, layout, name, items, position):
        groupbox = QtGui.QGroupBox (name)
        groupbox._position = position
        groupbox._parentLayout = layout
        self.groupboxes[name] = groupbox
        groupbox.setLayout (self.init_boxes (items))
        layout.insertWidget (position, groupbox)

    def show_groupboxes(self):
        for name in self.groupboxes.keys ():
            self.show_groupbox (name)
        self._advanced_cb.show ()

    def hide_groupboxes(self):
        for name in self.groupboxes.keys ():
            self.hide_groupbox (name)
        self._advanced_cb.hide ()

    def show_groupbox(self, name):
        if name in self.groupboxes.keys ():
            self.groupboxes[name].show ()
        else:
            print ('Groupbox {} not part of object.'.format (name))

    def hide_groupbox(self, name):
        if name in self.groupboxes.keys ():
            self.groupboxes[name].hide ()
        else:
            print ('Groupbox {} not part of object.'.format (name))

    def show_file_buttons(self):
        self.saveButton.show ()
        self.loadButton.show ()
        self.defaultsButton.show ()

    def hide_file_buttons(self):
        self.saveButton.hide ()
        self.loadButton.hide ()
        self.defaultsButton.hide ()

    def show_parameter(self, name=None):
        if not name:
            for name in self.boxes.keys ():
                self.show_parameter (name)
            return
        if name in self.boxes.keys () and name in self.labels.keys ():
            # comprising case type(self.boxes[name]) == list
            boxes = self.boxes[name]
            if not type (boxes) == list:
                boxes = [boxes]
            for box in boxes:
                box.show ()
            self.labels[name].show ()
        else:
            print ('Parameter {} not part of object.'.format (name))

    def hide_parameter(self, name=None):
        if not name:
            for name in self.boxes.keys ():
                self.hide_parameter (name)
            return
        if name in self.boxes.keys () and name in self.labels.keys ():
            # comprising case type(self.boxes[name]) == list
            boxes = self.boxes[name]
            if not type (boxes) == list:
                boxes = [boxes]
            for box in boxes:
                box.hide ()
            self.labels[name].hide ()
        else:
            print ('Parameter {} not part of object.'.format (name))

    def params_from_gui(self):
        for param in self.parameter.get_all_para_names ():
            box = self.boxes[param]
            value = self.getValue (box)
            self.parameter.checkValue (param, value)
            self.parameter.setParamKV (param, value)
        return self.parameter

    def params_to_gui(self, tuneMode=False):
        for param in self.parameter.get_all_para_names ():
            if param == 'eventID':
                if tuneMode:
                    continue
            box = self.boxes[param]
            value = self.parameter[param]
            # self.parameter.checkValue(param, value)
            self.setValue (box, value)

    def setValue(self, box, value):
        if type (box) == QtGui.QLineEdit:
            box.setText (str (value))
        elif type (box) == QtGui.QSpinBox or type (box) == QtGui.QDoubleSpinBox:
            if not value:
                value = 0.
            box.setValue (value)
        elif type (box) == QtGui.QCheckBox:
            if value == 'True':
                value = True
            if value == 'False':
                value = False
            box.setChecked (value)
        elif type (box) == list:
            for index, b in enumerate (box):
                self.setValue (b, value[index])

    def getValue(self, box):
        if type (box) == QtGui.QLineEdit:
            value = str (box.text ())
        elif type (box) == QtGui.QSpinBox or type (box) == QtGui.QDoubleSpinBox:
            value = box.value ()
        elif type (box) == QtGui.QCheckBox:
            value = box.isChecked ()
        elif type (box) == list:
            value = []
            for b in box:
                value.append (self.getValue (b))
            value = tuple (value)
        return value

    def openFile(self):
        fd = QtGui.QFileDialog ()
        fname = fd.getOpenFileName (self, 'Browse for settings file.',
                                    filter='PyLoT input file (*.in)')
        if fname[0]:
            try:
                self.parameter.from_file (fname[0])
                self.params_to_gui (tuneMode=True)
            except Exception as e:
                self._warn ('Could not open file {}:\n{}'.format (fname[0], e))
                return

    def saveFile(self):
        fd = QtGui.QFileDialog ()
        fname = fd.getSaveFileName (self, 'Browse for settings file.',
                                    filter='PyLoT input file (*.in)')
        if fname[0]:
            try:
                self.params_from_gui ()
                self.parameter.export2File (fname[0])
            except Exception as e:
                self._warn ('Could not save file {}:\n{}'.format (fname[0], e))
                return

    def restoreDefaults(self):
        try:
            self.parameter.reset_defaults ()
            self.params_to_gui (tuneMode=True)
        except Exception as e:
            self._warn ('Could not restore defaults:\n{}'.format (e))
            return

    def show(self):
        self.refresh ()
        self.show_parameter ()
        if hasattr (self, '_exclusive_dialog'):
            self._exclusive_dialog.close ()
        self._exclusive_widgets = []
        QtGui.QWidget.show (self)

    def close(self):
        self.rejected.emit ('reject')
        QtGui.QWidget.close (self)

    def accept(self):
        self.accepted.emit ('accept')

    def _warn(self, message):
        self.qmb = QtGui.QMessageBox (QtGui.QMessageBox.Icon.Warning,
                                      'Warning', message)
        self.qmb.show ()


class PropertiesDlg (QDialog):
    def __init__(self, parent=None, infile=None, inputs=None):
        super (PropertiesDlg, self).__init__ (parent)

        self.infile = infile
        self.inputs = inputs

        self.setWindowTitle ("PyLoT Properties")
        self.tabWidget = QTabWidget ()
        self.tabWidget.addTab (InputsTab (self), "Inputs")
        # self.tabWidget.addTab(OutputsTab(self), "Outputs")
        self.tabWidget.addTab (PhasesTab (self, inputs), "Phases")
        self.tabWidget.addTab (GraphicsTab (self), "Graphics")
        # self.tabWidget.addTab(LocalisationTab(self), "Loc. Tools")
        self.tabWidget.addTab (LocalisationTab (self), "NonLinLoc")
        self.tabWidget.addTab (ChannelOrderTab (self), "Channel Order")
        self.buttonBox = QDialogButtonBox (QDialogButtonBox.Ok |
                                           QDialogButtonBox.Apply |
                                           QDialogButtonBox.Close |
                                           QDialogButtonBox.RestoreDefaults)

        layout = QVBoxLayout ()
        layout.addWidget (self.tabWidget)
        layout.addWidget (self.buttonBox)
        self.setLayout (layout)
        self.setFixedWidth (700)

        self.buttonBox.accepted.connect (self.accept)
        self.buttonBox.rejected.connect (self.reject)
        self.buttonBox.button (QDialogButtonBox.Apply).clicked.connect (self.apply)
        self.buttonBox.button (QDialogButtonBox.RestoreDefaults).clicked.connect (self.restore)

    def getinfile(self):
        return self.infile

    def accept(self, *args, **kwargs):
        self.apply ()
        QDialog.accept (self)

    def apply(self):
        for widint in range (self.tabWidget.count ()):
            curwid = self.tabWidget.widget (widint)
            values = curwid.getValues ()
            if values is not None:
                self.setValues (values)

    def close(self):
        self.reset_current ()
        QDialog.close (self)

    def show(self):
        self.keep_current ()
        QDialog.show (self)

    def restore(self):
        for widint in range (self.tabWidget.count ()):
            curwid = self.tabWidget.widget (widint)
            values = curwid.resetValues (self.getinfile ())
            if values is not None:
                self.setValues (values)

    def keep_current(self):
        self._current_values = []
        for widint in range (self.tabWidget.count ()):
            curwid = self.tabWidget.widget (widint)
            values = curwid.getValues ()
            if values is not None:
                self._current_values.append (values)

    def reset_current(self):
        for values in self._current_values ():
            self.setValues (values)

    @staticmethod
    def setValues(tabValues):
        settings = QSettings ()
        compclass = settings.value ('compclass')
        if not compclass:
            print ('Warning: No settings for channel components found. Using default')
            compclass = SetChannelComponents ()

        for setting, value in tabValues.items ():
            settings.setValue (setting, value)
            if value is not None:
                if setting.startswith ('Channel Z'):
                    component = 'Z'
                    compclass.setCompPosition (value, component, False)
                elif setting.startswith ('Channel E'):
                    component = 'E'
                    compclass.setCompPosition (value, component, False)
                elif setting.startswith ('Channel N'):
                    component = 'N'
                    compclass.setCompPosition (value, component, False)

        settings.sync ()


class PropTab (QWidget):
    def __init__(self, parent=None):
        super (PropTab, self).__init__ (parent)

    def getValues(self):
        return None

    def resetValues(self, infile=None):
        return None


class InputsTab (PropTab):
    def __init__(self, parent, infile=None):
        super (InputsTab, self).__init__ (parent)

        settings = QSettings ()
        pylot_user = getpass.getuser ()
        fulluser = settings.value ("user/FullName")
        login = settings.value ("user/Login")

        fullNameLabel = QLabel ("Full name for user '{0}': ".format (pylot_user))

        # get the full name of the actual user
        self.fullNameEdit = QLineEdit ()
        try:
            self.fullNameEdit.setText (fulluser)
        except TypeError as e:
            self.fullNameEdit.setText (fulluser[0])

        # information about data structure
        dataroot = settings.value ("data/dataRoot")
        curstructure = settings.value ("data/Structure")
        dataDirLabel = QLabel ("data root directory: ")
        self.dataDirEdit = QLineEdit ()
        self.dataDirEdit.setText (dataroot)
        self.dataDirEdit.selectAll ()
        structureLabel = QLabel ("data structure: ")
        self.structureSelect = QComboBox ()

        from pylot.core.util.structure import DATASTRUCTURE

        self.structureSelect.addItems (list (DATASTRUCTURE.keys ()))

        dsind = findComboBoxIndex (self.structureSelect, curstructure)

        self.structureSelect.setCurrentIndex (dsind)

        layout = QGridLayout ()
        layout.addWidget (dataDirLabel, 0, 0)
        layout.addWidget (self.dataDirEdit, 0, 1)
        layout.addWidget (fullNameLabel, 1, 0)
        layout.addWidget (self.fullNameEdit, 1, 1)
        layout.addWidget (structureLabel, 2, 0)
        layout.addWidget (self.structureSelect, 2, 1)

        self.setLayout (layout)

    def getValues(self):
        values = {"data/dataRoot": self.dataDirEdit.text (),
                  "user/FullName": self.fullNameEdit.text (),
                  "data/Structure": self.structureSelect.currentText ()}
        return values

    def resetValues(self, infile):
        para = PylotParameter (infile)
        datstruct = para.get ('datastructure')
        if datstruct == 'SeisComp':
            index = 0
        else:
            index = 2
        datapath = para.get ('datapath')
        rootpath = para.get ('rootpath')
        database = para.get ('database')
        if isinstance (database, int):
            database = str (database)
        path = os.path.join (os.path.expanduser ('~'), rootpath, datapath, database)
        values = {"data/dataRoot": self.dataDirEdit.setText ("%s" % path),
                  "user/FullName": self.fullNameEdit.text (),
                  "data/Structure": self.structureSelect.setCurrentIndex (index)}
        return values


class OutputsTab (PropTab):
    def __init__(self, parent=None, infile=None):
        super (OutputsTab, self).__init__ (parent)

        settings = QSettings ()
        curval = settings.value ("output/Format", None)

        eventOutputLabel = QLabel ("event/picks output format")
        self.eventOutputComboBox = QComboBox ()
        eventoutputformats = OUTPUTFORMATS.keys ()
        self.eventOutputComboBox.addItems (eventoutputformats)

        ind = findComboBoxIndex (self.eventOutputComboBox, curval)

        self.eventOutputComboBox.setCurrentIndex (ind)
        layout = QGridLayout ()
        layout.addWidget (eventOutputLabel, 0, 0)
        layout.addWidget (self.eventOutputComboBox, 0, 1)

        self.setLayout (layout)

    def getValues(self):
        values = {"output/Format": self.eventOutputComboBox.currentText ()}
        return values

    def resetValues(self, infile):
        values = {"output/Format": self.eventOutputComboBox.setCurrentIndex (1)}
        return values


class PhasesTab (PropTab):
    def __init__(self, parent=None, inputs=None):
        super (PhasesTab, self).__init__ (parent)
        self.inputs = inputs

        self.Pphases = 'P, Pg, Pn, PmP, P1, P2, P3'
        self.Sphases = 'S, Sg, Sn, SmS, S1, S2, S3'

        self.PphasesEdit = QLineEdit ()
        self.SphasesEdit = QLineEdit ()

        self.pickDefaultsButton = QtGui.QPushButton ('Choose default phases...')
        PphasesLabel = QLabel ("P Phases to pick")
        SphasesLabel = QLabel ("S Phases to pick")

        settings = QSettings ()
        Pphases = settings.value ('p_phases')
        Sphases = settings.value ('s_phases')

        self.PphasesEdit.setText ("%s" % Pphases)
        self.SphasesEdit.setText ("%s" % Sphases)

        self.main_layout = QtGui.QHBoxLayout ()
        layout = QGridLayout ()
        layout.addWidget (PphasesLabel, 0, 0)
        layout.addWidget (SphasesLabel, 1, 0)

        layout.addWidget (self.PphasesEdit, 0, 1)
        layout.addWidget (self.SphasesEdit, 1, 1)
        self.main_layout.addLayout (layout)
        self.main_layout.addWidget (self.pickDefaultsButton)
        self.setLayout (self.main_layout)

        self.connectSignals ()

    def connectSignals(self):
        self.pickDefaultsButton.clicked.connect (self.get_defaults)

    def get_defaults(self):
        phases = [p.strip () for p in self.Pphases.split (',')] + [s.strip () for s in self.Sphases.split (',')]
        p_current = [p.strip () for p in self.PphasesEdit.text ().split (',')]
        s_current = [s.strip () for s in self.SphasesEdit.text ().split (',')]

        current_phases = p_current + s_current
        if self.inputs:
            parameter = self.inputs
            if parameter.get ('extent') == 'global':
                # get all default phase names known to obspy.taup
                # in a list and remove duplicates
                phases = list(set(get_phase_names ('ttall')))
                phases.sort()

        phaseDefaults = PhaseDefaults (self, phase_defaults=phases,
                                       current_phases=current_phases)
        if phaseDefaults.exec_ ():
            phase_dict = self.sortPhases (phaseDefaults.selected_phases)
            p_phases = ''
            s_phases = ''
            for index, p in enumerate (phase_dict['P']):
                p_phases += p
                if not index == len (phase_dict['P']) - 1:
                    p_phases += ', '
            for index, s in enumerate (phase_dict['S']):
                s_phases += s
                if not index == len (phase_dict['S']) - 1:
                    s_phases += ', '
            self.PphasesEdit.setText (p_phases)
            self.SphasesEdit.setText (s_phases)

    def sortPhases(self, phases):
        sorted_phases = {'P': [],
                         'S': []}
        for phase in phases:
            idf_phase = loopIdentifyPhase (phase)
            if idf_phase:
                sorted_phases[identifyPhase (idf_phase)].append (phase)
        return sorted_phases

    def getValues(self):
        p_phases = self.PphasesEdit.text ()
        s_phases = self.SphasesEdit.text ()
        values = {'p_phases': p_phases,
                  's_phases': s_phases}
        return values

    def resetValues(self, infile=None):
        values = {'p_phases': self.PphasesEdit.setText (self.Pphases),
                  's_phases': self.SphasesEdit.setText (self.Sphases)}
        return values


class GraphicsTab (PropTab):
    def __init__(self, parent=None):
        super (GraphicsTab, self).__init__ (parent)
        self.init_layout ()
        self.add_pg_cb ()
        self.add_nth_sample ()
        self.setLayout (self.main_layout)

    def init_layout(self):
        self.main_layout = QGridLayout ()

    def add_nth_sample(self):
        settings = QSettings ()
        nth_sample = settings.value ("nth_sample")
        if not nth_sample:
            nth_sample = 1

        self.spinbox_nth_sample = QtGui.QSpinBox ()
        label = QLabel ('nth sample')
        label.setToolTip ('Plot every nth sample (to speed up plotting)')
        self.spinbox_nth_sample.setMinimum (1)
        self.spinbox_nth_sample.setMaximum (10e3)
        self.spinbox_nth_sample.setValue (int (nth_sample))
        self.main_layout.addWidget (label, 1, 0)
        self.main_layout.addWidget (self.spinbox_nth_sample, 1, 1)

    def add_pg_cb(self):
        text = {True: 'Use pyqtgraphic library for plotting',
                False: 'Cannot use library: pyqtgraphic not found on system'}
        label = QLabel ('PyQt graphic')
        label.setToolTip (text[bool (pg)])
        label.setEnabled (bool (pg))
        self.checkbox_pg = QtGui.QCheckBox ()
        self.checkbox_pg.setEnabled (bool (pg))
        self.checkbox_pg.setChecked (bool (pg))
        self.main_layout.addWidget (label, 0, 0)
        self.main_layout.addWidget (self.checkbox_pg, 0, 1)

    def getValues(self):
        values = {'nth_sample': self.spinbox_nth_sample.value (),
                  'pyqtgraphic': self.checkbox_pg.isChecked ()}
        return values

    def resetValues(self, infile=None):
        values = {'nth_sample': self.spinbox_nth_sample.setValue (1),
                  'pyqtgraphic': self.checkbox_pg.setChecked (True)}
        return values


class ChannelOrderTab (PropTab):
    def __init__(self, parent=None, infile=None):
        super (ChannelOrderTab, self).__init__ (parent)

        settings = QSettings ()
        compclass = settings.value ('compclass')
        if not compclass:
            print ('Warning: No settings for channel components found. Using default')
            compclass = SetChannelComponents ()

        ChannelOrderLabelZ = QLabel ("Channel Z [up/down, default=3]")
        ChannelOrderLabelN = QLabel ("Channel N [north/south, default=1]")
        ChannelOrderLabelE = QLabel ("Channel E [east/west, default=2]")
        self.ChannelOrderZEdit = QLineEdit ()
        self.ChannelOrderZEdit.setMaxLength (1)
        self.ChannelOrderZEdit.setFixedSize (20, 20)
        self.ChannelOrderNEdit = QLineEdit ()
        self.ChannelOrderNEdit.setMaxLength (1)
        self.ChannelOrderNEdit.setFixedSize (20, 20)
        self.ChannelOrderEEdit = QLineEdit ()
        self.ChannelOrderEEdit.setMaxLength (1)
        self.ChannelOrderEEdit.setFixedSize (20, 20)
        # get channel order settings
        zcomp = compclass.getCompPosition ('Z')
        ncomp = compclass.getCompPosition ('N')
        ecomp = compclass.getCompPosition ('E')
        self.ChannelOrderZEdit.setText ("%s" % zcomp)
        self.ChannelOrderNEdit.setText ("%s" % ncomp)
        self.ChannelOrderEEdit.setText ("%s" % ecomp)

        layout = QGridLayout ()
        layout.addWidget (ChannelOrderLabelZ, 0, 0)
        layout.addWidget (ChannelOrderLabelN, 1, 0)
        layout.addWidget (ChannelOrderLabelE, 2, 0)
        layout.addWidget (self.ChannelOrderZEdit, 0, 1)
        layout.addWidget (self.ChannelOrderNEdit, 1, 1)
        layout.addWidget (self.ChannelOrderEEdit, 2, 1)

        self.setLayout (layout)
        self.connectSignals ()

    def connectSignals(self):
        self.ChannelOrderZEdit.textEdited.connect (self.checkDoubleZ)
        self.ChannelOrderNEdit.textEdited.connect (self.checkDoubleN)
        self.ChannelOrderEEdit.textEdited.connect (self.checkDoubleE)

    def checkDoubleZ(self, text):
        self.checkDouble (text, 'Z')

    def checkDoubleN(self, text):
        self.checkDouble (text, 'N')

    def checkDoubleE(self, text):
        self.checkDouble (text, 'E')

    def checkDouble(self, text, comp):
        channelOrderEdits = {
            'Z': self.ChannelOrderZEdit,
            'N': self.ChannelOrderNEdit,
            'E': self.ChannelOrderEEdit
        }
        for key in channelOrderEdits.keys ():
            if key == comp:
                continue
            if str (channelOrderEdits[key].text ()) == str (text):
                channelOrderEdits[key].setText ('')

    def getValues(self):
        values = {"Channel Z [up/down, default=3]": int (self.ChannelOrderZEdit.text ()),
                  "Channel N [north/south, default=1]": int (self.ChannelOrderNEdit.text ()),
                  "Channel E [east/west, default=2]": int (self.ChannelOrderEEdit.text ())}
        return values

    def resetValues(self, infile=None):
        Zdefault = 3
        Ndefault = 1
        Edefault = 2
        values = {"Channel Z [up/down, default=3]": self.ChannelOrderZEdit.setText ("%d" % Zdefault),
                  "Channel N [north/south, default=1]": self.ChannelOrderNEdit.setText ("%d" % Ndefault),
                  "Channel E [east/west, default=2]": self.ChannelOrderEEdit.setText ("%d" % Edefault)}
        return values

        # MP MP: No idea why this function exists!?
        # def getComponents(self):
        #     self.CompName = dict(Z='10', N='11', E='12')


class LocalisationTab (PropTab):
    def __init__(self, parent=None, infile=None):
        super (LocalisationTab, self).__init__ (parent)

        settings = QSettings ()
        curtool = settings.value ("loc/tool", None)

        # loctoollabel = QLabel("location tool")
        self.locToolComboBox = QComboBox ()
        # loctools = LOCTOOLS.keys()
        # self.locToolComboBox.addItems(loctools)

        # toolind = findComboBoxIndex(self.locToolComboBox, curtool)

        # self.locToolComboBox.setCurrentIndex(toolind)

        curroot = settings.value ("{0}/rootPath".format (curtool), None)
        curbin = settings.value ("{0}/binPath".format (curtool), None)

        self.rootlabel = QLabel ("root directory")
        self.binlabel = QLabel ("bin directory")

        self.rootedit = QLineEdit ('')
        self.binedit = QLineEdit ('')

        if curroot is not None:
            self.rootedit.setText (curroot)
        if curbin is not None:
            self.binedit.setText (curbin)

        rootBrowse = QPushButton ('...', self)
        rootBrowse.clicked.connect (lambda: self.selectDirectory (self.rootedit))

        binBrowse = QPushButton ('...', self)
        binBrowse.clicked.connect (lambda: self.selectDirectory (self.binedit))

        # self.locToolComboBox.currentIndexChanged.connect(self.updateUi)

        self.updateUi ()

        layout = QGridLayout ()
        # layout.addWidget(loctoollabel, 0, 0)
        # layout.addWidget(self.locToolComboBox, 0, 1)
        layout.addWidget (self.rootlabel, 1, 0)
        layout.addWidget (self.rootedit, 1, 1)
        layout.addWidget (rootBrowse, 1, 2)
        layout.addWidget (self.binlabel, 2, 0)
        layout.addWidget (self.binedit, 2, 1)
        layout.addWidget (binBrowse, 2, 2)

        self.setLayout (layout)

    def updateUi(self):
        curtool = self.locToolComboBox.currentText ()
        # if curtool is not None:
        self.rootlabel.setText ("{0} root directory".format (curtool))
        self.binlabel.setText ("{0} bin directory".format (curtool))

    def selectDirectory(self, edit):
        selected_directory = QFileDialog.getExistingDirectory ()
        # check if string is empty
        if selected_directory:
            edit.setText (selected_directory)

    def getValues(self):
        loctool = self.locToolComboBox.currentText ()
        values = {"{0}/rootPath".format (loctool): self.rootedit.text (),
                  "{0}/binPath".format (loctool): self.binedit.text ()}
        # "loc/tool": loctool}
        return values

    def resetValues(self, infile):
        para = PylotParameter (infile)
        nllocroot = para.get ('nllocroot')
        nllocbin = para.get ('nllocbin')
        loctool = self.locToolComboBox.setCurrentIndex (3)
        values = {"nll/rootPath": self.rootedit.setText ("%s" % nllocroot),
                  "nll/binPath": self.binedit.setText ("%s" % nllocbin)}


class NewEventDlg (QDialog):
    def __init__(self, parent=None, titleString="Create a new event"):
        """
        QDialog object utilized to create a new event manually.
        """
        super (NewEventDlg, self).__init__ ()

        self.setupUI ()

        now = datetime.datetime.now ()
        self.eventTimeEdit.setDateTime (now)
        # event dates in the future are forbidden
        self.eventTimeEdit.setMaximumDateTime (now)

        self.latEdit.setText ("51.0000")
        self.lonEdit.setText ("7.0000")
        self.depEdit.setText ("10.0")

        self.buttonBox.accepted.connect (self.accept)
        self.buttonBox.rejected.connect (self.reject)

    def getValues(self):
        return {'origintime': self.eventTimeEdit.dateTime ().toPython (),
                'latitude': self.latEdit.text (),
                'longitude': self.lonEdit.text (),
                'depth': self.depEdit.text ()}

    def setupUI(self):
        # create widget objects
        timeLabel = QLabel ()
        timeLabel.setText ("Select time: ")
        self.eventTimeEdit = QDateTimeEdit ()
        latLabel = QLabel ()
        latLabel.setText ("Latitude: ")
        self.latEdit = QLineEdit ()
        lonLabel = QLabel ()
        lonLabel.setText ("Longitude: ")
        self.lonEdit = QLineEdit ()
        depLabel = QLabel ()
        depLabel.setText ("Depth: ")
        self.depEdit = QLineEdit ()

        self.buttonBox = QDialogButtonBox (QDialogButtonBox.Ok |
                                           QDialogButtonBox.Cancel)

        grid = QGridLayout ()
        grid.addWidget (timeLabel, 0, 0)
        grid.addWidget (self.eventTimeEdit, 0, 1)
        grid.addWidget (latLabel, 1, 0)
        grid.addWidget (self.latEdit, 1, 1)
        grid.addWidget (lonLabel, 2, 0)
        grid.addWidget (self.lonEdit, 2, 1)
        grid.addWidget (depLabel, 3, 0)
        grid.addWidget (self.depEdit, 3, 1)
        grid.addWidget (self.buttonBox, 4, 1)

        self.setLayout (grid)


class FilterOptionsDialog (QDialog):
    def __init__(self, parent=None, titleString="Filter options",
                 filterOptions=None):
        """
        PyLoT widget FilterOptionsDialog is a QDialog object. It is an UI to
        adjust parameters for filtering seismic data.
        """
        super (FilterOptionsDialog, self).__init__ ()
        if parent is not None and parent.getFilters ():
            self.filterOptions = parent.getFilters ()
        # elif filterOptions is not None:
        #     self.filterOptions = filterOptions
        else:
            self.filterOptions = {'P': FilterOptions (),
                                  'S': FilterOptions ()}

        self.setWindowTitle (titleString)
        self.filterOptionWidgets = {'P': FilterOptionsWidget (self.filterOptions['P']),
                                    'S': FilterOptionsWidget (self.filterOptions['S'])}
        self.setupUi ()
        self.updateUi ()
        self.connectButtons ()

    def setupUi(self):
        self.main_layout = QtGui.QVBoxLayout ()
        self.filter_layout = QtGui.QHBoxLayout ()
        self.groupBoxes = {'P': QtGui.QGroupBox ('P Filter'),
                           'S': QtGui.QGroupBox ('S Filter')}

        self.buttonBox = QDialogButtonBox (QDialogButtonBox.Ok |
                                           QDialogButtonBox.Cancel)

        for key in ['P', 'S']:
            groupbox = self.groupBoxes[key]
            box_layout = QtGui.QVBoxLayout ()
            groupbox.setLayout (box_layout)

            self.filter_layout.addWidget (groupbox)
            box_layout.addWidget (self.filterOptionWidgets[key])

        self.main_layout.addLayout (self.filter_layout)
        self.main_layout.addWidget (self.buttonBox)
        self.setLayout (self.main_layout)

    def connectButtons(self):
        self.buttonBox.accepted.connect (self.accept)
        self.buttonBox.rejected.connect (self.reject)

    def accept(self):
        self.updateUi ()
        QDialog.accept (self)

    def updateUi(self):
        returnvals = []
        for foWidget in self.filterOptionWidgets.values ():
            foWidget.updateUi ()

    def getFilterOptions(self):
        filteroptions = {'P': self.filterOptionWidgets['P'].getFilterOptions (),
                         'S': self.filterOptionWidgets['S'].getFilterOptions ()}
        return filteroptions


class FilterOptionsWidget (QWidget):
    def __init__(self, filterOptions):
        super (FilterOptionsWidget, self).__init__ ()
        self.filterOptions = filterOptions

        _enable = True
        if self.getFilterOptions ().getFilterType () is None:
            _enable = False

        self.freqminLabel = QLabel ()
        self.freqminLabel.setText ("minimum:")
        self.freqminSpinBox = QDoubleSpinBox ()
        self.freqminSpinBox.setRange (5e-7, 1e6)
        self.freqminSpinBox.setDecimals (2)
        self.freqminSpinBox.setSingleStep (0.01)
        self.freqminSpinBox.setSuffix (' Hz')
        self.freqminSpinBox.setEnabled (_enable)

        self.freqmaxLabel = QLabel ()
        self.freqmaxLabel.setText ("maximum:")
        self.freqmaxSpinBox = QDoubleSpinBox ()
        self.freqmaxSpinBox.setRange (5e-7, 1e6)
        self.freqmaxSpinBox.setDecimals (2)
        self.freqmaxSpinBox.setSingleStep (0.01)
        self.freqmaxSpinBox.setSuffix (' Hz')

        # if _enable:
        #     self.freqminSpinBox.setValue(self.getFilterOptions().getFreq()[0])
        #     if self.getFilterOptions().getFilterType() in ['bandpass',
        #                                                    'bandstop']:
        #         self.freqmaxSpinBox.setValue(
        #             self.getFilterOptions().getFreq()[1])
        # else:
        try:
            self.freqminSpinBox.setValue (self.getFilterOptions ().getFreq ()[0])
            self.freqmaxSpinBox.setValue (self.getFilterOptions ().getFreq ()[1])
        except TypeError as e:
            print (e)
            self.freqmaxSpinBox.setValue (1.)
            self.freqminSpinBox.setValue (.1)

        typeOptions = [None, "bandpass", "bandstop", "lowpass", "highpass"]

        self.orderLabel = QLabel ()
        self.orderLabel.setText ("Order:")
        self.orderSpinBox = QSpinBox ()
        self.orderSpinBox.setRange (2, 10)
        self.orderSpinBox.setEnabled (_enable)
        self.selectTypeLabel = QLabel ()
        self.selectTypeLabel.setText ("Select filter type:")
        self.selectTypeCombo = QComboBox ()
        self.selectTypeCombo.addItems (typeOptions)
        self.selectTypeCombo.setCurrentIndex (typeOptions.index (self.getFilterOptions ().getFilterType ()))
        self.selectTypeLayout = QVBoxLayout ()
        self.selectTypeLayout.addWidget (self.orderLabel)
        self.selectTypeLayout.addWidget (self.orderSpinBox)
        self.selectTypeLayout.addWidget (self.selectTypeLabel)
        self.selectTypeLayout.addWidget (self.selectTypeCombo)

        self.freqGroupBox = QGroupBox ("Frequency range")
        self.freqGroupLayout = QGridLayout ()
        self.freqGroupLayout.addWidget (self.freqminLabel, 0, 0)
        self.freqGroupLayout.addWidget (self.freqminSpinBox, 0, 1)
        self.freqGroupLayout.addWidget (self.freqmaxLabel, 1, 0)
        self.freqGroupLayout.addWidget (self.freqmaxSpinBox, 1, 1)
        self.freqGroupBox.setLayout (self.freqGroupLayout)

        self.freqmaxSpinBox.setEnabled (_enable)

        try:
            self.orderSpinBox.setValue (self.getFilterOptions ().getOrder ())
        except:
            self.orderSpinBox.setValue (2)

        grid = QGridLayout ()
        grid.addWidget (self.freqGroupBox, 0, 2, 1, 2)
        grid.addLayout (self.selectTypeLayout, 1, 2, 1, 2)

        self.setLayout (grid)

        self.freqminSpinBox.valueChanged.connect (self.checkMin)
        self.freqmaxSpinBox.valueChanged.connect (self.checkMax)
        self.orderSpinBox.valueChanged.connect (self.updateUi)
        self.selectTypeCombo.currentIndexChanged.connect (self.updateUi)

    def checkMin(self):
        if not self.freqminSpinBox.value () <= self.freqmaxSpinBox.value ():
            self.freqmaxSpinBox.setValue (self.freqminSpinBox.value ())

    def checkMax(self):
        if not self.freqminSpinBox.value () <= self.freqmaxSpinBox.value ():
            self.freqminSpinBox.setValue (self.freqmaxSpinBox.value ())

    def updateUi(self):
        type = self.selectTypeCombo.currentText ()
        _enable = type in ['bandpass', 'bandstop']
        freq = [self.freqminSpinBox.value (), self.freqmaxSpinBox.value ()]
        self.freqmaxLabel.setEnabled (True)
        self.freqmaxSpinBox.setEnabled (True)
        self.freqminLabel.setEnabled (True)
        self.freqminSpinBox.setEnabled (True)
        self.freqminLabel.setText ("minimum:")
        self.freqmaxLabel.setText ("maximum:")

        if not _enable:
            if type == 'highpass':
                self.freqminLabel.setText ("cutoff:")
                self.freqmaxLabel.setEnabled (False)
                self.freqmaxSpinBox.setEnabled (False)
            elif type == 'lowpass':
                self.freqmaxLabel.setText ("cutoff:")
                self.freqminLabel.setEnabled (False)
                self.freqminSpinBox.setEnabled (False)
        else:
            if not isSorted (freq):
                QMessageBox.warning (self, "Value error",
                                     "Maximum frequency must be at least the "
                                     "same value as minimum frequency (notch)! "
                                     "Adjusted maximum frequency automatically!")
                freq[1] = freq[0]
                self.freqmaxSpinBox.setValue (freq[1])
                self.freqmaxSpinBox.selectAll ()
                self.freqmaxSpinBox.setFocus ()

        self.getFilterOptions ().setFilterType (type)
        self.getFilterOptions ().setFreq (freq)
        self.getFilterOptions ().setOrder (self.orderSpinBox.value ())

    def getFilterOptions(self):
        return self.filterOptions

    @staticmethod
    def getFilterObject():
        dlg = FilterOptionsDialog ()
        if dlg.exec_ ():
            return dlg.getFilterOptions ()
        return None


class LoadDataDlg (QDialog):
    def __init__(self, parent=None):
        super (LoadDataDlg, self).__init__ (parent)

        pass


class HelpForm (QDialog):
    def __init__(self, page=QUrl ('https://ariadne.geophysik.rub.de/trac/PyLoT'),
                 parent=None):
        super (HelpForm, self).__init__ (parent)
        self.setAttribute (Qt.WA_DeleteOnClose)
        self.setAttribute (Qt.WA_GroupLeader)

        backAction = QAction (QIcon (":/back.png"), "&Back", self)
        backAction.setShortcut (QKeySequence.Back)
        homeAction = QAction (QIcon (":/home.png"), "&Home", self)
        homeAction.setShortcut ("Home")
        self.pageLabel = QLabel ()

        toolBar = QToolBar ()
        toolBar.addAction (backAction)
        toolBar.addAction (homeAction)
        toolBar.addWidget (self.pageLabel)
        self.webBrowser = QWebView ()
        self.webBrowser.load (page)

        layout = QVBoxLayout ()
        layout.addWidget (toolBar)
        layout.addWidget (self.webBrowser, 1)
        self.setLayout (layout)

        self.connect (backAction, Signal ("triggered()"),
                      self.webBrowser, Slot ("backward()"))
        self.connect (homeAction, Signal ("triggered()"),
                      self.webBrowser, Slot ("home()"))
        self.connect (self.webBrowser, Signal ("sourceChanged(QUrl)"),
                      self.updatePageTitle)

        self.resize (400, 600)
        self.setWindowTitle ("{0} Help".format (QApplication.applicationName ()))

    def updatePageTitle(self):
        self.pageLabel.setText (self.webBrowser.documentTitle ())


if __name__ == '__main__':
    import doctest

    doctest.testmod ()
