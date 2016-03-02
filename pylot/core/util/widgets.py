# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:27:35 2014

@author: sebastianw
"""

import datetime
import numpy as np

from matplotlib.figure import Figure
try:
    from matplotlib.backends.backend_qt4agg import FigureCanvas
except ImportError:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.widgets import MultiCursor
from PySide.QtGui import QAction, QApplication, QComboBox, QDateTimeEdit, \
    QDialog, QDialogButtonBox, QDoubleSpinBox, QGroupBox, QGridLayout, \
    QIcon, QKeySequence, QLabel, QLineEdit, QMessageBox, QPixmap, QSpinBox, \
    QTabWidget, QToolBar, QVBoxLayout, QWidget, QPushButton, QFileDialog
from PySide.QtCore import QSettings, Qt, QUrl, Signal, Slot
from PySide.QtWebKit import QWebView
from obspy import Stream, UTCDateTime
from pylot.core.read.inputs import FilterOptions
from pylot.core.pick.utils import getSNR, earllatepicker, getnoisewin,\
    getResolutionWindow
from pylot.core.util.defaults import OUTPUTFORMATS, FILTERDEFAULTS, LOCTOOLS
from pylot.core.util.utils import prepTimeAxis, getGlobalTimes, scaleWFData, \
    demeanTrace, isSorted, findComboBoxIndex


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
        self.figure.set_facecolor((.92, .92, .92))
        # attribute plotdict is an dictionary connecting position and a name
        self.plotdict = dict()
        # create axes
        self.axes = self.figure.add_subplot(111)
        # clear axes each time plot is called
        self.axes.hold(True)
        # initialize super class
        FigureCanvas.__init__(self, self.figure)
        # add an cursor for station selection
        self.multiCursor = MultiCursor(self.figure.canvas, (self.axes,),
                                       horizOn=True,
                                       color='m', lw=1)
        # update labels of the entire widget
        self.updateWidget(xlabel, ylabel, title)

    def getPlotDict(self):
        return self.plotdict

    def setPlotDict(self, key, value):
        self.plotdict[key] = value

    def clearPlotDict(self):
        self.plotdict = dict()

    def getParent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    def plotWFData(self, wfdata, title=None, zoomx=None, zoomy=None,
                   noiselevel=None, scaleddata=False):
        self.getAxes().cla()
        self.clearPlotDict()
        wfstart, wfend = getGlobalTimes(wfdata)
        for n, trace in enumerate(wfdata):
            channel = trace.stats.channel
            station = trace.stats.station
            msg = 'plotting %s channel of station %s' % (channel, station)
            print(msg)
            stime = trace.stats.starttime - wfstart
            time_ax = prepTimeAxis(stime, trace)
            if not scaleddata:
                trace.normalize(np.max(np.abs(trace.data)) * 2)
            self.getAxes().plot(time_ax, trace.data + n, 'k')
            if noiselevel is not None:
                for level in noiselevel:
                    self.getAxes().plot([time_ax[0], time_ax[-1]],
                                        [level, level], '--k')
            self.setPlotDict(n, (station, channel))
        xlabel = 'seconds since {0}'.format(wfstart)
        ylabel = ''
        self.updateWidget(xlabel, ylabel, title)
        self.setXLims([0, wfend - wfstart])
        self.setYLims([-0.5, n + 0.5])
        if zoomx is not None:
            self.setXLims(zoomx)
        if zoomy is not None:
            self.setYLims(zoomy)
        self.draw()

    def getAxes(self):
        return self.axes

    def getXLims(self):
        return self.getAxes().get_xlim()

    def getYLims(self):
        return self.getAxes().get_ylim()

    def setXLims(self, lims):
        self.getAxes().set_xlim(lims)

    def setYLims(self, lims):
        self.getAxes().set_ylim(lims)

    def setYTickLabels(self, pos, labels):
        self.getAxes().set_yticks(pos)
        self.getAxes().set_yticklabels(labels)
        self.draw()

    def updateXLabel(self, text):
        self.getAxes().set_xlabel(text)
        self.draw()

    def updateYLabel(self, text):
        self.getAxes().set_ylabel(text)
        self.draw()

    def updateTitle(self, text):
        self.getAxes().set_title(text)
        self.draw()

    def updateWidget(self, xlabel, ylabel, title):
        self.updateXLabel(xlabel)
        self.updateYLabel(ylabel)
        self.updateTitle(title)

    def insertLabel(self, pos, text):
        pos = pos / max(self.getAxes().ylim)
        axann = self.getAxes().annotate(text, xy=(.03, pos),
                                   xycoords='axes fraction')
        axann.set_bbox(dict(facecolor='lightgrey', alpha=.6))

class PickDlg(QDialog):
    def __init__(self, parent=None, data=None, station=None, picks=None,
                 rotate=False):
        super(PickDlg, self).__init__(parent)

        # initialize attributes
        self.station = station
        self.rotate = rotate
        self.components = 'ZNE'
        if picks:
            self.picks = picks
        else:
            self.picks = {}
        self.filteroptions = FILTERDEFAULTS

        # initialize panning attributes
        self.press = None
        self.xpress = None
        self.ypress = None
        self.cur_xlim = None
        self.cur_ylim = None

        # set attribute holding data
        if data is None:
            try:
                data = parent.getData().getWFData().copy()
                self.data = data.select(station=station)
            except AttributeError as e:
                errmsg = 'You either have to put in a data or an appropriate ' \
                         'parent (PyLoT MainWindow) object: {0}'.format(e)
                raise Exception(errmsg)
        else:
            self.data = data

        self.stime, self.etime = getGlobalTimes(self.getWFData())

        # initialize plotting widget
        self.multicompfig = MPLWidget(self)

        # setup ui
        self.setupUi()

        # plot data
        self.getPlotWidget().plotWFData(wfdata=self.getWFData(),
                                        title=self.getStation())

        xlims = self.getPlotWidget().getXLims()
        ylims = self.getPlotWidget().getYLims()

        self.limits = {'x': xlims,
                       'y': ylims}

        self.updateCurrentLimits()

        # set plot labels
        self.setPlotLabels()

        # draw picks if present
        self.drawPicks()

        # connect button press event to an action
        self.cidpress = self.connectPressEvent(self.panPress)
        self.cidmotion = self.connectMotionEvent(self.panMotion)
        self.cidrelease = self.connectReleaseEvent(self.panRelease)
        self.cidscroll = self.connectScrollEvent(self.scrollZoom)

    def setupUi(self):

        # create matplotlib toolbar to inherit functionality
        self.figToolBar = NavigationToolbar2QT(self.getPlotWidget(), self)
        self.figToolBar.hide()

        # create icons
        filter_icon = QIcon()
        filter_icon.addPixmap(QPixmap(':/icons/filter.png'))
        zoom_icon = QIcon()
        zoom_icon.addPixmap(QPixmap(':/icons/zoom_in.png'))
        home_icon = QIcon()
        home_icon.addPixmap(QPixmap(':/icons/zoom_0.png'))
        del_icon = QIcon()
        del_icon.addPixmap(QPixmap(':/icons/delete.png'))

        # create actions
        self.filterAction = createAction(parent=self, text='Filter',
                                         slot=self.filterWFData,
                                         icon=filter_icon,
                                         tip='Toggle filtered/original'
                                             ' waveforms',
                                         checkable=True)
        self.zoomAction = createAction(parent=self, text='Zoom',
                                       slot=self.zoom, icon=zoom_icon,
                                       tip='Zoom into waveform',
                                       checkable=True)
        self.resetZoomAction = createAction(parent=self, text='Home',
                                        slot=self.resetZoom, icon=home_icon,
                                        tip='Reset zoom to original limits')
        self.resetPicksAction = createAction(parent=self, text='Delete Picks',
                                             slot=self.delPicks, icon=del_icon,
                                             tip='Delete current picks.')

        # create other widget elements
        self.selectPhase = QComboBox()
        phaseitems = [None] + FILTERDEFAULTS.keys()
        self.selectPhase.addItems(phaseitems)

        # layout the outermost appearance of the Pick Dialog
        _outerlayout = QVBoxLayout()
        _dialtoolbar = QToolBar()

        # fill toolbar with content
        _dialtoolbar.addAction(self.filterAction)
        _dialtoolbar.addWidget(self.selectPhase)
        _dialtoolbar.addAction(self.zoomAction)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addAction(self.resetZoomAction)
        _dialtoolbar.addSeparator()
        _dialtoolbar.addAction(self.resetPicksAction)

        # layout the innermost widget
        _innerlayout = QVBoxLayout()
        _innerlayout.addWidget(self.multicompfig)

        # add button box to the dialog
        _buttonbox = QDialogButtonBox(QDialogButtonBox.Ok |
                                      QDialogButtonBox.Cancel)

        # merge widgets and layouts to establish the dialog
        _innerlayout.addWidget(_buttonbox)
        _outerlayout.addWidget(_dialtoolbar)
        _outerlayout.addLayout(_innerlayout)

        # connect widget element signals with slots (methods to the dialog
        # object
        self.selectPhase.currentIndexChanged.connect(self.verifyPhaseSelection)
        _buttonbox.accepted.connect(self.accept)
        _buttonbox.rejected.connect(self.reject)

        # finally layout the entire dialog
        self.setLayout(_outerlayout)

    def disconnectPressEvent(self):
        widget = self.getPlotWidget()
        widget.mpl_disconnect(self.cidpress)
        self.cidpress = None

    def connectPressEvent(self, slot):
        widget = self.getPlotWidget()
        return widget.mpl_connect('button_press_event', slot)

    def disconnectScrollEvent(self):
        widget = self.getPlotWidget()
        widget.mpl_disconnect(self.cidscroll)
        self.cidscroll = None

    def connectScrollEvent(self, slot):
        widget = self.getPlotWidget()
        return widget.mpl_connect('scroll_event', slot)

    def disconnectMotionEvent(self):
        widget = self.getPlotWidget()
        widget.mpl_disconnect(self.cidmotion)
        self.cidmotion = None

    def connectMotionEvent(self, slot):
        widget = self.getPlotWidget()
        return widget.mpl_connect('motion_notify_event', slot)

    def disconnectReleaseEvent(self):
        widget = self.getPlotWidget()
        widget.mpl_disconnect(self.cidrelease)
        self.cidrelease = None

    def connectReleaseEvent(self, slot):
        widget = self.getPlotWidget()
        return widget.mpl_connect('button_release_event', slot)

    def verifyPhaseSelection(self):
        phase = self.selectPhase.currentText()
        self.updateCurrentLimits()
        if phase:
            if self.zoomAction.isChecked():
                self.zoomAction.trigger()
            self.disconnectReleaseEvent()
            self.disconnectScrollEvent()
            self.disconnectMotionEvent()
            self.disconnectPressEvent()
            self.cidpress = self.connectPressEvent(self.setIniPick)
            self.filterWFData()
        else:
            self.disconnectPressEvent()
            self.cidpress = self.connectPressEvent(self.panPress)
            self.cidmotion = self.connectMotionEvent(self.panMotion)
            self.cidrelease = self.connectReleaseEvent(self.panRelease)
            self.cidscroll = self.connectScrollEvent(self.scrollZoom)

    def getStartTime(self):
        return self.stime

    def getEndTime(self):
        return self.etime

    def getComponents(self):
        return self.components

    def getStation(self):
        return self.station

    def getPlotWidget(self):
        return self.multicompfig

    def getChannelID(self, key):
        return self.getPlotWidget().getPlotDict()[int(key)][1]

    def getTraceID(self, channels):
        plotDict = self.getPlotWidget().getPlotDict()
        traceIDs = []
        for channel in channels:
            channel = channel.upper()
            for traceID, channelID in plotDict.items():
                if channelID[1].upper().endswith(channel):
                    traceIDs.append(traceID)
        return traceIDs

    def getFilterOptions(self, phase):
        options = self.filteroptions[phase]
        return FilterOptions(**options)

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
        self.setXLims(self.getPlotWidget().getXLims())
        self.setYLims(self.getPlotWidget().getYLims())

    def getWFData(self):
        return self.data

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
        elif component == 'Z':
            wfdata = self.getWFData().select(component=component)
        return wfdata

    def getPicks(self):
        return self.picks

    def resetPicks(self):
        self.picks = {}

    def delPicks(self):
        self.resetPicks()
        self.resetPlot()

    def setIniPick(self, gui_event):

        trace_number = round(gui_event.ydata)

        channel = self.getChannelID(trace_number)
        wfdata = self.selectWFData(channel)

        self.disconnectScrollEvent()
        self.disconnectPressEvent()
        self.disconnectReleaseEvent()
        self.disconnectMotionEvent()
        self.cidpress = self.connectPressEvent(self.setPick)

        if self.selectPhase.currentText().upper().startswith('P'):
            self.setIniPickP(gui_event, wfdata, trace_number)
        elif self.selectPhase.currentText().upper().startswith('S'):
            self.setIniPickS(gui_event, wfdata)

        self.zoomAction.setEnabled(False)

        # reset labels
        self.setPlotLabels()
        self.draw()

    def setIniPickP(self, gui_event, wfdata, trace_number):

        ini_pick = gui_event.xdata

        settings = QSettings()

        nfac = settings.value('picking/nfac_P', 1.3)
        noise_win = settings.value('picking/noise_win_P', 5.)
        gap_win = settings.value('picking/gap_win_P', .2)
        signal_win = settings.value('picking/signal_win_P', 3.)
        itrace = int(trace_number)

        while itrace > len(wfdata) - 1:
            itrace -= 1

        result = getSNR(wfdata, (noise_win, gap_win, signal_win), ini_pick, itrace)

        snr = result[0]
        noiselevel = result[2] * nfac

        x_res = getResolutionWindow(snr)

        # remove mean noise level from waveforms
        wfdata = self.getWFData().copy()
        for trace in wfdata:
            t = prepTimeAxis(trace.stats.starttime - self.getStartTime(), trace)
            inoise = getnoisewin(t, ini_pick, noise_win, gap_win)
            trace = demeanTrace(trace=trace, window=inoise)

        self.setXLims([ini_pick - x_res, ini_pick + x_res])
        self.setYLims(np.array([-noiselevel * 2.5, noiselevel * 2.5]) +
                      trace_number)
        self.getPlotWidget().plotWFData(wfdata=wfdata,
                                        title=self.getStation() +
                                              ' picking mode',
                                        zoomx=self.getXLims(),
                                        zoomy=self.getYLims(),
                                        noiselevel=(trace_number + noiselevel,
                                                    trace_number - noiselevel))

    def setIniPickS(self, gui_event, wfdata):

        ini_pick = gui_event.xdata

        settings = QSettings()

        nfac = settings.value('picking/nfac_P', 1.5)
        noise_win = settings.value('picking/noise_win_P', 5.)
        gap_win = settings.value('picking/gap_win_P', .2)
        signal_win = settings.value('picking/signal_win_P', 3.)

        result = getSNR(wfdata, (noise_win, gap_win, signal_win), ini_pick)

        snr = result[0]
        noiselevel = result[2] * nfac

        data = self.getWFData().copy()

        phase = self.selectPhase.currentText()
        filteroptions = self.getFilterOptions(phase).parseFilterOptions()
        data.filter(**filteroptions)

        for trace in data:
            t = prepTimeAxis(trace.stats.starttime - self.getStartTime(), trace)
            inoise = getnoisewin(t, ini_pick, noise_win, gap_win)
            trace = demeanTrace(trace, inoise)

        horiz_comp = ('n', 'e')

        data = scaleWFData(data, noiselevel * 2.5, horiz_comp)

        x_res = getResolutionWindow(snr)

        self.setXLims(tuple([ini_pick - x_res, ini_pick + x_res]))
        traces = self.getTraceID(horiz_comp)
        traces.sort()
        self.setYLims(tuple(np.array([-0.5, +0.5]) +
                      np.array(traces)))
        noiselevels = [trace + 1 / (2.5 * 2) for trace in traces] +\
                      [trace - 1 / (2.5 * 2) for trace in traces]

        self.getPlotWidget().plotWFData(wfdata=data,
                                        title=self.getStation() +
                                              ' picking mode',
                                        zoomx=self.getXLims(),
                                        zoomy=self.getYLims(),
                                        noiselevel=noiselevels,
                                        scaleddata=True)

    def setPick(self, gui_event):

        # get axes limits
        self.updateCurrentLimits()

        # setting pick
        pick = gui_event.xdata  # get pick time relative to the traces timeaxis not to the global
        channel = self.getChannelID(round(gui_event.ydata))

        wfdata = self.getWFData().copy().select(channel=channel)
        stime = self.getStartTime()
        # get earliest and latest possible pick and symmetric pick error
        [epp, lpp, spe] = earllatepicker(wfdata, 1.5, (5., .5, 2.), pick)

        # get name of phase actually picked
        phase = self.selectPhase.currentText()

        # return absolute time values for phases
        epp = stime + epp
        mpp = stime + pick
        lpp = stime + lpp

        # save pick times for actual phase
        phasepicks = {'epp': epp, 'lpp': lpp, 'mpp': mpp, 'spe': spe}

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
                     latest possible pick: {olpp}\n""".format(phase=phase,
                                                              epp=epp,
                                                              mpp=pick,
                                                              lpp=lpp,
                                                              oepp=oepp,
                                                              ompp=ompp,
                                                              olpp=olpp)
        self.getPlotWidget().plotWFData(wfdata=self.getWFData(),
                                        title=self.getStation())
        self.drawPicks()
        self.disconnectPressEvent()
        self.zoomAction.setEnabled(True)
        self.selectPhase.setCurrentIndex(-1)
        self.setPlotLabels()

    def drawPicks(self, phase=None):
        # plotting picks
        ax = self.getPlotWidget().axes
        ylims = self.getGlobalLimits('y')
        phase_col = {'P': ('c', 'c--', 'b-'),
                     'S': ('m', 'm--', 'r-')}
        if self.getPicks():
            if phase is not None and type(self.getPicks()[phase]) is dict:
                picks = self.getPicks()[phase]
                colors = phase_col[phase[0].upper()]
            elif phase is None:
                for phase in self.getPicks():
                    self.drawPicks(phase)
                return
            else:
                return
        else:
            return

        mpp = picks['mpp'] - self.getStartTime()
        epp = picks['epp'] - self.getStartTime()
        lpp = picks['lpp'] - self.getStartTime()
        spe = picks['spe']

        ax.fill_between([epp, lpp], ylims[0], ylims[1],
                        alpha=.5, color=colors[0])
        ax.plot([mpp - spe, mpp - spe], ylims, colors[1],
                [mpp, mpp], ylims, colors[2],
                [mpp + spe, mpp + spe], ylims, colors[1])

    def panPress(self, gui_event):
        ax = self.getPlotWidget().axes
        if gui_event.inaxes != ax: return
        self.cur_xlim = ax.get_xlim()
        self.cur_ylim = ax.get_ylim()
        self.press = gui_event.xdata, gui_event.ydata
        self.xpress, self.ypress = self.press

    def panRelease(self, gui_event):
        ax = self.getPlotWidget().axes
        self.press = None
        ax.figure.canvas.draw()

    def panMotion(self, gui_event):
        if self.press is None: return
        ax = self.getPlotWidget().axes
        if gui_event.inaxes != ax: return
        dx = gui_event.xdata - self.xpress
        dy = gui_event.ydata - self.ypress
        self.cur_xlim -= dx
        self.cur_ylim -= dy
        ax.set_xlim(self.cur_xlim)
        ax.set_ylim(self.cur_ylim)

        ax.figure.canvas.draw()

    def filterWFData(self):
        self.updateCurrentLimits()
        data = self.getWFData().copy()
        old_title = self.getPlotWidget().getAxes().get_title()
        title = None
        phase = self.selectPhase.currentText()
        filtoptions = None
        if phase:
            filtoptions = self.getFilterOptions(phase).parseFilterOptions()
        if self.filterAction.isChecked():
            if not phase:
                filtoptions = FilterOptionsDialog.getFilterObject()
                filtoptions = filtoptions.parseFilterOptions()
        if filtoptions is not None:
            data.filter(**filtoptions)
            if old_title.endswith(')'):
                title = old_title[:-1] + ', filtered)'
            else:
                title = old_title + ' (filtered)'
        else:
            if old_title.endswith(' (filtered)'):
                title = old_title.replace(' (filtered)', '')
        if title is None:
            title = old_title
        self.getPlotWidget().plotWFData(wfdata=data, title=title,
                                        zoomx=self.getXLims(),
                                        zoomy=self.getYLims())
        self.setPlotLabels()
        self.drawPicks()
        self.draw()

    def resetPlot(self):
        self.updateCurrentLimits()
        data = self.getWFData().copy()
        title = self.getPlotWidget().getAxes().get_title()
        self.getPlotWidget().plotWFData(wfdata=data, title=title,
                                        zoomx=self.getXLims(),
                                        zoomy=self.getYLims())
        self.setPlotLabels()
        self.drawPicks()
        self.draw()


    def setPlotLabels(self):

        # get channel labels
        pos = self.getPlotWidget().getPlotDict().keys()
        labels = [self.getPlotWidget().getPlotDict()[key][1] for key in pos]

        # set channel labels
        self.getPlotWidget().setYTickLabels(pos, labels)
        self.getPlotWidget().setXLims(self.getXLims())
        self.getPlotWidget().setYLims(self.getYLims())

    def zoom(self):
        if self.zoomAction.isChecked():
            self.disconnectPressEvent()
            self.disconnectMotionEvent()
            self.disconnectReleaseEvent()
            self.disconnectScrollEvent()
            self.figToolBar.zoom()
        else:
            self.figToolBar.zoom()
            self.cidpress = self.connectPressEvent(self.panPress)
            self.cidmotion = self.connectMotionEvent(self.panMotion)
            self.cidrelease = self.connectReleaseEvent(self.panRelease)
            self.cidscroll = self.connectScrollEvent(self.scrollZoom)

    def scrollZoom(self, gui_event, factor=2.):

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
                   scale_factor * (gui_event.xdata - self.getXLims())
        new_ylim = gui_event.ydata - \
                   scale_factor * (gui_event.ydata - self.getYLims())

        new_xlim.sort()
        global_x = self.getGlobalLimits('x')
        global_y = self.getGlobalLimits('y')
        new_xlim[0] = max(new_xlim[0], global_x[0])
        new_xlim[1] = min(new_xlim[1], global_x[1])
        new_ylim.sort()
        new_ylim[0] = max(new_ylim[0], global_y[0])
        new_ylim[1] = min(new_ylim[1], global_y[1])

        self.getPlotWidget().setXLims(new_xlim)
        self.getPlotWidget().setYLims(new_ylim)
        self.draw()

    def resetZoom(self):
        self.getPlotWidget().setXLims(self.getGlobalLimits('x'))
        self.getPlotWidget().setYLims(self.getGlobalLimits('y'))
        self.draw()

    def draw(self):
        self.getPlotWidget().draw()

    def apply(self):
        picks = self.getPicks()
        for pick in picks:
            print(pick, picks[pick])

    def accept(self):
        self.apply()
        QDialog.accept(self)


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
        self.tabWidget.addTab(LocalisationTab(self), "Loc Tools")
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok |
                                          QDialogButtonBox.Apply |
                                          QDialogButtonBox.Close)

        layout = QVBoxLayout()
        layout.addWidget(self.tabWidget)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(
            self.apply)

    def accept(self, *args, **kwargs):
        self.apply()
        QDialog.accept(self)

    def apply(self):
        for widint in range(self.tabWidget.count()):
            curwid = self.tabWidget.widget(widint)
            values = curwid.getValues()
            if values is not None:
                self.setValues(values)

    @staticmethod
    def setValues(tabValues):
        settings = QSettings()
        for setting, value in tabValues.items():
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
        try:
            self.fullNameEdit.setText(fulluser)
        except TypeError as e:
            self.fullNameEdit.setText(fulluser[0])

        # information about data structure
        dataroot = settings.value("data/dataRoot")
        curstructure = settings.value("data/Structure")
        dataDirLabel = QLabel("data root directory: ")
        self.dataDirEdit = QLineEdit()
        self.dataDirEdit.setText(dataroot)
        self.dataDirEdit.selectAll()
        structureLabel = QLabel("data structure: ")
        self.structureSelect = QComboBox()

        from pylot.core.util.structure import DATASTRUCTURE

        self.structureSelect.addItems(DATASTRUCTURE.keys())

        dsind = findComboBoxIndex(self.structureSelect, curstructure)

        self.structureSelect.setCurrentIndex(dsind)

        layout = QGridLayout()
        layout.addWidget(dataDirLabel, 0, 0)
        layout.addWidget(self.dataDirEdit, 0, 1)
        layout.addWidget(fullNameLabel, 1, 0)
        layout.addWidget(self.fullNameEdit, 1, 1)
        layout.addWidget(structureLabel, 2, 0)
        layout.addWidget(self.structureSelect, 2, 1)

        self.setLayout(layout)

    def getValues(self):
        values = {"data/dataRoot": self.dataDirEdit.text(),
                  "user/FullName": self.fullNameEdit.text(),
                  "data/Structure": self.structureSelect.currentText()}
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

        ind = findComboBoxIndex(self.eventOutputComboBox, curval)

        self.eventOutputComboBox.setCurrentIndex(ind)
        layout = QGridLayout()
        layout.addWidget(eventOutputLabel, 0, 0)
        layout.addWidget(self.eventOutputComboBox, 0, 1)

        self.setLayout(layout)

    def getValues(self):
        values = {"output/Format": self.eventOutputComboBox.currentText()}
        return values


class PhasesTab(PropTab):
    def __init__(self, parent=None):
        super(PhasesTab, self).__init__(parent)

        pass


class GraphicsTab(PropTab):
    def __init__(self, parent=None):
        super(GraphicsTab, self).__init__(parent)

        pass


class LocalisationTab(PropTab):
    def __init__(self, parent=None):
        super(LocalisationTab, self).__init__(parent)

        settings = QSettings()
        curtool = settings.value("loc/tool", None)

        loctoollabel = QLabel("location tool")
        self.locToolComboBox = QComboBox()
        loctools = LOCTOOLS.keys()
        self.locToolComboBox.addItems(loctools)

        toolind = findComboBoxIndex(self.locToolComboBox, curtool)

        self.locToolComboBox.setCurrentIndex(toolind)

        curroot = settings.value("%s/rootPath".format(curtool), None)
        curbin = settings.value("%s/binPath".format(curtool), None)

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

        self.locToolComboBox.currentIndexChanged.connect(self.updateUi)

        self.updateUi()

        layout = QGridLayout()
        layout.addWidget(loctoollabel, 0, 0)
        layout.addWidget(self.locToolComboBox, 0, 1)
        layout.addWidget(self.rootlabel, 1, 0)
        layout.addWidget(self.rootedit, 1, 1)
        layout.addWidget(rootBrowse, 1, 2)
        layout.addWidget(self.binlabel, 2, 0)
        layout.addWidget(self.binedit, 2, 1)
        layout.addWidget(binBrowse, 2, 2)

        self.setLayout(layout)

    def updateUi(self):
        curtool = self.locToolComboBox.currentText()
        if curtool is not None:
            self.rootlabel.setText("{0} root directory".format(curtool))
            self.binlabel.setText("{0} bin directory".format(curtool))

    def selectDirectory(self, edit):
        selected_directory =  QFileDialog.getExistingDirectory()
        edit.setText(selected_directory)

    def getValues(self):
        loctool = self.locToolComboBox.currentText()
        values = {"%s/rootPath".format(loctool): self.rootedit.text(),
                  "%s/binPath".format(loctool): self.binedit.text(),
                  "loc/tool": loctool}
        return values



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
        super(FilterOptionsDialog, self).__init__()

        if parent is not None and parent.getFilterOptions():
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
            if self.getFilterOptions().getFilterType() in ['bandpass',
                                                           'bandstop']:
                self.freqmaxSpinBox.setValue(
                    self.getFilterOptions().getFreq()[1])
        else:
            try:
                self.freqmaxSpinBox.setValue(self.getFilterOptions().getFreq())
                self.freqminSpinBox.setValue(self.getFilterOptions().getFreq())
            except TypeError as e:
                print(e)
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
        self.selectTypeCombo.setCurrentIndex(typeOptions.index(self.getFilterOptions().getFilterType()))
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

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok |
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
        type = self.selectTypeCombo.currentText()
        _enable = type in ['bandpass', 'bandstop']
        freq = [self.freqminSpinBox.value(), self.freqmaxSpinBox.value()]
        self.freqmaxLabel.setEnabled(_enable)
        self.freqmaxSpinBox.setEnabled(_enable)

        if not _enable:
            self.freqminLabel.setText("cutoff:")
            self.freqmaxSpinBox.setValue(freq[0])
            freq.remove(freq[1])
        else:
            self.freqminLabel.setText("minimum:")
            if not isSorted(freq):
                QMessageBox.warning(self, "Value error",
                                    "Maximum frequency must be at least the "
                                    "same value as minimum frequency (notch)!")
                self.freqmaxSpinBox.setValue(freq[0])
                self.freqmaxSpinBox.selectAll()
                self.freqmaxSpinBox.setFocus()
                return

        self.getFilterOptions().setFilterType(type)
        self.getFilterOptions().setFreq(freq)
        self.getFilterOptions().setOrder(self.orderSpinBox.value())

    def getFilterOptions(self):
        return self.filterOptions

    @staticmethod
    def getFilterObject():
        dlg = FilterOptionsDialog()
        if dlg.exec_():
            return dlg.getFilterOptions()
        return None

    def accept(self):
        self.updateUi()
        QDialog.accept(self)


class LoadDataDlg(QDialog):
    def __init__(self, parent=None):
        super(LoadDataDlg, self).__init__(parent)

        pass


class HelpForm(QDialog):
    def __init__(self, page=QUrl('https://ariadne.geophysik.rub.de/trac/PyLoT'),
                 parent=None):
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()
