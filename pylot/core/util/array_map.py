#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import matplotlib

from PySide2 import QtCore, QtGui, QtWidgets
from PySide2.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.patheffects as PathEffects

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cf
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter

import traceback
import obspy
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import griddata
from pylot.core.util.widgets import PickDlg
from pylot.core.pick.utils import get_quality_class

matplotlib.use('Qt5Agg')


class MplCanvas(FigureCanvas):

    def __init__(self, parent=None, extern_axes=None, width=5, height=4, dpi=100):
        if extern_axes is None:
            self.fig = plt.figure(figsize=(width, height), dpi=dpi)
            self.axes = self.fig.add_subplot(111)
        else:
            self.fig = extern_axes.figure
            self.axes = extern_axes

        super(MplCanvas, self).__init__(self.fig)


class Array_map(QtWidgets.QWidget):
    def __init__(self, parent, metadata, parameter=None, axes=None, annotate=True, pointsize=25.,
                 linewidth=1.5, width=5e6, height=2e6):
        QtWidgets.QWidget.__init__(self)

        assert (parameter is not None or parent is not None), 'either parent or parameter has to be set'
        # set properties
        self._parent = parent
        self.metadata = metadata
        self.pointsize = pointsize
        self.linewidth = linewidth
        self.extern_plot_axes = axes
        self.width = width
        self.height = height
        self.annotate = annotate
        self.picks = None
        self.picks_dict = None
        self.uncertainties = None
        self.autopicks_dict = None
        self.hybrids_dict = None
        self.eventLoc = None
        self.parameter = parameter if parameter else parent._inputs

        self.picks_rel = {}
        self.marked_stations = []
        self.highlighted_stations = []

        # call functions to draw everything
        self.init_graphics()
        self.init_stations()
        self.init_crtpyMap()
        self.init_map()
        # set original map limits to fall back on when home button is pressed
        self.org_xlim = self.canvas.axes.get_xlim()
        self.org_ylim = self.canvas.axes.get_ylim()
        # initial map without event
        self.canvas.axes.set_xlim(self.org_xlim[0], self.org_xlim[1])
        self.canvas.axes.set_ylim(self.org_ylim[0], self.org_ylim[1])

        self._style = None if not hasattr(parent, '_style') else parent._style

        self.show()

    def init_map(self):
        self.init_colormap()
        self.connectSignals()
        self.draw_everything()

    def init_graphics(self):
        """
        Initializes all GUI components and figure elements to be populated by other functions
        """
        # initialize figure elements

        if self.extern_plot_axes is None:
            self.canvas = MplCanvas(self)
            self.plotWidget = FigureCanvas(self.canvas.fig)
        else:
            self.canvas = MplCanvas(self, extern_axes=self.extern_plot_axes)
            self.plotWidget = FigureCanvas(self.canvas.fig)

        # initialize GUI elements
        self.status_label = QtWidgets.QLabel()
        self.map_reset_button = QtWidgets.QPushButton('Reset Map View')
        self.save_map_button = QtWidgets.QPushButton('Save Map')
        self.go2eq_button = QtWidgets.QPushButton('Go to Event Location')

        self.main_box = QtWidgets.QVBoxLayout()
        self.setLayout(self.main_box)

        self.top_row = QtWidgets.QHBoxLayout()
        self.main_box.addLayout(self.top_row, 1)

        self.comboBox_phase = QtWidgets.QComboBox()
        self.comboBox_phase.insertItem(0, 'P')
        self.comboBox_phase.insertItem(1, 'S')

        self.comboBox_am = QtWidgets.QComboBox()
        self.comboBox_am.insertItem(0, 'hybrid (prefer manual)')
        self.comboBox_am.insertItem(1, 'manual')
        self.comboBox_am.insertItem(2, 'auto')

        self.annotations_box = QtWidgets.QCheckBox('Annotate')
        self.annotations_box.setChecked(True)
        self.auto_refresh_box = QtWidgets.QCheckBox('Automatic refresh')
        self.auto_refresh_box.setChecked(True)
        self.refresh_button = QtWidgets.QPushButton('Refresh')
        self.cmaps_box = QtWidgets.QComboBox()
        self.cmaps_box.setMaxVisibleItems(20)
        [self.cmaps_box.addItem(map_name) for map_name in sorted(plt.colormaps())]
        # try to set to hsv as default
        self.cmaps_box.setCurrentIndex(self.cmaps_box.findText('hsv'))

        self.top_row.addWidget(QtWidgets.QLabel('Select a phase: '))
        self.top_row.addWidget(self.comboBox_phase)
        self.top_row.setStretch(1, 1)  # set stretch of item 1 to 1
        self.top_row.addWidget(QtWidgets.QLabel('Pick type: '))
        self.top_row.addWidget(self.comboBox_am)
        self.top_row.setStretch(3, 1)  # set stretch of item 1 to 1
        self.top_row.addWidget(self.cmaps_box)
        self.top_row.addWidget(self.annotations_box)
        self.top_row.addWidget(self.auto_refresh_box)
        self.top_row.addWidget(self.refresh_button)

        self.main_box.addWidget(self.plotWidget, 1)

        self.bot_row = QtWidgets.QHBoxLayout()
        self.main_box.addLayout(self.bot_row, 0.3)
        self.bot_row.addWidget(self.map_reset_button)
        self.bot_row.addWidget(self.go2eq_button)
        self.bot_row.addWidget(self.save_map_button)
        self.bot_row.addWidget(self.status_label)

    def init_colormap(self):
        self.init_lat_lon_dimensions()
        self.init_lat_lon_grid()

    def init_crtpyMap(self):
        self.canvas.axes.cla()
        self.canvas.axes = plt.axes(projection=ccrs.PlateCarree())
        self.canvas.axes.add_feature(cf.LAND)
        self.canvas.axes.add_feature(cf.OCEAN)
        self.canvas.axes.add_feature(cf.COASTLINE, linewidth=1, edgecolor='gray')
        self.canvas.axes.add_feature(cf.BORDERS, alpha=0.7)
        self.canvas.axes.add_feature(cf.LAKES, alpha=0.7)
        self.canvas.axes.add_feature(cf.RIVERS, linewidth=1)

        # parallels and meridians
        self.add_merid_paral()

        self.canvas.fig.tight_layout()

    def add_merid_paral(self):
        self.gridlines = self.canvas.axes.gridlines(draw_labels=False, alpha=0.6, color='gray', linewidth=self.linewidth/2, zorder=7)
        # TODO: current cartopy version does not support label removal. Devs are working on it.
        #  Should be fixed in coming cartopy versions
        # self.gridlines.xformatter = LONGITUDE_FORMATTER
        # self.gridlines.yformatter = LATITUDE_FORMATTER

    def remove_merid_paral(self):
        if len(self.gridlines.xline_artists):
            self.gridlines.xline_artists[0].remove()
            self.gridlines.yline_artists[0].remove()

    def org_map_view(self):
        self.canvas.axes.set_xlim(self.org_xlim[0], self.org_xlim[1])
        self.canvas.axes.set_ylim(self.org_ylim[0], self.org_ylim[1])
        # parallels and meridians
        self.remove_merid_paral()
        self.add_merid_paral()

        self.canvas.axes.figure.canvas.draw_idle()

    def go2eq(self):
        if self.eventLoc:
            lats, lons = self.eventLoc
            self.canvas.axes.set_xlim(lons - 10, lons + 10)
            self.canvas.axes.set_ylim(lats - 5, lats + 5)
            # parallels and meridians
            self.remove_merid_paral()
            self.add_merid_paral()

            self.canvas.axes.figure.canvas.draw_idle()

        else:
            self.status_label.setText('No event information available')

    def connectSignals(self):
        self.comboBox_phase.currentIndexChanged.connect(self._refresh_drawings)
        self.comboBox_am.currentIndexChanged.connect(self._refresh_drawings)
        self.cmaps_box.currentIndexChanged.connect(self._refresh_drawings)
        self.annotations_box.stateChanged.connect(self.switch_annotations)
        self.refresh_button.clicked.connect(self._refresh_drawings)
        self.map_reset_button.clicked.connect(self.org_map_view)
        self.go2eq_button.clicked.connect(self.go2eq)

        self.plotWidget.mpl_connect('motion_notify_event', self.mouse_moved)
        self.plotWidget.mpl_connect('scroll_event', self.mouse_scroll)
        self.plotWidget.mpl_connect('button_press_event', self.mouseLeftPress)
        self.plotWidget.mpl_connect('button_release_event', self.mouseLeftRelease)

    # set mouse events -----------------------------------------------------
    def mouse_moved(self, event):
        if not event.inaxes == self.canvas.axes:
            return
        lat = event.ydata
        lon = event.xdata
        self.status_label.setText('Latitude: {}, Longitude: {}'.format(lat, lon))

    def mouse_scroll(self, event):
        if not event.inaxes == self.canvas.axes:
            return

        zoom = {'up': 1. / 2., 'down': 2.}

        if event.button in zoom:
            xlim = self.canvas.axes.get_xlim()
            ylim = self.canvas.axes.get_ylim()

            x, y = event.xdata, event.ydata

            factor = zoom[event.button]
            xdiff = (xlim[1] - xlim[0]) * factor
            xl = x - 0.5 * xdiff
            xr = x + 0.5 * xdiff
            ydiff = (ylim[1] - ylim[0]) * factor
            yb = y - 0.5 * ydiff
            yt = y + 0.5 * ydiff

            self.canvas.axes.set_xlim(xl, xr)
            self.canvas.axes.set_ylim(yb, yt)
            # parallels and meridians
            self.remove_merid_paral()
            self.add_merid_paral()

            self.canvas.axes.figure.canvas.draw_idle()

    def mouseLeftPress(self, event):
        if not event.inaxes == self.canvas.axes:
            return
        self.map_x = event.xdata
        self.map_y = event.ydata
        self.map_xlim = self.canvas.axes.get_xlim()
        self.map_ylim = self.canvas.axes.get_ylim()

    def mouseLeftRelease(self, event):
        if not event.inaxes == self.canvas.axes:
            return
        new_x = event.xdata
        new_y = event.ydata

        dx = new_x - self.map_x
        dy = new_y - self.map_y

        self.canvas.axes.set_xlim((self.map_xlim[0] - dx, self.map_xlim[1] - dx))
        self.canvas.axes.set_ylim(self.map_ylim[0] - dy, self.map_ylim[1] - dy)
        # parallels and meridians
        self.remove_merid_paral()
        self.add_merid_paral()

        self.canvas.axes.figure.canvas.draw_idle()

    def onpick(self, event):
        ind = event.ind
        button = event.mouseevent.button
        if ind == []:
            return
        if button == 1:
            self.openPickDlg(ind)
        elif button == 2:
            self.deletePick(ind)
        elif button == 3:
            self.pickInfo(ind)

    # data handling -----------------------------------------------------
    def update_hybrids_dict(self):
        self.hybrids_dict = self.picks_dict.copy()
        for station, pick in self.autopicks_dict.items():
            if not station in self.hybrids_dict.keys():
                self.hybrids_dict[station] = pick
        return self.hybrids_dict

    def deletePick(self, ind):
        self.update_hybrids_dict()
        for index in ind:
            network, station = self._station_onpick_ids[index].split('.')[:2]
            try:
                phase = self.comboBox_phase.currentText()
                picks = self.current_picks_dict()[station]
                pick = picks.get(phase)
                if pick:
                    picker = pick['picker']
                    message = 'Deleted {} pick for phase {}, station {}.{} at timestamp {}'
                    message = message.format(picker, phase, network, station,
                                             pick['mpp'])
                    if picker == 'auto':
                        del (self.autopicks_dict[station])
                    elif picker == 'manual':
                        del (self.picks_dict[station])
                    else:
                        raise TypeError('Unknown "picker" {}'.format(picker))
                    print(message)
                    pyl_mw = self._parent
                    pyl_mw.deletePicks(station, pick, type=picker)
                    pyl_mw.setDirty(True)
                    pyl_mw.update_status(message)
                    if self.auto_refresh_box.isChecked():
                        self._refresh_drawings()
                    else:
                        self.highlight_station(network, station, color='red')
                    pyl_mw.drawPicks(station)
                    pyl_mw.draw()
            except Exception as e:
                print('Could not delete pick for station {}.{}: {}'.format(network, station, e))

    def pickInfo(self, ind):
        self.update_hybrids_dict()
        for index in ind:
            network, station = self._station_onpick_ids[index].split('.')[:2]
            dic = self.current_picks_dict()[station]
            for phase, picks in dic.items():
                # because of wadati...
                if phase == 'SPt':
                    continue
                print('{} - Pick:'.format(phase))
                for key, info in picks.items():
                    print('{}: {}'.format(key, info))

    def _from_dict(self, function, key):
        return function(self.stations_dict.values(), key=lambda x: x[key])[key]

    def get_min_from_stations(self, key):
        return self._from_dict(min, key)

    def get_max_from_stations(self, key):
        return self._from_dict(max, key)

    def get_min_from_picks(self):
        return min(self.picks_rel.values())

    def get_max_from_picks(self):
        return max(self.picks_rel.values())

    def current_picks_dict(self):
        picktype = self.comboBox_am.currentText().split(' ')[0]
        auto_manu = {'auto': self.autopicks_dict,
                     'manual': self.picks_dict,
                     'hybrid': self.hybrids_dict}
        return auto_manu[picktype]

    def init_stations(self):
        self.stations_dict = self.metadata.get_all_coordinates()
        self.latmin = self.get_min_from_stations('latitude')
        self.lonmin = self.get_min_from_stations('longitude')
        self.latmax = self.get_max_from_stations('latitude')
        self.lonmax = self.get_max_from_stations('longitude')

    def init_picks(self):
        def get_picks(station_dict):
            self.update_hybrids_dict()
            picks = {}
            uncertainties = {}
            # selected phase
            phase = self.comboBox_phase.currentText()
            for st_id in station_dict.keys():
                try:
                    station_name = st_id.split('.')[-1]
                    # current_picks_dict: auto or manual
                    pick = self.current_picks_dict()[station_name][phase]
                    if pick['picker'] == 'auto':
                        if not pick['spe']:
                            continue
                    picks[st_id] = pick['mpp']
                    uncertainties[st_id] = pick['spe']
                except KeyError:
                    continue
                except Exception as e:
                    print('Cannot display pick for station {}. Reason: {}'.format(station_name, e))
            return picks, uncertainties

        def get_picks_rel(picks):
            picks_rel = {}
            picks_utc = []
            for pick in picks.values():
                if type(pick) is obspy.core.utcdatetime.UTCDateTime:
                    picks_utc.append(pick)
            if picks_utc:
                self._earliest_picktime = min(picks_utc)
                for st_id, pick in picks.items():
                    if type(pick) is obspy.core.utcdatetime.UTCDateTime:
                        pick -= self._earliest_picktime
                    picks_rel[st_id] = pick
            return picks_rel

        self.picks, self.uncertainties = get_picks(self.stations_dict)
        self.picks_rel = get_picks_rel(self.picks)

    def init_lat_lon_dimensions(self):
        # init minimum and maximum lon and lat dimensions
        self.londim = self.lonmax - self.lonmin
        self.latdim = self.latmax - self.latmin

    def init_lat_lon_grid(self, nstep=250):
        # create a regular grid to display colormap
        lataxis = np.linspace(self.latmin, self.latmax, nstep)
        lonaxis = np.linspace(self.lonmin, self.lonmax, nstep)
        self.longrid, self.latgrid = np.meshgrid(lonaxis, lataxis)

    def init_picksgrid(self):
        picks, uncertainties, lats, lons = self.get_picks_lat_lon()
        try:
            self.picksgrid_active = griddata((lats, lons), picks, (self.latgrid, self.longrid), method='linear')
        except Exception as e:
            self._warn('Could not init picksgrid: {}'.format(e))

    def get_st_lat_lon_for_plot(self):
        stations = []
        latitudes = []
        longitudes = []
        for st_id, coords in self.stations_dict.items():
            stations.append(st_id)
            latitudes.append(coords['latitude'])
            longitudes.append(coords['longitude'])
        return stations, latitudes, longitudes

    def get_picks_lat_lon(self):
        picks = []
        uncertainties = []
        latitudes = []
        longitudes = []
        for st_id, pick in self.picks_rel.items():
            picks.append(pick)
            uncertainties.append(self.uncertainties.get(st_id))
            latitudes.append(self.stations_dict[st_id]['latitude'])
            longitudes.append(self.stations_dict[st_id]['longitude'])
        return picks, uncertainties, latitudes, longitudes

    # plotting -----------------------------------------------------
    def highlight_station(self, network, station, color):
        stat_dict = self.stations_dict['{}.{}'.format(network, station)]
        lat = stat_dict['latitude']
        lon = stat_dict['longitude']
        self.highlighted_stations.append(self.canvas.axes.scatter(lon, lat, s=self.pointsize, edgecolors=color,
                                                                  facecolors='none', zorder=12,
                                                                  transform=ccrs.PlateCarree(), label='deleted'))

    def openPickDlg(self, ind):
        data = self._parent.get_data().getWFData()
        for index in ind:
            network, station = self._station_onpick_ids[index].split('.')[:2]
            pyl_mw = self._parent
            try:
                data = data.select(station=station)
                if not data:
                    self._warn('No data for station {}'.format(station))
                    return
                pickDlg = PickDlg(self._parent, parameter=self.parameter,
                                  data=data, network=network, station=station,
                                  picks=self._parent.get_current_event().getPick(station),
                                  autopicks=self._parent.get_current_event().getAutopick(station),
                                  filteroptions=self._parent.filteroptions, metadata=self.metadata,
                                  event=pyl_mw.get_current_event())
            except Exception as e:
                message = 'Could not generate Plot for station {st}.\n {er}'.format(st=station, er=e)
                self._warn(message)
                print(message, e)
                print(traceback.format_exc())
                return
            try:
                if pickDlg.exec_():
                    pyl_mw.setDirty(True)
                    pyl_mw.update_status('picks accepted ({0})'.format(station))
                    pyl_mw.addPicks(station, pickDlg.getPicks(picktype='manual'), type='manual')
                    pyl_mw.addPicks(station, pickDlg.getPicks(picktype='auto'), type='auto')
                    if self.auto_refresh_box.isChecked():
                        self._refresh_drawings()
                    else:
                        self.highlight_station(network, station, color='yellow')
                    pyl_mw.drawPicks(station)
                    pyl_mw.draw()
                else:
                    pyl_mw.update_status('picks discarded ({0})'.format(station))
            except Exception as e:
                message = 'Could not save picks for station {st}.\n{er}'.format(st=station, er=e)
                self._warn(message)
                print(message, e)
                print(traceback.format_exc())

    def draw_contour_filled(self, nlevel=50):
        levels = np.linspace(self.get_min_from_picks(), self.get_max_from_picks(), nlevel)

        self.contourf = self.canvas.axes.contourf(self.longrid, self.latgrid, self.picksgrid_active, levels,
                                                  linewidths=self.linewidth*5, transform=ccrs.PlateCarree(),
                                                  alpha=0.4, zorder=8, cmap=self.get_colormap())

    def get_colormap(self):
        return plt.get_cmap(self.cmaps_box.currentText())

    def scatter_all_stations(self):
        stations, lats, lons = self.get_st_lat_lon_for_plot()

        self.sc = self.canvas.axes.scatter(lons, lats, s=self.pointsize * 3, facecolor='none', marker='.',
                                           zorder=10, picker=True, edgecolor='0.5', label='Not Picked',
                                           transform=ccrs.PlateCarree())

        self.cid = self.plotWidget.mpl_connect('pick_event', self.onpick)
        self._station_onpick_ids = stations
        if self.eventLoc:
            lats, lons = self.eventLoc
            self.sc_event = self.canvas.axes.scatter(lons, lats, s=5 * self.pointsize, facecolor='red', zorder=11,
                                                     label='Event (might be outside map region)', marker='*',
                                                     edgecolors='black',
                                                     transform=ccrs.PlateCarree())

    def scatter_picked_stations(self):
        picks, uncertainties, lats, lons = self.get_picks_lat_lon()
        if len(lons) < 1 and len(lats) < 1:
            return

        phase = self.comboBox_phase.currentText()
        timeerrors = self.parameter['timeerrors{}'.format(phase)]
        sizes = np.array([self.pointsize * (5. - get_quality_class(uncertainty, timeerrors))
                          for uncertainty in uncertainties])

        cmap = self.get_colormap()
        self.sc_picked = self.canvas.axes.scatter(lons, lats, s=sizes, edgecolors='white', cmap=cmap,
                                                  c=picks, zorder=11, label='Picked', transform=ccrs.PlateCarree())

    def annotate_ax(self):
        self.annotations = []
        stations, ys, xs = self.get_st_lat_lon_for_plot()
        # MP MP testing station highlighting if they have high impact on mean gradient of color map
        # if self.picks_rel:
        #    self.test_gradient()
        color_marked = {True: 'red',
                        False: 'white'}
        for st, x, y in zip(stations, xs, ys):
            if st in self.picks_rel:
                color = 'white'
            else:
                color = 'lightgrey'
            if st in self.marked_stations:
                color = 'red'
            self.annotations.append(
                self.canvas.axes.annotate(' %s' % st, xy=(x + 0.003, y + 0.003), fontsize=self.pointsize / 4.,
                                          fontweight='semibold', color=color, alpha=0.8,
                                          transform=ccrs.PlateCarree(), zorder=14,
                                          path_effects=[PathEffects.withStroke(
                                              linewidth=self.pointsize / 15., foreground='k')]))

        self.legend = self.canvas.axes.legend(loc=1, framealpha=1)
        self.legend.set_zorder(100)
        self.legend.get_frame().set_facecolor((1, 1, 1, 0.95))

    def add_cbar(self, label):
        self.cbax_bg = inset_axes(self.canvas.axes, width="6%", height="75%", loc=5)
        cbax = inset_axes(self.canvas.axes, width='2%', height='70%', loc=5)
        cbar = self.canvas.axes.figure.colorbar(self.sc_picked, cax=cbax)
        cbar.set_label(label)
        cbax.yaxis.tick_left()
        cbax.yaxis.set_label_position('left')
        for spine in self.cbax_bg.spines.values():
            spine.set_visible(False)
        self.cbax_bg.yaxis.set_ticks([])
        self.cbax_bg.xaxis.set_ticks([])
        self.cbax_bg.patch.set_facecolor((1, 1, 1, 0.75))
        return cbar

    # handle drawings -----------------------------------------------------
    def refresh_drawings(self, picks=None, autopicks=None):
        self.picks_dict = picks
        self.autopicks_dict = autopicks
        self._refresh_drawings()

    def _refresh_drawings(self):
        self.remove_drawings()
        self.init_stations()
        self.init_colormap()
        self.draw_everything()

    def switch_annotations(self):
        if self.annotations_box.isChecked():
            self.annotate = True
        else:
            self.annotate = False
        self._refresh_drawings()

    def draw_everything(self):
        picktype = self.comboBox_am.currentText()
        picks_available = (self.picks_dict and picktype == 'manual') \
                          or (self.autopicks_dict and picktype == 'auto') \
                          or ((self.autopicks_dict or self.picks_dict) and picktype.startswith('hybrid'))

        if picks_available:
            self.init_picks()
            if len(self.picks) >= 3:
                self.init_picksgrid()
                self.draw_contour_filled()
        self.scatter_all_stations()
        if picks_available:
            self.scatter_picked_stations()
            if hasattr(self, 'sc_picked'):
                self.cbar = self.add_cbar(label='Time relative to first onset ({}) [s]'.format(self._earliest_picktime))
            self.comboBox_phase.setEnabled(True)
        else:
            self.comboBox_phase.setEnabled(False)
        if self.annotate:
            self.annotate_ax()
        self.plotWidget.draw_idle()

    def remove_drawings(self):
        self.remove_annotations()
        for item in reversed(self.highlighted_stations):
            item.remove()
            self.highlighted_stations.remove(item)
        if hasattr(self, 'cbar'):
            try:
                self.cbar.remove()
                self.cbax_bg.remove()
            except Exception as e:
                print('Warning: could not remove color bar or color bar bg.\nReason: {}'.format(e))
            del (self.cbar, self.cbax_bg)
        if hasattr(self, 'sc_picked'):
            self.sc_picked.remove()
            del self.sc_picked
        if hasattr(self, 'sc_event'):
            self.sc_event.remove()
            del self.sc_event
        if hasattr(self, 'contourf'):
            self.remove_contourf()
            del self.contourf
        if hasattr(self, 'cid'):
            self.plotWidget.mpl_disconnect(self.cid)
            del self.cid
        try:
            self.sc.remove()
        except Exception as e:
            print('Warning: could not remove station scatter plot.\nReason: {}'.format(e))
        try:
            self.legend.remove()
        except Exception as e:
            print('Warning: could not remove legend. Reason: {}'.format(e))
        self.plotWidget.draw_idle()

    def remove_contourf(self):
        for item in self.contourf.collections:
            item.remove()

    def remove_annotations(self):
        for annotation in self.annotations:
            annotation.remove()
        self.annotations = []

    def _warn(self, message):
        self.qmb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Warning, 'Warning', message)
        self.qmb.show()
