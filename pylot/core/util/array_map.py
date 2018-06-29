#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import obspy
from PySide import QtGui
from mpl_toolkits.basemap import Basemap
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pylot.core.util.widgets import PickDlg, PylotCanvas
from scipy.interpolate import griddata

plt.interactive(False)


class Array_map(QtGui.QWidget):
    def __init__(self, parent, figure=None):
        '''
        Create a map of the array.
        :param parent: PyLoT Mainwindow class
        :param figure:
        '''
        QtGui.QWidget.__init__(self)
        self._parent = parent
        self.metadata_type = parent.metadata[0]
        self.metadata = parent.metadata[1]
        self.picks = None
        self.picks_dict = None
        self.autopicks_dict = None
        self.eventLoc = None
        self.figure = figure
        self.init_graphics()
        self.init_stations()
        self.init_basemap(resolution='l')
        self.init_map()
        self._style = parent._style
        # self.show()

    def init_map(self):
        self.init_lat_lon_dimensions()
        self.init_lat_lon_grid()
        self.init_x_y_dimensions()
        self.connectSignals()
        self.draw_everything()
        self.canvas.setZoomBorders2content()

    def onpick(self, event):
        ind = event.ind
        button = event.mouseevent.button
        if ind == [] or not button == 1:
            return
        data = self._parent.get_data().getWFData()
        for index in ind:
            station = str(self._station_onpick_ids[index].split('.')[-1])
            try:
                data = data.select(station=station)
                if not data:
                    self._warn('No data for station {}'.format(station))
                    return
                pickDlg = PickDlg(self._parent, parameter=self._parent._inputs,
                                  data=data,
                                  station=station,
                                  picks=self._parent.get_current_event().getPick(station),
                                  autopicks=self._parent.get_current_event().getAutopick(station),
                                  filteroptions=self._parent.filteroptions)
            except Exception as e:
                message = 'Could not generate Plot for station {st}.\n {er}'.format(st=station, er=e)
                self._warn(message)
                print(message, e)
                return
            pyl_mw = self._parent
            try:
                if pickDlg.exec_():
                    pyl_mw.setDirty(True)
                    pyl_mw.update_status('picks accepted ({0})'.format(station))
                    replot = pyl_mw.get_current_event().setPick(station, pickDlg.getPicks())
                    self._refresh_drawings()
                    if replot:
                        pyl_mw.plotWaveformData()
                        pyl_mw.drawPicks()
                        pyl_mw.draw()
                    else:
                        pyl_mw.drawPicks(station)
                        pyl_mw.draw()
                else:
                    pyl_mw.update_status('picks discarded ({0})'.format(station))
            except Exception as e:
                message = 'Could not save picks for station {st}.\n{er}'.format(st=station, er=e)
                self._warn(message)
                print(message, e)

    def connectSignals(self):
        self.comboBox_phase.currentIndexChanged.connect(self._refresh_drawings)
        self.comboBox_am.currentIndexChanged.connect(self._refresh_drawings)
        self.canvas.mpl_connect('motion_notify_event', self.mouse_moved)
        #self.zoom_id = self.basemap.ax.figure.canvas.mpl_connect('scroll_event', self.zoom)

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

    def mouse_moved(self, event):
        if not event.inaxes == self.main_ax:
            return
        x = event.xdata
        y = event.ydata
        lat, lon = self.basemap(x, y, inverse=True)
        self.status_label.setText('Latitude: {}, Longitude: {}'.format(lat, lon))

    def current_picks_dict(self):
        picktype = self.comboBox_am.currentText()
        auto_manu = {'auto': self.autopicks_dict,
                     'manual': self.picks_dict}
        return auto_manu[picktype]

    def init_graphics(self):
        if not self.figure:
            self.figure = Figure()

        self.status_label = QtGui.QLabel()

        self.main_ax = self.figure.add_subplot(111)
        self.canvas = PylotCanvas(self.figure, parent=self._parent, multicursor=True,
                                  panZoomX=False, panZoomY=False)

        self.main_box = QtGui.QVBoxLayout()
        self.setLayout(self.main_box)

        self.top_row = QtGui.QHBoxLayout()
        self.main_box.addLayout(self.top_row, 1)

        self.comboBox_phase = QtGui.QComboBox()
        self.comboBox_phase.insertItem(0, 'P')
        self.comboBox_phase.insertItem(1, 'S')

        self.comboBox_am = QtGui.QComboBox()
        self.comboBox_am.insertItem(0, 'auto')
        self.comboBox_am.insertItem(1, 'manual')

        self.top_row.addWidget(QtGui.QLabel('Select a phase: '))
        self.top_row.addWidget(self.comboBox_phase)
        self.top_row.setStretch(1, 1)  # set stretch of item 1 to 1
        self.top_row.addWidget(QtGui.QLabel('Pick type: '))
        self.top_row.addWidget(self.comboBox_am)
        self.top_row.setStretch(3, 1)  # set stretch of item 1 to 1

        self.main_box.addWidget(self.canvas, 1)
        self.main_box.addWidget(self.status_label, 0)


    def init_stations(self):
        def stat_info_from_parser(parser):
            stations_dict = {}
            for station in parser.stations:
                station_name = station[0].station_call_letters
                network_name = station[0].network_code
                if not station_name in stations_dict.keys():
                    st_id = network_name + '.' + station_name
                    stations_dict[st_id] = {'latitude': station[0].latitude,
                                            'longitude': station[0].longitude}
            return stations_dict

        def stat_info_from_inventory(inventory):
            stations_dict = {}
            for network in inventory.networks:
                for station in network.stations:
                    station_name = station.code
                    network_name = network_name.code
                    if not station_name in stations_dict.keys():
                        st_id = network_name + '.' + station_name
                        stations_dict[st_id] = {'latitude': station[0].latitude,
                                                'longitude': station[0].longitude}
            return stations_dict

        read_stat = {'xml': stat_info_from_inventory,
                     'dless': stat_info_from_parser}

        self.stations_dict = read_stat[self.metadata_type](self.metadata)
        self.latmin = self.get_min_from_stations('latitude')
        self.lonmin = self.get_min_from_stations('longitude')
        self.latmax = self.get_max_from_stations('latitude')
        self.lonmax = self.get_max_from_stations('longitude')

    def init_picks(self):
        def get_picks(station_dict):
            picks = {}
            phase = self.comboBox_phase.currentText()
            for st_id in station_dict.keys():
                try:
                    station_name = st_id.split('.')[-1]
                    pick = self.current_picks_dict()[station_name][phase]
                    if pick['picker'] == 'auto':
                        if pick['weight'] > 3:
                            continue
                    picks[st_id] = pick['mpp']
                except KeyError:
                    continue
                except Exception as e:
                    print('Cannot display pick for station {}. Reason: {}'.format(station_name, e))
            return picks

        def get_picks_rel(picks):
            picks_rel = {}
            picks_utc = []
            for pick in picks.values():
                if type(pick) is obspy.core.utcdatetime.UTCDateTime:
                    picks_utc.append(pick)
            self._earliest_picktime = min(picks_utc)
            for st_id, pick in picks.items():
                if type(pick) is obspy.core.utcdatetime.UTCDateTime:
                    pick -= self._earliest_picktime
                picks_rel[st_id] = pick
            return picks_rel

        self.picks = get_picks(self.stations_dict)
        self.picks_rel = get_picks_rel(self.picks)

    def init_lat_lon_dimensions(self):
        # init minimum and maximum lon and lat dimensions
        self.londim = self.lonmax - self.lonmin
        self.latdim = self.latmax - self.latmin

    def init_x_y_dimensions(self):
        # transformation of lat/lon to ax coordinate system
        for st_id, coords in self.stations_dict.items():
            lat, lon = coords['latitude'], coords['longitude']
            coords['x'], coords['y'] = self.basemap(lon, lat)

        self.xdim = self.get_max_from_stations('x') - self.get_min_from_stations('x')
        self.ydim = self.get_max_from_stations('y') - self.get_min_from_stations('y')

    def init_basemap(self, resolution='l'):
        # basemap = Basemap(projection=projection, resolution = resolution, ax=self.main_ax)
        width = 5e6
        height = 2e6
        basemap = Basemap(projection='lcc', resolution=resolution, ax=self.main_ax,
                          width=width, height=height,
                          lat_0=(self.latmin + self.latmax) / 2.,
                          lon_0=(self.lonmin + self.lonmax) / 2.)

        # basemap.fillcontinents(color=None, lake_color='aqua',zorder=1)
        basemap.drawmapboundary(zorder=2)  # fill_color='darkblue')
        basemap.shadedrelief(zorder=3)
        basemap.drawcountries(zorder=4)
        basemap.drawstates(zorder=5)
        basemap.drawcoastlines(zorder=6)
        # labels = [left,right,top,bottom]
        parallels = np.arange(-90, 90, 5.)
        parallels_small = np.arange(-90, 90, 2.5)
        basemap.drawparallels(parallels_small, linewidth=0.5, zorder=7)
        basemap.drawparallels(parallels, zorder=7)
        meridians = np.arange(-180, 180, 5.)
        meridians_small = np.arange(-180, 180, 2.5)
        basemap.drawmeridians(meridians_small, linewidth=0.5, zorder=7)
        basemap.drawmeridians(meridians, zorder=7)
        self.basemap = basemap
        self.figure._tight = True
        self.figure.tight_layout()

    def init_lat_lon_grid(self, nstep=250):
        # create a regular grid to display colormap
        lataxis = np.linspace(self.latmin, self.latmax, nstep)
        lonaxis = np.linspace(self.lonmin, self.lonmax, nstep)
        self.longrid, self.latgrid = np.meshgrid(lonaxis, lataxis)

    def init_picksgrid(self):
        picks, lats, lons = self.get_picks_lat_lon()
        try:
            self.picksgrid_active = griddata((lats, lons), picks, (self.latgrid, self.longrid),
                                             method='linear')
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

    def get_st_x_y_for_plot(self):
        stations = []
        xs = []
        ys = []
        for st_id, coords in self.stations_dict.items():
            stations.append(st_id)
            xs.append(coords['x'])
            ys.append(coords['y'])
        return stations, xs, ys

    def get_picks_lat_lon(self):
        picks = []
        latitudes = []
        longitudes = []
        for st_id, pick in self.picks_rel.items():
            picks.append(pick)
            latitudes.append(self.stations_dict[st_id]['latitude'])
            longitudes.append(self.stations_dict[st_id]['longitude'])
        return picks, latitudes, longitudes

    def draw_contour_filled(self, nlevel='50'):
        levels = np.linspace(self.get_min_from_picks(), self.get_max_from_picks(), nlevel)
        self.contourf = self.basemap.contourf(self.longrid, self.latgrid, self.picksgrid_active,
                                              levels, latlon=True, zorder=9, alpha=0.5)

    def scatter_all_stations(self):
        stations, lats, lons = self.get_st_lat_lon_for_plot()
        self.sc = self.basemap.scatter(lons, lats, s=50, facecolor='none', latlon=True,
                                       zorder=10, picker=True, edgecolor='m', label='Not Picked')
        self.cid = self.canvas.mpl_connect('pick_event', self.onpick)
        self._station_onpick_ids = stations
        if self.eventLoc:
            lats, lons = self.eventLoc
            self.sc_event = self.basemap.scatter(lons, lats, s=100, facecolor='red',
                                                 latlon=True, zorder=11, label='Event (might be outside map region)')

    def scatter_picked_stations(self):
        picks, lats, lons = self.get_picks_lat_lon()
        # workaround because of an issue with latlon transformation of arrays with len <3
        if len(lons) <= 2 and len(lats) <= 2:
            self.sc_picked = self.basemap.scatter(lons[0], lats[0], s=50, facecolor='white',
                                                  c=picks[0], latlon=True, zorder=11, label='Picked')
        if len(lons) == 2 and len(lats) == 2:
            self.sc_picked = self.basemap.scatter(lons[1], lats[1], s=50, facecolor='white',
                                                  c=picks[1], latlon=True, zorder=11)
        else:
            self.sc_picked = self.basemap.scatter(lons, lats, s=50, facecolor='white',
                                                  c=picks, latlon=True, zorder=11, label='Picked')

    def annotate_ax(self):
        self.annotations = []
        stations, xs, ys = self.get_st_x_y_for_plot()
        for st, x, y in zip(stations, xs, ys):
            self.annotations.append(self.main_ax.annotate(' %s' % st, xy=(x, y),
                                                          fontsize='x-small', color='white', zorder=12))
        self.legend = self.main_ax.legend(loc=1)
        self.legend.get_frame().set_facecolor((1, 1, 1, 0.75))

    def add_cbar(self, label):
        self.cbax_bg = inset_axes(self.main_ax, width="6%", height="75%", loc=5)
        cbax = inset_axes(self.main_ax, width='2%', height='70%', loc=5)
        cbar = self.main_ax.figure.colorbar(self.sc_picked, cax = cbax)
        cbar.set_label(label)
        cbax.yaxis.tick_left()
        cbax.yaxis.set_label_position('left')
        for spine in self.cbax_bg.spines.values():
            spine.set_visible(False)
        self.cbax_bg.yaxis.set_ticks([])
        self.cbax_bg.xaxis.set_ticks([])
        self.cbax_bg.patch.set_facecolor((1, 1, 1, 0.75))
        return cbar

    def refresh_drawings(self, picks=None, autopicks=None):
        self.picks_dict = picks
        self.autopicks_dict = autopicks
        self._refresh_drawings()

    def _refresh_drawings(self):
        self.remove_drawings()
        self.draw_everything()

    def draw_everything(self):
        if self.picks_dict or self.autopicks_dict:
            self.init_picks()
            if len(self.picks) >= 3:
                self.init_picksgrid()
                self.draw_contour_filled()
        self.scatter_all_stations()
        if self.picks_dict or self.autopicks_dict:
            self.scatter_picked_stations()
            self.cbar = self.add_cbar(label='Time relative to first onset ({}) [s]'.format(self._earliest_picktime))
            self.comboBox_phase.setEnabled(True)
        else:
            self.comboBox_phase.setEnabled(False)
        self.annotate_ax()
        self.canvas.draw()

    def remove_drawings(self):
        if hasattr(self, 'cbar'):
            self.cbar.remove()
            self.cbax_bg.remove()
            del (self.cbar, self.cbax_bg)
        if hasattr(self, 'sc_picked'):
            self.sc_picked.remove()
            del (self.sc_picked)
        if hasattr(self, 'sc_event'):
            self.sc_event.remove()
            del (self.sc_event)
        if hasattr(self, 'contourf'):
            self.remove_contourf()
            del (self.contourf)
        if hasattr(self, 'cid'):
            self.canvas.mpl_disconnect(self.cid)
            del (self.cid)
        try:
            self.sc.remove()
        except Exception as e:
            print('Warning: could not remove station scatter plot.\nReason: {}'.format(e))
        try:
            self.legend.remove()
        except Exception as e:
            print('Warning: could not remove legend. Reason: {}'.format(e))
        self.canvas.draw()

    def remove_contourf(self):
        for item in self.contourf.collections:
            item.remove()

    def remove_annotations(self):
        for annotation in self.annotations:
            annotation.remove()

    def zoom(self, event):
        map = self.basemap
        xlim = map.ax.get_xlim()
        ylim = map.ax.get_ylim()
        x, y = event.xdata, event.ydata
        zoom = {'up': 1. / 2.,
                'down': 2.}

        if not event.xdata or not event.ydata:
            return

        if event.button in zoom:
            factor = zoom[event.button]
            xdiff = (xlim[1] - xlim[0]) * factor
            xl = x - 0.5 * xdiff
            xr = x + 0.5 * xdiff
            ydiff = (ylim[1] - ylim[0]) * factor
            yb = y - 0.5 * ydiff
            yt = y + 0.5 * ydiff

            if xl < map.xmin or yb < map.ymin or xr > map.xmax or yt > map.ymax:
                xl, xr = map.xmin, map.xmax
                yb, yt = map.ymin, map.ymax
            map.ax.set_xlim(xl, xr)
            map.ax.set_ylim(yb, yt)
            map.ax.figure.canvas.draw()

    def _warn(self, message):
        self.qmb = QtGui.QMessageBox(QtGui.QMessageBox.Icon.Warning,
                                     'Warning', message)
        self.qmb.show()
