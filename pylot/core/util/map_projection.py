from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import obspy
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PySide import QtCore, QtGui

from pylot.core.util.widgets import PickDlg

plt.interactive(False)

class map_projection(QtGui.QWidget):
    def __init__(self, parent, figure=None):
        '''
        :param: picked, can be False, auto, manual
        :value: str
        '''
        QtGui.QWidget.__init__(self)
        self._parent = parent
        self.parser = parent.metadata[1]
        self.picks = None
        self.picks_dict = None
        self.figure = figure
        self.init_graphics()
        self.init_basemap(projection='mill', resolution='l')
        self.init_map()
        #self.show()

    def init_map(self):
        self.init_stations()
        self.init_lat_lon_dimensions()
        self.init_lat_lon_grid()
        self.init_x_y_dimensions()
        self.connectSignals()
        self.draw_everything()
        
    def onpick(self, event):
        ind = event.ind
        button = event.mouseevent.button
        if ind == [] or not button == 1:
            return
        data = self._parent.get_data().getWFData()
        for index in ind:
            station=str(self.station_names[index])
            try:
                pickDlg = PickDlg(self, parameter=self._parent._inputs,
                                  data=data.select(station=station),
                                  station=station,
                                  picks=self._parent.get_current_event().getPick(station),
                                  autopicks=self._parent.get_current_event().getAutopick(station))
            except Exception as e:
                message = 'Could not generate Plot for station {st}.\n{er}'.format(st=station, er=e)
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
        
    def init_graphics(self):
        if not self.figure:
            if not hasattr(self._parent, 'am_figure'):
                self.figure = plt.figure()
                self.toolbar = NavigationToolbar(self.figure.canvas, self)
            else:
                self.figure = self._parent.am_figure
                self.toolbar = self._parent.am_toolbar
                
        self.main_ax = self.figure.add_subplot(111)
        self.canvas = self.figure.canvas

        self.main_box = QtGui.QVBoxLayout()
        self.setLayout(self.main_box)

        self.top_row = QtGui.QHBoxLayout()
        self.main_box.addLayout(self.top_row)

        self.comboBox_phase = QtGui.QComboBox()
        self.comboBox_phase.insertItem(0, 'P')
        self.comboBox_phase.insertItem(1, 'S')

        self.comboBox_am = QtGui.QComboBox()
        self.comboBox_am.insertItem(0, 'auto')
        self.comboBox_am.insertItem(1, 'manual')        
        
        self.top_row.addWidget(QtGui.QLabel('Select a phase: '))
        self.top_row.addWidget(self.comboBox_phase)
        self.top_row.setStretch(1,1) #set stretch of item 1 to 1        

        self.main_box.addWidget(self.canvas)
        self.main_box.addWidget(self.toolbar)
        
    def init_stations(self):
        def get_station_names_lat_lon(parser):
            station_names=[]
            lat=[]
            lon=[]
            for station in parser.stations:
                station_name=station[0].station_call_letters
                if not station_name in station_names:
                    station_names.append(station_name)
                    lat.append(station[0].latitude)
                    lon.append(station[0].longitude)
            return station_names, lat, lon
        
        station_names, lat, lon = get_station_names_lat_lon(self.parser)
        self.station_names = station_names
        self.lat = lat
        self.lon = lon

    def init_picks(self):
        phase = self.comboBox_phase.currentText()
        def get_picks(station_names):
            picks=[]
            for station in station_names:
                try:
                    picks.append(self.picks_dict[station][phase]['mpp'])
                except:
                    picks.append(np.nan)
            return picks

        def get_picks_rel(picks):
            picks_rel=[]
            picks_utc = []
            for pick in picks:
                if type(pick) is obspy.core.utcdatetime.UTCDateTime:
                    picks_utc.append(pick)
            minp = min(picks_utc)
            for pick in picks:
                if type(pick) is obspy.core.utcdatetime.UTCDateTime:                
                    pick -= minp
                picks_rel.append(pick)
            return picks_rel
        
        self.picks = get_picks(self.station_names)
        self.picks_rel = get_picks_rel(self.picks)

    def init_picks_active(self):
        def remove_nan_picks(picks):
            picks_no_nan=[]
            for pick in picks:
                if not np.isnan(pick):
                    picks_no_nan.append(pick)
            return picks_no_nan
        
        self.picks_no_nan = remove_nan_picks(self.picks_rel)

    def init_stations_active(self):
        def remove_nan_lat_lon(picks, lat, lon):
            lat_no_nan=[]
            lon_no_nan=[]
            for index, pick in enumerate(picks):
                if not np.isnan(pick):
                    lat_no_nan.append(lat[index])
                    lon_no_nan.append(lon[index])
            return lat_no_nan, lon_no_nan
        
        self.lat_no_nan, self.lon_no_nan = remove_nan_lat_lon(self.picks_rel, self.lat, self.lon)

    def init_lat_lon_dimensions(self):
        def get_lon_lat_dim(lon, lat):
            londim = max(lon) - min(lon)
            latdim = max(lat) - min(lat)
            return londim, latdim
        
        self.londim, self.latdim = get_lon_lat_dim(self.lon, self.lat)

    def init_x_y_dimensions(self):
        def get_x_y_dim(x, y):
            xdim = max(x) - min(x)
            ydim = max(y) - min(y)
            return xdim, ydim
        
        self.x, self.y = self.basemap(self.lon, self.lat)
        self.xdim, self.ydim = get_x_y_dim(self.x, self.y)

    def init_basemap(self, projection, resolution='l'):
        basemap = Basemap(projection=projection, resolution = resolution, ax=self.main_ax)
        basemap.drawmapboundary(fill_color='darkblue')
        basemap.drawcountries()
        basemap.drawstates()
        basemap.fillcontinents(color='grey', lake_color='aqua')
        basemap.drawcoastlines()
        self.basemap = basemap
        self.figure.tight_layout()
        
    def init_lat_lon_grid(self):
        def get_lat_lon_axis(lat, lon):
            steplat = (max(lat)-min(lat))/250
            steplon = (max(lon)-min(lon))/250

            lataxis = np.arange(min(lat), max(lat), steplat)
            lonaxis = np.arange(min(lon), max(lon), steplon)
            return lataxis, lonaxis

        def get_lat_lon_grid(lataxis, lonaxis):
            longrid, latgrid = np.meshgrid(lonaxis, lataxis)
            return latgrid, longrid

        self.lataxis, self.lonaxis = get_lat_lon_axis(self.lat, self.lon)
        self.latgrid, self.longrid = get_lat_lon_grid(self.lataxis, self.lonaxis)

    def init_picksgrid(self):
        self.picksgrid_no_nan = griddata((self.lat_no_nan, self.lon_no_nan),
                                         self.picks_no_nan, (self.latgrid, self.longrid), method='linear') ##################

    def draw_contour_filled(self, nlevel='50'):
        levels = np.linspace(min(self.picks_rel), max(self.picks_rel), nlevel)
        self.contourf = self.basemap.contourf(self.longrid, self.latgrid, self.picksgrid_no_nan,
                                              levels, latlon=True, zorder=9)

    def scatter_all_stations(self):
        self.sc = self.basemap.scatter(self.lon, self.lat, s=50, facecolor='none', latlon=True,
                                  zorder=10, picker=True, edgecolor='m', label='Not Picked')
        self.cid = self.canvas.mpl_connect('pick_event', self.onpick)

    def scatter_picked_stations(self):
        lon = self.lon_no_nan
        lat = self.lat_no_nan

        #workaround because of an issue with latlon transformation of arrays with len <3
        if len(lon) <= 2 and len(lat) <= 2:
            self.sc_picked = self.basemap.scatter(lon[0], lat[0], s=50, facecolor='white',
                                                  c=self.picks_no_nan[0], latlon=True, zorder=11, label='Picked')
        if len(lon) == 2 and len(lat) == 2:            
            self.sc_picked = self.basemap.scatter(lon[1], lat[1], s=50, facecolor='white',
                                                  c=self.picks_no_nan[1], latlon=True, zorder=11)
        else:
            self.sc_picked = self.basemap.scatter(lon, lat, s=50, facecolor='white',
                                                  c=self.picks_no_nan, latlon=True, zorder=11, label='Picked')

    def annotate_ax(self):
        self.annotations=[]
        for index, name in enumerate(self.station_names):
            self.annotations.append(self.main_ax.annotate(' %s' % name, xy=(self.x[index], self.y[index]),
                                                          fontsize='x-small', color='white', zorder=12))
        self.legend=self.main_ax.legend()

    def add_cbar(self, label):
        cbar = self.main_ax.figure.colorbar(self.sc_picked)
        cbar.set_label(label)
        return cbar

    def refresh_drawings(self, picks=None):
        self.picks_dict = picks
        self.remove_drawings()
        self.draw_everything()

    def _refresh_drawings(self):
        self.remove_drawings()
        self.draw_everything()

    def draw_everything(self):
        if self.picks_dict:
            self.init_picks()
            self.init_picks_active()
            self.init_stations_active()
            if len(self.picks_no_nan) >= 3:
                self.init_picksgrid()
                self.draw_contour_filled()
        self.scatter_all_stations()
        if self.picks_dict:
            self.scatter_picked_stations()
            self.cbar = self.add_cbar(label='Time relative to first onset [s]')
            self.comboBox_phase.setEnabled(True)
        else:
            self.comboBox_phase.setEnabled(False)
        self.annotate_ax()
        self.canvas.draw()

    def remove_drawings(self):
        if hasattr(self, 'sc_picked'):
            self.sc_picked.remove()
            del(self.sc_picked)
        if hasattr(self, 'cbar'):
            self.cbar.remove()
            del(self.cbar)
        if hasattr(self, 'contourf'):
            self.remove_contourf()
            del(self.contourf)
        if hasattr(self, 'cid'):
            self.canvas.mpl_disconnect(self.cid)
            del(self.cid)
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
    
    def _warn(self, message):
        self.qmb = QtGui.QMessageBox(QtGui.QMessageBox.Icon.Warning,
                                     'Warning', message)
        self.qmb.show()        
            
    
