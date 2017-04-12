from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import obspy
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PySide import QtCore, QtGui

from pylot.core.util.dataprocessing import read_metadata
from pylot.core.util.widgets import PickDlg

plt.interactive(False)

class map_projection(QtGui.QWidget):
    def __init__(self, mainwindow):
        QtGui.QWidget.__init__(self)
        self.pyl_mainwindow = mainwindow
        self.parser = self.get_metadata('/data/Geothermie/Insheim/STAT_INFO/MAGS2_net.dless')
        self.init_graphics()
        self.init_stations()
        self.init_lat_lon_dimensions()
        self.init_lat_lon_grid()
        self.init_basemap(projection='mill', resolution='l')
        self.init_x_y_dimensions()
        self.connectSignals()
        self.draw_everything()
        self.show()
        
    def onpick(self, event):
        ind = event.ind
        if ind == []:
            return
        data = self.pyl_mainwindow.get_data().getWFData()
        for index in ind:
            station=str(self.station_names[index])
            try:
                pickDlg = PickDlg(self, infile=self.pyl_mainwindow.getinfile(), 
                                  data=data.select(station=station),
                                  station=station,
                                  picks=self.pyl_mainwindow.getPicksOnStation(station, 'manual'),
                                  autopicks=self.pyl_mainwindow.getPicksOnStation(station, 'auto'))
                pyl_mw = self.pyl_mainwindow
                if pickDlg.exec_():
                    pyl_mw.setDirty(True)
                    pyl_mw.update_status('picks accepted ({0})'.format(station))
                    replot = pyl_mw.addPicks(station, pickDlg.getPicks())
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
                print('Could not generate Plot for station {st}.\n{er}'.format(st=station, er=e))

    def get_metadata(self, path):
        metadata=read_metadata(path)
        parser=metadata[1]
        return parser

    def connectSignals(self):
        self.combobox.currentIndexChanged.connect(self.refresh_drawings)
        
    def init_graphics(self):
        self.main_box = QtGui.QVBoxLayout()
        self.setLayout(self.main_box)

        self.top_row = QtGui.QHBoxLayout()
        self.main_box.addLayout(self.top_row)

        self.combobox = QtGui.QComboBox()
        self.combobox.insertItem(0, 'P')
        self.combobox.insertItem(1, 'S')
        self.top_row.addWidget(QtGui.QLabel('Select a phase: '))        
        self.top_row.addWidget(self.combobox)
        
        fig = plt.figure()
        self.main_ax = fig.add_subplot(111)
        self.canvas = fig.canvas
        self.main_box.addWidget(self.canvas)

        self.toolbar = NavigationToolbar(self.canvas, self)
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
        phase = self.combobox.currentText()
        def get_picks(station_names):
            picks=[]
            for station in station_names:
                try:
                    picks.append(self.pyl_mainwindow.autopicks[station][phase]['mpp'])
                except:
                    picks.append(np.nan)
            return picks

        def get_picks_rel(picks):
            picks_rel=[]
            minp = min(picks)    
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
        self.sc_picked = self.basemap.scatter(self.lon_no_nan, self.lat_no_nan, s=50,
                                              c=self.picks_no_nan, latlon=True, zorder=11, label='Picked')

    def annotate_ax(self):
        self.annotations=[]
        for index, name in enumerate(self.station_names):
            self.annotations.append(self.main_ax.annotate(' %s' % name, xy=(self.x[index], self.y[index]),
                                                          fontsize='x-small', zorder=12))
        self.legend=self.main_ax.legend()

    def add_cbar(self, label):
        cbar = self.main_ax.figure.colorbar(self.sc_picked)
        cbar.set_label(label)
        return cbar

    def refresh_drawings(self):
        self.remove_drawings()
        self.draw_everything()
    
    def draw_everything(self):
        self.init_picks()
        self.init_picks_active()
        self.init_stations_active()
        self.init_picksgrid()
        self.draw_contour_filled()
        self.scatter_all_stations()
        self.scatter_picked_stations()
        self.annotate_ax()
        self.cbar = self.add_cbar(label='Time relative to first onset [s]')
        self.canvas.draw()

    def remove_drawings(self):
        self.sc_picked.remove()
        self.sc.remove()
        self.cbar.remove()
        self.remove_annotations()
        self.legend.remove()
        self.remove_contourf()
        self.canvas.mpl_disconnect(self.cid)
        self.canvas.draw()

    def remove_contourf(self):
        for item in self.contourf.collections:
            item.remove()

    def remove_annotations(self):
        for annotation in self.annotations:
            annotation.remove()
    
    
