from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import obspy
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from PySide import QtCore, QtGui

#import QtPyLoT
from pylot.core.util.dataprocessing import read_metadata
from pylot.core.util.widgets import PickDlg


#pf=QtPyLoT.main()

def onpick(event):
    ind = event.ind
    if ind == []:
        return
    data = pf.get_data().getWFData()
    for index in ind:
        station=str(station_names[index])
        try:
            pickDlg = PickDlg(pf, infile=pf.getinfile(), 
                              data=data.select(station=station),
                              station=station,
                              picks=pf.getPicksOnStation(station, 'manual'),
                              autopicks=pf.getPicksOnStation(station, 'auto'))
            pickDlg.exec_()
        except Exception as e:
            print('Could not generate Plot for station {st}.\n{er}'.format(st=station, er=e))
        

        
def get_metadata(path):
    metadata=read_metadata(path)
    parser=metadata[1]
    return parser

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

def get_picks(pf, station_names):
    picks=[]
    for station in station_names:
        try:
            picks.append(pf.autopicks[station]['P']['mpp'])
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

def remove_nan_picks(picks):
    picks_no_nan=[]
    for pick in picks:
        if not np.isnan(pick):
            picks_no_nan.append(pick)
    return picks_no_nan
            
def remove_nan_lat_lon(picks, lat, lon):
    lat_no_nan=[]
    lon_no_nan=[]
    for index, pick in enumerate(picks):
        if not np.isnan(pick):
            lat_no_nan.append(lat[index])
            lon_no_nan.append(lon[index])
    return lat_no_nan, lon_no_nan
    
def get_lon_lat_dim(lon, lat):
    londim = max(lon) - min(lon)
    latdim = max(lat) - min(lat)
    return londim, latdim

def get_x_y_dim(x, y):
    xdim = max(x) - min(x)
    ydim = max(y) - min(y)
    return xdim, ydim
    
def init_map(projection, ax, resolution='l'):
    m = Basemap(projection=projection, resolution = resolution, ax=ax)
    m.drawmapboundary(fill_color='darkblue')
    m.drawcountries()
    m.drawstates()
    m.fillcontinents(color='grey', lake_color='aqua')
    m.drawcoastlines()
    return m
    
def get_lat_lon_axis(lat, lon):
    steplat = (max(lat)-min(lat))/250
    steplon = (max(lon)-min(lon))/250
    
    lataxis = np.arange(min(lat), max(lat), steplat)
    lonaxis = np.arange(min(lon), max(lon), steplon)
    return lataxis, lonaxis

def get_lat_lon_grid(lataxis, lonaxis):
    longrid, latgrid = np.meshgrid(lonaxis, lataxis)
    return latgrid, longrid

def draw_contour_filled(picks, longrid, latgrid, picksgrid, levels='50'):
    levels = np.linspace(min(picks), max(picks), 50)
    contourf = m.contourf(longrid, latgrid, picksgrid, levels, latlon=True, zorder=9)
    return contourf

def annotate_ax(ax, x, y, station_names):
    for index, name in enumerate(station_names):
        ax.annotate(' %s' % name, xy=(x[index], y[index]), fontsize='x-small', zorder=12)

def connect_pick(ax, onpick):
    ax.figure.canvas.mpl_connect('pick_event', onpick)

def add_cbar(ax, scatter, label):
    cbar = ax.figure.colorbar(scatter)
    cbar.set_label(label)
    return cbar

qwidget = QtGui.QWidget()
QtGui.QVBoxLayout()
main_box = QtGui.QVBoxLayout()
qwidget.setLayout(main_box)

fig = plt.figure()
main_ax = fig.add_subplot(111)
canvas = main_ax.figure.canvas
main_box.addWidget(canvas)

toolbar = NavigationToolbar(canvas, qwidget)
main_box.addWidget(toolbar)


parser = get_metadata('/data/Geothermie/Insheim/STAT_INFO/MAGS2_net.dless')

station_names, lat, lon = get_station_names_lat_lon(parser)

picks = get_picks(pf, station_names)
picks_rel = get_picks_rel(picks)

picks_no_nan = remove_nan_picks(picks_rel)
lat_no_nan, lon_no_nan = remove_nan_lat_lon(picks_rel, lat, lon)

londim, latdim = get_lon_lat_dim(lon, lat)
x, y = m(lon, lat)
xdim, ydim = get_x_y_dim(x, y)

m = init_map(projection='mill', ax=main_ax, resolution='l')

lataxis, lonaxis = get_lat_lon_axis(lat, lon)
latgrid, longrid = get_lat_lon_grid(lataxis, lonaxis)

picksgrid_no_nan = griddata((lat_no_nan, lon_no_nan), picks_no_nan, (latgrid, longrid), method='linear')

contourf = draw_contour_filled(picks_no_nan, longrid, latgrid, picksgrid_no_nan)

sc = m.scatter(lon, lat, s=50, facecolor='none', latlon=True, zorder=10, picker=True, edgecolor='m', label='Not Picked')
sc_picked = m.scatter(lon_no_nan, lat_no_nan, s=50, c=picks_no_nan, latlon=True, zorder=11, label='Picked')

annotate_ax(main_ax, x, y, station_names)

main_ax.legend()

connect_pick(main_ax, onpick)

cbar = add_cbar(main_ax, sc_picked, label='Time relative to first onset [s]')

canvas.draw()

qwidget.show()

# ax.set_xlim(min(x)-0.5*xdim, max(x)+0.5*xdim)
# ax.set_ylim(min(y)-0.5*ydim, max(y)+0.5*ydim)



#plt.show()
