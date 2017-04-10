from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import obspy

#import QtPyLoT
from pylot.core.util.dataprocessing import read_metadata

#pf=QtPyLoT.main()

def onpick(event):
    ind = event.ind
    print(ind)


metadata=read_metadata('/data/Geothermie/Insheim/STAT_INFO/MAGS2_net.dless')
parser=metadata[1]

lat=[]
lon=[]
stcl=[]
picks=[]
picks_rel=[]

for station in parser.stations:
    stcl.append(station[0].station_call_letters)
    lat.append(station[0].latitude)
    lon.append(station[0].longitude)

for station in stcl:
    try:
        picks.append(pf.autopicks[station]['P']['mpp'])
    except:
        picks.append(np.nan)

minp = min(picks)

for pick in picks:
    if type(pick) is obspy.core.utcdatetime.UTCDateTime:
        pick -= minp
    picks_rel.append(pick)


m = Basemap(projection='hammer',lon_0=0,lat_0=0)
m.drawmapboundary(fill_color='aqua')
m.drawcountries()
m.drawstates()
m.fillcontinents(color='grey', lake_color='aqua')
m.drawcoastlines()
m.scatter(lon, lat, s=50, c=picks_rel, latlon=True, zorder=10, picker=True)

ax = plt.gca() # IMPROVE!!!!

ax.figure.canvas.mpl_connect('pick_event', onpick)

plt.show()
