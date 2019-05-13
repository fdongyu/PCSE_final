from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import utm
import numpy as np

import pdb


def convert_to_latlon(x, y):
    x = x.flatten()
    y = y.flatten()
    nn = len(x)
    lat = np.zeros_like(x)
    lon = np.zeros_like(y)
    
    for i in range(len(x)):
        lat[i], lon[i] = utm.to_latlon(x[i], y[i], 15, 'U')

    return lat, lon



west = -95.218; east = -94.412
south = 28.979; north = 29.789


filename='Output/Euler/xy_Euler.nc'
nc = Dataset(filename, 'r')
print nc
time = nc.variables['time'][:]
x = nc.variables['X'][:]
y = nc.variables['Y'][:]

filename2='Output/RK4/xy_RK4.nc'
nc2 = Dataset(filename2, 'r')
print nc2
time2 = nc2.variables['time'][:]
x2 = nc2.variables['X'][:]
y2 = nc2.variables['Y'][:]


lat, lon = convert_to_latlon(x, y)
lat2, lon2 = convert_to_latlon(x2, y2)

lat_tem, lon_tem = convert_to_latlon(x[:,0], y[:,0])

fig = plt.figure(figsize=(10,9.5))
ax = fig.add_subplot(111)
basemap = Basemap(projection='merc',llcrnrlat=south,urcrnrlat=north,\
                llcrnrlon=west,urcrnrlon=east, resolution='h')
        
basemap.drawcoastlines()
basemap.fillcontinents()
basemap.drawcountries()
basemap.drawstates()

llons, llats=basemap(*(lon,lat))
llons2, llats2=basemap(*(lon2,lat2))
llons_tem, llats_tem=basemap(*(lon_tem,lat_tem))

basemap.plot(llons,llats, '.', color='r', markersize=4)    # Euler
basemap.plot(llons2,llats2, '.', color='b', markersize=4)  # RK4

basemap.plot(llons_tem,llats_tem, '.', color='y', markersize=8)

## add coordination
lats = np.linspace(south, north,4)
lons = np.linspace(west, east,4)
lonsnew, latsnew = basemap(lons, lats)
    
sw = utm.from_latlon(south, west)[:2]
ne = utm.from_latlon(north, east)[:2]
xlims = np.asarray([sw[0], ne[0]])
ylims = np.asarray([sw[1], ne[1]])
originx = 327500
originy = 3244000
xs = (np.linspace(xlims[0], xlims[1], 4) - originx)/1000.
ys = (np.linspace(ylims[0], ylims[1], 4) - originy)/1000.
xlabel = []
ylabel = []
for i in range(len(xs)):
    xlabel.append(str(round(xs[i],1)))
    ylabel.append(str(round(ys[i],1)))
    
ax.set_xticks((lonsnew))
ax.set_yticks((latsnew))
ax.set_xticklabels(xlabel, fontsize=22)
ax.set_yticklabels(ylabel, fontsize=22)
ax.set_aspect('equal')
ax.set_xlabel('Easting (km)', fontsize=22)
ax.set_ylabel('Northing (km)', fontsize=22)

basemap.plot(llons[0],llats[0], '.', color='r', markersize=8, label='Euler')    # Euler
basemap.plot(llons2[0],llats2[0], '.', color='b', markersize=8, label='RK4')  # RK4
ax.legend(fontsize=20,numpoints=1)
ax.grid()
plt.tight_layout()

#plt.plot(x,y, 'ob',markersize=2.5)
#plt.plot(x2,y2, 'or',markersize=2.0)
#plt.plot(x[:,0],y[:,0], 'or', markersize=5)
#plt.show()
#plt.savefig('figures/Eular.png')
#plt.savefig('figures/RK4.png')
plt.savefig('comparison_basemap_new.png')
plt.close()
#pdb.set_trace()
