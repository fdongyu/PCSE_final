from netCDF4 import Dataset
import matplotlib.pyplot as plt

import pdb


filename='DATA/xy_Euler.nc'
nc = Dataset(filename, 'r')
print nc
time = nc.variables['time'][:]
x = nc.variables['X'][:]
y = nc.variables['Y'][:]


pdb.set_trace()
#filename2='src_RK4/DATA/xy_RK4.nc'
#nc2 = Dataset(filename2, 'r')
#print nc2
#time2 = nc2.variables['time'][:]
#x2 = nc2.variables['X'][:]
#y2 = nc2.variables['Y'][:]

plt.plot(x,y, 'ob',markersize=2.5)
#plt.plot(x2,y2, 'or',markersize=2.0)
plt.plot(x[:,0],y[:,0], 'or', markersize=5)
plt.show()
#plt.savefig('figures/Eular.png')
#plt.savefig('figures/RK4.png')
#plt.savefig('figures/comparison.png')
plt.close()
#pdb.set_trace()
