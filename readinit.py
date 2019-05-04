import numpy as np
import matplotlib.pyplot as plt

import pdb

def readdat(filename):
    #filename = 'trajectory_360.dat'
    output = np.loadtxt(filename)
    x = output[:,0]
    y = output[:,1]
    return x, y

x1, y1 = readdat('DATA/init.dat')
#x2, y2 = readdat('tracks/trajectory_3600.dat')
#x3, y3 = readdat('tracks/trajectory_720.dat')

plt.plot(x1, y1, '.r', markersize=8)
#plt.plot(x2, y2, '.b', markersize=10)
#plt.plot(x3, y3, '.k', markersize=5)
plt.show()

pdb.set_trace()
