#!/usr/local/bin/python

from netCDF4 import Dataset,date2num,num2date
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

##############################
# user defined
##############################
infile = 'input/latlon.smaller'
outfile = './output/grid_map.png'

##############################
# plot map
##############################
latlon = np.loadtxt(infile)

fig = plt.figure(figsize=(8,8))
ax = plt.axes([0, 0.08, 1, 0.75])
m = Basemap(projection='mill', llcrnrlat=32, urcrnrlat=44,\
            llcrnrlon=-126, urcrnrlon=-112, resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 5.), labels=[True,True,False,False], fontsize=20)
m.drawmeridians(np.arange(-180., 181., 5.), labels=[False,False,True,True], fontsize=20)
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

x, y = m(latlon[:,1], latlon[:,0])

m.scatter(x, y, s=10, c='b', marker='s', linewidths=0)

fig.savefig(outfile, format='png')
plt.show()




