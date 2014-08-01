#!/usr/local/bin/python

from netCDF4 import Dataset,date2num,num2date
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

##############################
# user defined
##############################
infile1 = 'input/latlon.smaller.higher_third'
infile2 = 'input/latlon.smaller.middle_third'
infile3 = 'input/latlon.smaller.lower_third'
outfile = './output/grid_map_elevation_bands.png'

##############################
# plot map
##############################
latlon1 = np.loadtxt(infile1)
latlon2 = np.loadtxt(infile2)
latlon3 = np.loadtxt(infile3)

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

x, y = m(latlon1[:,1], latlon1[:,0])
m.scatter(x, y, s=10, c='b', marker='s', linewidths=0)

x, y = m(latlon2[:,1], latlon2[:,0])
m.scatter(x, y, s=10, c='r', marker='s', linewidths=0)

x, y = m(latlon3[:,1], latlon3[:,0])
m.scatter(x, y, s=10, c='g', marker='s', linewidths=0)

fig.savefig(outfile, format='png')
plt.show()




