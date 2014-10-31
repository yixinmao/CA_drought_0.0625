#!/usr/local/bin/python

from netCDF4 import Dataset,date2num,num2date
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

######################## load latlon data ###########################
CDEC_snotel = np.loadtxt('/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_std_latlon/latlon_CDEC_snotel.txt')
CDEC_snow_course = np.loadtxt('/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_std_latlon/latlon_CDEC_snow_course.txt')
NRCS_snotel = np.loadtxt('/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_std_latlon/latlon_NRCS_snotel.txt')
NRCS_snow_course = np.loadtxt('/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_std_latlon/latlon_NRCS_snow_course.txt')

######################## plot stations #############################
# CDEC stations
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

x, y = m(-CDEC_snotel[:,1], CDEC_snotel[:,0])
m.scatter(x, y, s=10, c='b', marker='o', label='SNOTEL')
x, y = m(CDEC_snow_course[:,1], CDEC_snow_course[:,0])
m.scatter(x, y, s=10, c='r', marker='o', label='Snow course')
plt.legend()
plt.text(0.5, 1.08, 'CDEC snow stations', horizontalalignment='center', \
         fontsize=20, transform = ax.transAxes)
fig.savefig('/raid2/ymao/other/CA_drought_0.0625/snow_observation/output/CDEC_std.png', format='png')

# NRCS stations
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

x, y = m(NRCS_snotel[:,1], NRCS_snotel[:,0])
m.scatter(x, y, s=10, c='b', marker='o', label='SNOTEL')
x, y = m(NRCS_snow_course[:,1], NRCS_snow_course[:,0])
m.scatter(x, y, s=10, c='r', marker='o', label='Snow course')
plt.legend()
plt.text(0.5, 1.08, 'NRCS snow stations', horizontalalignment='center', \
         fontsize=20, transform = ax.transAxes)
fig.savefig('/raid2/ymao/other/CA_drought_0.0625/snow_observation/output/NRCS_std.png', format='png')




