#!/usr/local/bin/python

import numpy as np
import datetime as dt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib as mpl
import matplotlib.pyplot as plt

######################## user defined ############################
ori_vic_forcing_path = '/home/raid2/ymao/other/CA_drought_0.0625/run_vic/vic_forcing_1920-20140430'
latlonlist_path = '/raid2/ymao/other/CA_drought_0.0625/input/latlon.smaller'
output_dir = '/raid2/ymao/other/CA_drought_0.0625/output'

start_date = dt.datetime(year=1920, month=1, day=1)
end_date = dt.datetime(year=2014, month=4, day=30)

nday = (end_date-start_date).days + 1
nyear = end_date.year-start_date.year + 1
start_year = start_date.year
end_year = end_date.year

latlonlist = np.loadtxt(latlonlist_path)
nfile = np.shape(latlonlist)[0]

#nfile = 20

####################################################################################
################ calculate trend and adjust T for summer and/or winter #############
################ winter: Nov-Mar; summer: Apr-Oct ##################################
####################################################################################
trend_Tmax_winter = np.empty(nfile)
trend_Tmin_winter = np.empty(nfile)
trend_Tmax_summer = np.empty(nfile)
trend_Tmin_summer = np.empty(nfile)
############# calculate trend and adjust T ###################
for i in range(nfile):
	print 'Grid cell %d' %(i+1)
	################ read in original forcing file
	ori_data = np.loadtxt('%s/data_%.5f_%.5f' %(ori_vic_forcing_path,latlonlist[i,0],latlonlist[i,1]))

	################ winter #################
	# calculate annual average T
	Tmax_year = np.zeros(nyear)
	Tmin_year = np.zeros(nyear)
	nday_year = np.zeros(nyear)
	for t in range(nday):
		date = start_date + dt.timedelta(days=t)
		year = date.year
		month = date.month
		day = date.day
		Tmax = ori_data[t,1]
		Tmin = ori_data[t,2]
		if month>=11: # if Nov-Dec, add to the next water year
			Tmax_year[year-start_year+1] = Tmax_year[year-start_year+1] + Tmax 
			Tmin_year[year-start_year+1] = Tmin_year[year-start_year+1] + Tmin 
			nday_year[year-start_year+1] = nday_year[year-start_year+1] + 1
		elif month<=3:  # if Jan-Mar, add to this water year
			Tmax_year[year-start_year] = Tmax_year[year-start_year] + Tmax 
			Tmin_year[year-start_year] = Tmin_year[year-start_year] + Tmin 
			nday_year[year-start_year] = nday_year[year-start_year] + 1
	Tmax_year = Tmax_year / nday_year
	Tmax_year = (np.array([np.arange(start_year+1,end_year+1), Tmax_year[1:nyear]])).T  # only analyze 1921-2014 for winter [year; T]
	Tmin_year = Tmin_year / nday_year
	Tmin_year = (np.array([np.arange(start_year+1,end_year+1), Tmin_year[1:nyear]])).T  # only analyze 1921-2014 for winter [year; T]
	# calculate trend using linear regression
	x = Tmax_year[:,0]
	A = np.array([x, np.ones(np.shape(x)[0])])
	y = Tmax_year[:,1]
	w = np.linalg.lstsq(A.T, y)[0]
	trend_Tmax_winter[i] = w[0]  # Tmax; degC/year
	x = Tmin_year[:,0]
	A = np.array([x, np.ones(np.shape(x)[0])])
	y = Tmin_year[:,1]
	w = np.linalg.lstsq(A.T, y)[0]
	trend_Tmin_winter[i] = w[0]  # Tmax; degC/year

	################ summer #################
	# calculate annual average T
	Tmax_year = np.zeros(nyear)
	Tmin_year = np.zeros(nyear)
	nday_year = np.zeros(nyear)
	for t in range(nday):
		date = start_date + dt.timedelta(days=t)
		year = date.year
		month = date.month
		day = date.day
		Tmax = ori_data[t,1]
		Tmin = ori_data[t,2]
		if month>=4 and month<=10:  # if Apr-Oct, add to this water year
			Tmax_year[year-start_year] = Tmax_year[year-start_year] + Tmax 
			Tmin_year[year-start_year] = Tmin_year[year-start_year] + Tmin 
			nday_year[year-start_year] = nday_year[year-start_year] + 1
	Tmax_year = Tmax_year / nday_year
	Tmax_year = (np.array([np.arange(start_year+1,end_year+1), Tmax_year[1:nyear]])).T  # only analyze 1921-2014 for winter [year; T]
	Tmin_year = Tmin_year / nday_year
	Tmin_year = (np.array([np.arange(start_year+1,end_year+1), Tmin_year[1:nyear]])).T  # only analyze 1921-2014 for winter [year; T]
	# calculate trend using linear regression
	x = Tmax_year[:,0]
	A = np.array([x, np.ones(np.shape(x)[0])])
	y = Tmax_year[:,1]
	w = np.linalg.lstsq(A.T, y)[0]
	trend_Tmax_summer[i] = w[0]  # Tmax; degC/year
	x = Tmin_year[:,0]
	A = np.array([x, np.ones(np.shape(x)[0])])
	y = Tmin_year[:,1]
	w = np.linalg.lstsq(A.T, y)[0]
	trend_Tmin_summer[i] = w[0]  # Tmax; degC/year

################################### plot map ##############################
## winter, Tmax
#fig = plt.figure(figsize=(8,8))
#ax = plt.axes([0, 0.08, 1, 0.75])
#m = Basemap(projection='mill', llcrnrlat=32, urcrnrlat=44,\
#            llcrnrlon=-126, urcrnrlon=-112, resolution='l')
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[True,True,False,False], fontsize=16)
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[False,False,True,True], fontsize=16)
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#x, y = m(latlonlist[:,1], latlonlist[:,0])
#data = trend_Tmax_winter * (nyear-1)  # T change [deg C]
#cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu_r', vmax=2, vmin=-2, marker='s', linewidth=0)
#cbar = plt.colorbar(cs)
#cbar.set_label('T change (mm)')
#plt.text(0.5, 1.1, "Winter Tmax change over 1920-2014", \
#             horizontalalignment='center', \
#             fontsize=16, transform = ax.transAxes)
#fig.savefig('%s/trend_map_winter_Tmax.png' %output_dir, format='png')
#
## winter, Tmin
#fig = plt.figure(figsize=(8,8))
#ax = plt.axes([0, 0.08, 1, 0.75])
#m = Basemap(projection='mill', llcrnrlat=32, urcrnrlat=44,\
#            llcrnrlon=-126, urcrnrlon=-112, resolution='l')
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[True,True,False,False], fontsize=16)
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[False,False,True,True], fontsize=16)
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#x, y = m(latlonlist[:,1], latlonlist[:,0])
#data = trend_Tmin_winter * (nyear-1)  # T change [deg C]
#cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu_r', vmax=2, vmin=-2, marker='s', linewidth=0)
#cbar = plt.colorbar(cs)
#cbar.set_label('T change (mm)')
#plt.text(0.5, 1.1, "Winter Tmin change over 1920-2014", \
#             horizontalalignment='center', \
#             fontsize=16, transform = ax.transAxes)
#fig.savefig('%s/trend_map_winter_Tmin.png' %output_dir, format='png')

# summer, Tmax
fig = plt.figure(figsize=(8,8))
ax = plt.axes([0, 0.08, 1, 0.75])
m = Basemap(projection='mill', llcrnrlat=32, urcrnrlat=44,\
            llcrnrlon=-126, urcrnrlon=-112, resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 5.), labels=[True,True,False,False], fontsize=16)
m.drawmeridians(np.arange(-180., 181., 5.), labels=[False,False,True,True], fontsize=16)
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

x, y = m(latlonlist[:,1], latlonlist[:,0])
data = trend_Tmax_summer * (nyear-1)  # T change [deg C]
cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu_r', vmax=2, vmin=-2, marker='s', linewidth=0)
cbar = plt.colorbar(cs)
cbar.set_label('T change (mm)')
plt.text(0.5, 1.1, "Summer Tmax change over 1920-2014", \
             horizontalalignment='center', \
             fontsize=16, transform = ax.transAxes)
fig.savefig('%s/trend_map_summer_Tmax.png' %output_dir, format='png')

# summer, Tmin
fig = plt.figure(figsize=(8,8))
ax = plt.axes([0, 0.08, 1, 0.75])
m = Basemap(projection='mill', llcrnrlat=32, urcrnrlat=44,\
            llcrnrlon=-126, urcrnrlon=-112, resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 5.), labels=[True,True,False,False], fontsize=16)
m.drawmeridians(np.arange(-180., 181., 5.), labels=[False,False,True,True], fontsize=16)
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

x, y = m(latlonlist[:,1], latlonlist[:,0])
data = trend_Tmin_summer * (nyear-1)  # T change [deg C]
cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu_r', vmax=2, vmin=-2, marker='s', linewidth=0)
cbar = plt.colorbar(cs)
cbar.set_label('T change (mm)')
plt.text(0.5, 1.1, "Summer Tmin change over 1920-2014", \
             horizontalalignment='center', \
             fontsize=16, transform = ax.transAxes)
fig.savefig('%s/trend_map_summer_Tmin.png' %output_dir, format='png')



