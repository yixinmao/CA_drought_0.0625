#!/usr/local/bin/python

import numpy as np
import datetime as dt
import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt

######################## user defined ############################
ori_vic_forcing_path = '/home/raid2/ymao/other/CA_drought_0.0625/run_vic/vic_forcing_1920-20140731'
latlonlist_path = '/raid2/ymao/other/CA_drought_0.0625/input/latlon.smaller'
output_map_dir = '/raid2/ymao/other/CA_drought_0.0625/output'
output_new_forcing_dir = '/raid2/ymao/other/CA_drought_0.0625/run_vic/vic_forcing_detrend_T_pivot1920_1920-Jul2014'

start_date = dt.datetime(year=1920, month=1, day=1)
end_date = dt.datetime(year=2014, month=7, day=31)

nday = (end_date-start_date).days + 1
nyear = end_date.year-start_date.year + 1
start_year = start_date.year
end_year = end_date.year

latlonlist = np.loadtxt(latlonlist_path)
nfile = np.shape(latlonlist)[0]

is_trend_winter_Tmax = 0  # 1 for having significant trend; 0 for not having significant trend
is_trend_winter_Tmin = 1
is_trend_summer_Tmax = 0
is_trend_summer_Tmin = 1

pivot_year = end_year


####################################################################################
######################## calculate trend  ##########################################
################ winter: Nov-Mar; summer: Apr-Oct ##################################
####################################################################################
trend_Tmax_winter = np.empty(nfile)
trend_Tmin_winter = np.empty(nfile)
trend_Tmax_summer = np.empty(nfile)
trend_Tmin_summer = np.empty(nfile)
############# calculate trend and adjust T ###################
for i in range(nfile):
	print 'Grid cell %d' %(i+1)
	######################################################################
	################ read in original forcing file #######################
	######################################################################
	ori_data = np.loadtxt('%s/data_%.5f_%.5f' %(ori_vic_forcing_path,latlonlist[i,0],latlonlist[i,1]))

	######################################################################
	##################### calculate trend ################################
	######################################################################

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

	#########################################################################
	############################ detrend T ##################################
	#########################################################################

	new_data = ori_data
	###################### detrend Tmin #######################
	# winter
	if is_trend_winter_Tmin==1: 
		for t in range(nday):
			date = start_date + dt.timedelta(days=t)
			year = date.year
			month = date.month
			day = date.day
			Tmin_ori = ori_data[t,2]
			if month>=11:  # if Nov-Dec, account for next water year
				Tmin_new = Tmin_ori + trend_Tmin_winter[i] * (pivot_year - (year+1))
				new_data[t,2] = Tmin_new
			elif month<=3:  # if Jan-Mar, account for this water year
				Tmin_new = Tmin_ori + trend_Tmin_winter[i] * (pivot_year - year)
				new_data[t,2] = Tmin_new
	# summer
	if is_trend_summer_Tmin==1: 
		for t in range(nday):
			date = start_date + dt.timedelta(days=t)
			year = date.year
			month = date.month
			day = date.day
			Tmin_ori = ori_data[t,2]
			if month>=4 and month<=10:  # if Apr-Oct, account for this water year
				Tmin_new = Tmin_ori + trend_Tmin_summer[i] * (pivot_year - year)
				new_data[t,2] = Tmin_new

	###################### detrend Tmax #######################
	# winter
	if is_trend_winter_Tmax==1: 
		for t in range(nday):
			date = start_date + dt.timedelta(days=t)
			year = date.year
			month = date.month
			day = date.day
			Tmax_ori = ori_data[t,1]
			if month>=11:  # if Nov-Dec, account for next water year
				Tmax_new = Tmax_ori + trend_Tmax_winter[i] * (pivot_year - (year+1))
				new_data[t,1] = Tmax_new
			elif month<=3:  # if Jan-Mar, account for this water year
				Tmax_new = Tmax_ori + trend_Tmax_winter[i] * (pivot_year - year)
				new_data[t,1] = Tmax_new
	# summer
	if is_trend_summer_Tmax==1: 
		for t in range(nday):
			date = start_date + dt.timedelta(days=t)
			year = date.year
			month = date.month
			day = date.day
			Tmax_ori = ori_data[t,1]
			if month>=4 and month<=10:  # if Apr-Oct, account for this water year
				Tmax_new = Tmax_ori + trend_Tmax_summer[i] * (pivot_year - year)
				new_data[t,1] = Tmax_new

	################################################################################
	########################### write new data into files ##########################
	f = open('%s/data_%.5f_%.5f' %(output_new_forcing_dir,latlonlist[i,0],latlonlist[i,1]), 'w')
	for t in range(nday):
		f.write('%.2f %.2f %.2f %.1f\n' %(new_data[t,0],new_data[t,1],new_data[t,2],new_data[t,3]))
	f.close()



############################################################################
############################## plot trend map ##############################
############################################################################
# winter, Tmax
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
data = trend_Tmax_winter * (nyear-1)  # T change [deg C]
cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu_r', vmax=2, vmin=-2, marker='s', linewidth=0)
cbar = plt.colorbar(cs)
cbar.set_label('T change (mm)')
plt.text(0.5, 1.1, "Winter Tmax change over 1920-2014", \
             horizontalalignment='center', \
             fontsize=16, transform = ax.transAxes)
fig.savefig('%s/trend_map_winter_Tmax.png' %output_map_dir, format='png')

# winter, Tmin
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
data = trend_Tmin_winter * (nyear-1)  # T change [deg C]
cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu_r', vmax=2, vmin=-2, marker='s', linewidth=0)
cbar = plt.colorbar(cs)
cbar.set_label('T change (mm)')
plt.text(0.5, 1.1, "Winter Tmin change over 1920-2014", \
             horizontalalignment='center', \
             fontsize=16, transform = ax.transAxes)
fig.savefig('%s/trend_map_winter_Tmin.png' %output_map_dir, format='png')

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
fig.savefig('%s/trend_map_summer_Tmax.png' %output_map_dir, format='png')

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
fig.savefig('%s/trend_map_summer_Tmin.png' %output_map_dir, format='png')



