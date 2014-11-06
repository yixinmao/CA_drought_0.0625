#!/usr/local/bin/python

import numpy as np
import datetime as dt
import argparse
import pdb
import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt

############################# set parameters #############################
MAXY = 200
MAXC = 500
cellsize = 0.0625

start_date = dt.datetime(year=1920, month=1, day=1)
end_date = dt.datetime(year=2014, month=9, day=30)

duration = end_date - start_date
nday = duration.days + 1
nyear = (end_date.year - start_date.year) + 1
start_year = start_date.year
end_year = end_date.year

parser = argparse.ArgumentParser()
parser.add_argument("--latlonlist", help="Latlon list")
parser.add_argument("--arealist", help="area list [lat; lon; area]")
parser.add_argument("--ind", help="input VIC output files directory")
parser.add_argument("--outf", help="output file")
parser.add_argument("--outmap", help="output map file")
args = parser.parse_args()

latlonlist = np.loadtxt(args.latlonlist)
arealist = np.loadtxt(args.arealist)
nfile = np.shape(latlonlist)[0]

############################# load data and calculate #######################
tot_area = 0
runoff_CT_ave = np.zeros(nyear)
trend_CT = np.empty(nfile)
for i in range(nfile):
	print 'Grid cell %d' %(i+1)
	if arealist[i,0]!=latlonlist[i,0] or arealist[i,1]!=latlonlist[i,1]:
		print "Error: area list does not match with latlon list!"
		exit()
	area = arealist[i,2]
	tot_area = tot_area + area

	# load data for this grid cell
	filename = 'fluxes_%.5f_%.5f' %(latlonlist[i,0], latlonlist[i,1])
	data = np.loadtxt('%s/%s' %(args.ind, filename))  # year; month; day; prec; evap; runoff; baseflow; airT; sm1; sm2; sm3; swe

	# calculation
	runoff_CT_year = np.zeros(nyear)
	runoff_year = np.zeros(nyear)
	for t in range(nday):
		date = start_date + dt.timedelta(days=t)
		year = date.year
		month = date.month
		day = date.day
		year_ind = year - start_year  # year index; starts from 0
		runoff = data[t,5]
		baseflow = data[t,6]

		if month>=10:  # if Oct-Dec, add it to the next water year
			ti = (date - dt.datetime(year=year, month=11, day=1)).days + 1 # time in days from the beginning of the water year, start from 1
			runoff_CT_year[year_ind+1] = runoff_CT_year[year_ind+1] + (runoff+baseflow)*ti
			runoff_year[year_ind+1] = runoff_year[year_ind+1] + (runoff+baseflow)
		elif month<=9:  # if Jan-Sep, add it to this water year
			ti = (date - dt.datetime(year=year-1, month=10, day=1)).days + 1 # time in days from the beginning of the water year, start from 1
			runoff_CT_year[year_ind] = runoff_CT_year[year_ind] + (runoff+baseflow)*ti
			runoff_year[year_ind] = runoff_year[year_ind] + (runoff+baseflow)

	for y in range(nyear):
		runoff_CT_year[y] = runoff_CT_year[y] / runoff_year[y]  # unit: day
		runoff_CT_ave[y] = runoff_CT_ave[y] + runoff_CT_year[y] * area

	# calculate linear trend of CT for this grid cell
	x = range(start_year+1, end_year+1)
	x = np.asarray(x).T
	A = np.array([x, np.ones(np.shape(x)[0])])
	y = runoff_CT_year[1:nyear]
	w = np.linalg.lstsq(A.T, y)[0]
	trend_CT[i] = w[0]  # day/year

runoff_CT_ave = runoff_CT_ave / tot_area

############################## plot trend map ##################################
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
data = trend_CT * (nyear-2)  # T change [deg C]
cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu', vmax=10, vmin=-10, marker='s', linewidth=0)
cbar = plt.colorbar(cs)
cbar.set_label('CT change (day)')
plt.text(0.5, 1.1, "Runoff CT change over 1921-2014", \
             horizontalalignment='center', \
             fontsize=16, transform = ax.transAxes)
fig.savefig(args.outmap, format='png')

############################## print out results ################################
f = open(args.outf, 'w')
for y in range(1, nyear):  # only print out 1921-2014
	f.write("%4d %.4f\n" %(start_year+y, runoff_CT_ave[y]))
f.close()



