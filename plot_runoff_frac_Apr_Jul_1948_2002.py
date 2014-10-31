#!/usr/local/bin/python

#########################################################################
################## plot Apr-Jul fractional flow map #####################
#########################################################################

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
end_date = dt.datetime(year=2014, month=7, day=31)

duration = end_date - start_date
nday = duration.days + 1
nyear = (end_date.year - start_date.year) + 1
start_year = start_date.year
end_year = end_date.year

nyear_plot = 55

parser = argparse.ArgumentParser()
parser.add_argument("--latlonlist", help="Latlon list")
parser.add_argument("--ind", help="input VIC output files directory")
parser.add_argument("--outmap", help="output map file")
args = parser.parse_args()

latlonlist = np.loadtxt(args.latlonlist)
nfile = np.shape(latlonlist)[0]

############################# load data and calculate #######################
trend_runoff_frac_Apr_Jul = np.empty(nfile)
for i in range(nfile):
	print 'Grid cell %d' %(i+1)

	# load data for this grid cell
	filename = 'fluxes_%.5f_%.5f' %(latlonlist[i,0], latlonlist[i,1])
	data = np.loadtxt('%s/%s' %(args.ind, filename))  # year; month; day; prec; evap; runoff; baseflow; airT; sm1; sm2; sm3; swe

	# calculation
	runoff_whole_year = np.zeros(nyear)
	runoff_Apr_Jul = np.zeros(nyear)
	for t in range(nday):
		date = start_date + dt.timedelta(days=t)
		year = date.year
		month = date.month
		day = date.day
		year_ind = year - start_year  # year index; starts from 0
		runoff = data[t,5]
		baseflow = data[t,6]

		if month>=10:  # if Oct-Dec, add it to the next water year
			runoff_whole_year[year_ind+1] = runoff_whole_year[year_ind+1] + (runoff+baseflow)
		elif month<=9:  # if Jan-Sep, add it to this water year
			runoff_whole_year[year_ind] = runoff_whole_year[year_ind] + (runoff+baseflow)
		if month>=4 and month<=7:  # add it to this water year
			runoff_Apr_Jul[year_ind] = runoff_Apr_Jul[year_ind] + (runoff+baseflow)

	runoff_frac_Apr_Jul = runoff_Apr_Jul / runoff_whole_year

	# calculate linear trend of fractional flow for this grid cell (only consider 1948-2002 water years)
	x = range(1948, 2003)
	x = np.asarray(x).T
	A = np.array([x, np.ones(np.shape(x)[0])])
	y = runoff_frac_Apr_Jul[28:83]
	w = np.linalg.lstsq(A.T, y)[0]
	trend_runoff_frac_Apr_Jul[i] = w[0]  # year-1


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
data = trend_runoff_frac_Apr_Jul * nyear_plot * 100  # fractional flow change [%/(nyear-2)]
cs = plt.scatter(x[0:nfile], y[0:nfile], s=10, c=data, cmap='RdBu', vmax=20, vmin=-20, marker='s', linewidth=0)
cbar = plt.colorbar(cs)
cbar.set_label('Fractional flow change (%/55 years)')
plt.text(0.5, 1.1, "Apr-Jul fractional flow over 1948-2002", \
             horizontalalignment='center', \
             fontsize=16, transform = ax.transAxes)
fig.savefig(args.outmap, format='png')



