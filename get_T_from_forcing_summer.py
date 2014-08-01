#!/usr/local/bin/python

import numpy as np
import datetime as dt
import argparse
import pdb

############################# set parameters #############################
MAXY = 200
MAXC = 500
cellsize = 0.0625

start_date = dt.datetime(year=1920, month=1, day=1)
end_date = dt.datetime(year=2014, month=4, day=30)

duration = end_date - start_date
nday = duration.days + 1
nyear = (end_date.year - start_date.year) + 1
start_year = start_date.year

parser = argparse.ArgumentParser()
parser.add_argument("--latlonlist", help="Latlon list")
parser.add_argument("--arealist", help="area list [lat; lon; area]")
parser.add_argument("--ind", help="input VIC forcing files directory")
parser.add_argument("--outfTmax", help="output file")
parser.add_argument("--outfTmin", help="output file")
args = parser.parse_args()

latlonlist = np.loadtxt(args.latlonlist)
arealist = np.loadtxt(args.arealist)
nfiles = np.shape(latlonlist)[0]

############################# load data and calculate #######################
tot_area = 0
Tmax_year = np.zeros(nyear)
Tmin_year = np.zeros(nyear)
for i in range(nfiles):
	print 'Grid cell %d' %(i+1)
	if arealist[i,0]!=latlonlist[i,0] or arealist[i,1]!=latlonlist[i,1]:
		print "Error: area list does not match with latlon list!"
		exit()
	area = arealist[i,2]
	tot_area = tot_area + area

	# load data for this grid cell
	filename = 'data_%.5f_%.5f' %(latlonlist[i,0], latlonlist[i,1])
	data = np.loadtxt('%s/%s' %(args.ind, filename))  # prec[mm/d]; Tmax [deg C]; Tmin [deg C]; wind speed [m/s]

	# calculation
	nday_year = np.zeros(nyear)
	for t in range(nday):
		date = start_date + dt.timedelta(days=t)
		year = date.year
		month = date.month
		day = date.day
		year_ind = year - start_year  # year index; starts from 0
		Tmax = data[t,1]
		Tmin = data[t,2]

		if month>=4 and month<=10:  # if Apr-Oct, add it to this water year
			Tmax_year[year_ind] = Tmax_year[year_ind] + area*Tmax
			Tmin_year[year_ind] = Tmin_year[year_ind] + area*Tmin
			nday_year[year_ind] = nday_year[year_ind] + 1

for y in range(nyear):
	Tmax_year[y] = Tmax_year[y] / nday_year[y] / tot_area
	Tmin_year[y] = Tmin_year[y] / nday_year[y] / tot_area

print nday_year

############################## print out results ################################
f = open(args.outfTmax, 'w')
for y in range(nyear-1):  # only print out 1920-2013
	f.write("%4d %.4f\n" %(start_year+y, Tmax_year[y]))
f.close()
f = open(args.outfTmin, 'w')
for y in range(nyear-1):  # only print out 1920-2013
	f.write("%4d %.4f\n" %(start_year+y, Tmin_year[y]))
f.close()




