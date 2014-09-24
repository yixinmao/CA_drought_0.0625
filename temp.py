#!/usr/local/bin/python

import numpy as np
import datetime as dt
import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt

######################## user defined ############################
latlonlist_path = '/raid2/ymao/other/CA_drought_0.0625/input/latlon.smaller'
old_forcing_dir = '/raid2/ymao/other/CA_drought_0.0625/run_vic/vic_forcing_detrend_T_P_pivot1920_1920-Jul2014'
new_forcing_dir = '/raid2/ymao/other/CA_drought_0.0625/run_vic/vic_forcing_detrend_T_P_pivot1920_1920-Jul2014_new'

latlonlist = np.loadtxt(latlonlist_path)
nfile = np.shape(latlonlist)[0]

nday = 34546

for i in range(nfile):
	print (i+1)
	old_data = np.loadtxt('%s/data_%.5f_%.5f' %(old_forcing_dir, latlonlist[i,0], latlonlist[i,1]))
	new_data = old_data
	for t in range(nday):
		if new_data[t,0]<0:
			new_data[t,0] = 0.0
	f = open('%s/data_%.5f_%.5f' %(new_forcing_dir, latlonlist[i,0], latlonlist[i,1]), 'w')
	for t in range(nday):
		f.write('%.2f %.2f %.2f %.1f\n' %(new_data[t,0], new_data[t,1], new_data[t,2], new_data[t,3]))
	f.close()


