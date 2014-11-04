#!/usr/local/bin/python

import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import os.path

###################################################################
########################### user defined ##########################
###################################################################
snow_obs_dir = '/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_data_from_Darrin/SWE2014/CA_DEC'
# 1st: name; 2st: abbr; 3nd: [elev(ft), lat, lon(no negative sign)]; 4th: Course Year Month SWE Meas.Date Name; from 5th: data 
vic_output_dir = '/raid2/ymao/other/CA_drought_0.0625/run_vic/model_output/fluxes_1920_20140930_old2014' # year; month; day; prec; evap; runoff; baseflow; airT; sm1; sm2; sm3; swe; shortwave
snow_obs_list_path = '/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_data_from_Darrin/SWE2014/CA_DEC/CATXTFiles.txt'
latlon_list_path = '/raid2/ymao/other/CA_drought_0.0625/input/latlon.smaller'
output_dir = '/raid2/ymao/other/CA_drought_0.0625/snow_observation/output'

start_date = dt.datetime(year=1920, month=1, day=1)
end_date = dt.datetime(year=2014, month=9, day=30)

duration = end_date - start_date
nday = duration.days + 1
nyear = (end_date.year - start_date.year) + 1
start_year = start_date.year
end_year = end_date.year

# load latlon list
latlon_list = np.loadtxt(latlon_list_path)
ncell = len(latlon_list)

# load site list
snow_obs_list = []
f = open(snow_obs_list_path, 'r')
while 1:
	line = f.readline().rstrip("\n")
	if line=="":
		break
	snow_obs_list.append(line)
f.close()
nobs = len(snow_obs_list)

###################################################################
################### load snow obs data ############################
###################################################################
print 'Loading snow obs data...'
data_obs = [] # [[site name, abbr, elev(ft), lat, lon], 
              # np.array(year, month, SWE(inch))]
for i in range(nobs):
	data_obs.append([])
	f = open('%s/%s' %(snow_obs_dir, snow_obs_list[i]))
	data_obs[i].append([])
	data_obs[i][0].append(f.readline().rstrip("\n")) # site name
	data_obs[i][0].append(f.readline().rstrip("\n")) # site abbr
	line = f.readline().rstrip("\n") 
	data_obs[i][0].append(line.split()[0]) # elev (ft)
	data_obs[i][0].append(float(line.split()[1])) # lat
	data_obs[i][0].append(-float(line.split()[2])) # lon
	f.close()
	data_obs[i].append([])
	data_obs[i][1] = np.loadtxt('%s/%s' %(snow_obs_dir, snow_obs_list[i]), skiprows=4, usecols=(1,2,3))

###################################################################
# find grid cell each site falls in, and load corresponding data ##
###################################################################
# data_obs:  # [[site name, abbr, elev(ft), lat, lon], 
             # np.array(year, month, SWE(inch)),
             # [lat_grid, lon_grid],
             # np.array from vic output (year, month, day, SWE(mm), ariT(degC))]
             #        (if grid cell not in the study domain, -1)
print 'Loading VIC output...'
for i in range(nobs):
	print 'Loading grid cell %d...' %(i+1)
	# find grid cell the site falls in
	lat = data_obs[i][0][3]
	lon = data_obs[i][0][4]
	lat_grid = round(16*lat-0.5)/16.0 + 0.03125
	lon_grid = round(16*lon-0.5)/16.0 + 0.03125
	data_obs[i].append([lat_grid, lon_grid])
	# load SWE data from VIC output, only the grid corresponding grid cell
	filepath = '%s/fluxes_%.5f_%.5f' %(vic_output_dir, lat_grid, lon_grid)
	if os.path.isfile(filepath)==False: # if grid cell not in the study domain, append -1
		data_obs[i].append(-1)
	else:
		data_obs[i].append(np.loadtxt('%s/fluxes_%.5f_%.5f' %(vic_output_dir, lat_grid, lon_grid), usecols=(0,1,2,11,7)))

###################################################################
################# check grid cells with multiple sites ############
###################################################################

###################################################################
####### calculate average Apr 1 SWE in conrresponding years #######
###################################################################
data_avg = [] # [[site name, abbr, elev(ft), lat, lon], 
             # mean obs Apr 1 SWE (mm)
             # [lat_grid, lon_grid],
             # mean Apr 1 SWE (mm) in corresponding years
             #        (if grid cell not in the study domain, -1)
             # mean Nov-Mar T (degC) in corresponding years

for i in range(nobs):
	data_avg.append([])
	data_avg[i].append(data_obs[i][0])

	# calculate mean obs Apr 1 SWE, and record measured years
	years = []
	swe_obs_avg = 0
	for j in range(len(data_obs[i][1])):
		if data_obs[i][1][j,1]==4:  # if the mesurement is on (or around) Apr 1
			years.append(int(data_obs[i][1][j,0]))
			swe_obs_avg = swe_obs_avg + data_obs[i][1][j,2] * 25.4 # convert unit to mm
	swe_obs_avg = swe_obs_avg / len(years)
	data_avg[i].append(swe_obs_avg)

	data_avg[i].append(data_obs[i][2])
	# calculate mean Apr 1 SWE from VIC in corresponding years
	if type(data_obs[i][3])==int: # if this grid cell not available
		data_avg[i].append(-1)
	else:
		swe_vic_avg = 0
		for year in years:
			month = 4
			day = 1
			date_ind = (dt.datetime(year=year, month=month, day=day)-start_date).days # index starts from 0
			swe = data_obs[i][3][date_ind, 3]
			swe_vic_avg = swe_vic_avg + swe
		swe_vic_avg = swe_vic_avg / len(years)
		data_avg[i].append(swe_vic_avg)

###################################################################
###################### making scatter plot  #######################
###################################################################
fig = plt.figure()
for i in range(nobs):
	if data_avg[i][3]!=-1:
		plt.plot(data_avg[i][1], data_avg[i][3], 'bo')
plt.plot(np.arange(0,2000), np.arange(0,2000), 'k--')
plt.xlabel('Observed Apr 1 SWE (mm)', fontsize=16)
plt.ylabel('VIC-simulated Apr 1 SWE (mm)', fontsize=16)
plt.xlim(0,2000)
plt.ylim(0,2000)
plt.axis('scaled')
fig.savefig('%s/snow_obs_scatter.png' %output_dir, format='png')

###################################################################
###################### making binned plot  ########################
###################################################################
