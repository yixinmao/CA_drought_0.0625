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
vic_output_dir = '/raid2/ymao/other/CA_drought_0.0625/run_vic/model_output/fluxes_1920_20140930_old2014_snowband' # year; month; day; prec; evap; runoff; baseflow; airT; sm1; sm2; sm3; swe; shortwave; swe for each snow band (5 columns)
snow_obs_list_path = '/raid2/ymao/other/CA_drought_0.0625/snow_observation/snow_data_from_Darrin/SWE2014/CA_DEC/CATXTFiles.txt'
latlon_list_path = '/raid2/ymao/other/CA_drought_0.0625/input/latlon.smaller'
cellnum_list_path = '/raid2/ymao/other/CA_drought_0.0625/input/cellnum_latlon_list' # cellnum; lat; lon; should be the same order with latlon_list_path
snowband_path = '/raid2/ymao/other/CA_drought_0.0625/run_vic/snow.band.smaller' # snow band file; should be the same order with cellnum_list_path and latlon_list_path
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

# load cellnum list
cellnum_list = np.loadtxt(cellnum_list_path)

# load snowband
snowband = np.loadtxt(snowband_path)

# load site list
snow_obs_list = []
f = open(snow_obs_list_path, 'r')
while 1:
	line = f.readline().rstrip("\n")
	if line=="":
		break
	snow_obs_list.append(line)
f.close()
nobs_orig = len(snow_obs_list)  # number of all sites


###################################################################
################### load snow obs data ############################
###################################################################
print 'Loading snow obs data...'
data_obs = [] # [[site name, abbr, elev(ft), lat, lon], 
              # np.array(year, month, SWE(inch))]
for i in range(nobs_orig):
	data_obs.append([])
	f = open('%s/%s' %(snow_obs_dir, snow_obs_list[i]))
	data_obs[i].append([])
	data_obs[i][0].append(f.readline().rstrip("\n")) # site name
	data_obs[i][0].append(f.readline().rstrip("\n")) # site abbr
	line = f.readline().rstrip("\n") 
	data_obs[i][0].append(float(line.split()[0])) # elev (ft)
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
             # np.array from vic output (year, month, day, SWE for the snowband with the closest elev(mm), ariT(degC), airT (degC)]
             #        (if grid cell not in the study domain, delete it)
print 'Loading VIC output...'
elev_diff = []
for i in range(nobs_orig):
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
		cell_ind = -1
		for j in range(len(cellnum_list)):
			if cellnum_list[j,1]==lat_grid and cellnum_list[j,2]==lon_grid:
				cell_ind = j  # index starts from 0
				break
		if cellnum_list[cell_ind,1]!=latlon_list[cell_ind,0] or cellnum_list[cell_ind,2]!=latlon_list[cell_ind,1] or cellnum_list[cell_ind,0]!=snowband[cell_ind,0]:  # make sure we are reading the correct snowband grid cell
			print 'Error in matching snowband file!'
			exit()
		else:
			# Read vic SWE output at the snowband with the closest elevation as the site
			area_fraction = snowband[cell_ind,1:6]
			elev = snowband[cell_ind,6:11] # unit: m
			elev[np.absolute(area_fraction)<0.001] = np.nan # set elev with area_fraction=0 to nan
			elev_site = data_obs[i][0][2] * 12.0*25.4/1000.0 # convert unit to m
			band_ind = np.nanargmin(np.absolute(elev_site-elev)) # closest band index, starting from 0
			print 'Elev_site - elev_band = %.1f' %(elev_site-elev[band_ind]), band_ind
			elev_diff.append(elev_site-elev[band_ind])
			data_obs[i].append(np.loadtxt('%s/fluxes_%.5f_%.5f' %(vic_output_dir, lat_grid, lon_grid), usecols=(0,1,2,13+band_ind, 7)))  # snowband swe: column 13-17 (ind starts from 0)

for i in range(nobs_orig-1, -1, -1):
	if type(data_obs[i][3])==int:
		del data_obs[i]
nobs = len(data_obs)  # number of sites falling in our study domain

# plot difference of elevation between snowband and site
fig = plt.figure()
plt.hist(elev_diff)
plt.xlabel('Elevation difference (Obs.-VIC (snowband))', fontsize=16)
plt.ylabel('Number of sites', fontsize=16)
fig.savefig('%s/hist_elev_diff.png' %output_dir, format='png')


###################################################################
################# check grid cells with multiple sites ############
###################################################################

###################################################################
####### calculate average Apr 1 SWE in conrresponding years #######
###################################################################
data_avg = [] # [[site name, abbr, elev(ft), lat, lon], 
              # mean obs Apr 1 SWE (mm)
              # [lat_grid, lon_grid],
              # mean Apr 1 SWE (mm) from VIC output in corresponding years
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
	swe_vic_avg = 0
	for year in years:
		month = 4
		day = 1
		date_ind = (dt.datetime(year=year, month=month, day=day)-start_date).days # index starts from 0
		swe = data_obs[i][3][date_ind, 3]
		swe_vic_avg = swe_vic_avg + swe
	swe_vic_avg = swe_vic_avg / len(years)
	data_avg[i].append(swe_vic_avg)

	# calculate mean Nov-Mar T from VIC in corresponding years
	T_winter_avg = 0
	count = 0
	for year in years:
		date1 = dt.datetime(year=year-1, month=11, day=1)
		date2 = dt.datetime(year=year, month=3, day=31)
		day_count = (date2-date1).days + 1
		for date in (date1+dt.timedelta(n) for n in range(day_count)): # loop over all days in Nov-Mar in this wate year
			date_ind = (date-start_date).days # index starts from 0
			T = data_obs[i][3][date_ind, 4]
			T_winter_avg = T_winter_avg + T
			count = count + 1
	T_winter_avg = T_winter_avg/count
	data_avg[i].append(T_winter_avg)

###################################################################
###################### making scatter plot  #######################
###################################################################
fig = plt.figure()
for i in range(nobs):
	plt.plot(data_avg[i][1], data_avg[i][3], 'bo')
plt.plot(np.arange(0,2000), np.arange(0,2000), 'k--')
plt.xlabel('Observed Apr 1 SWE (mm)', fontsize=16)
plt.ylabel('VIC-simulated Apr 1 SWE (mm)', fontsize=16)
plt.xlim(0,2000)
plt.ylim(0,2000)
plt.axis('scaled')
fig.savefig('%s/snow_obs_scatter_snowband.png' %output_dir, format='png')

###################################################################
####### making binned plot (binned by every 1 K winter T) #########
###################################################################
# put SWE data into T bins
T_winter_avg = np.asarray([row[4] for row in data_avg]) # deg C
swe_obs_avg = np.asarray([row[1] for row in data_avg]) # mm
swe_vic_avg = np.asarray([row[3] for row in data_avg]) # mm
T_winter_avg_bin = np.around(T_winter_avg).astype(int) # rounded T
T_min_bin = np.min(T_winter_avg_bin)
T_max_bin = np.max(T_winter_avg_bin)
nbin = T_max_bin - T_min_bin + 1
swe_obs_binned = []
swe_vic_binned = []
# binned data structure: [[swe1, swe2, swe3, ...] -> Tmin
#                         [swe1, swe2, swe3, ...] -> Tmin+1]
#                         ...
#                         [swe1, swe2, ...] > Tmax ]
for i in range(nbin):
	swe_obs_binned.append([])
	swe_vic_binned.append([])
for i in range(nobs):
	Tind = T_winter_avg_bin[i] - T_min_bin  # index starts from 0
	swe_obs_binned[Tind].append(swe_obs_avg[i])
	swe_vic_binned[Tind].append(swe_vic_avg[i])

# calculate statistics (for each bin, a dict: [n (number of data)] [min] [max] [10th] [90th] [median])
bin_swe_stat_obs = []
bin_swe_stat_vic = []
for i in range(nbin):
	# obs
	bin_swe_stat_obs.append({})
	bin_swe_stat_obs[i]['n'] = len(swe_obs_binned[i])  # number of data
	if bin_swe_stat_obs[i]['n']>=1:
		bin_swe_stat_obs[i]['min'] = np.min(swe_obs_binned[i]) # min
		bin_swe_stat_obs[i]['max'] = np.max(swe_obs_binned[i]) # max
		bin_swe_stat_obs[i]['10th'] = np.percentile(swe_obs_binned[i], 10) # 10th
		bin_swe_stat_obs[i]['90th'] = np.percentile(swe_obs_binned[i], 90) # 90th
		bin_swe_stat_obs[i]['median'] = np.median(swe_obs_binned[i]) # median
	# vic
	bin_swe_stat_vic.append({})
	bin_swe_stat_vic[i]['n'] = len(swe_vic_binned[i])  # number of data
	if bin_swe_stat_vic[i]['n']>=1:
		bin_swe_stat_vic[i]['min'] = np.min(swe_vic_binned[i]) # min
		bin_swe_stat_vic[i]['max'] = np.max(swe_vic_binned[i]) # max
		bin_swe_stat_vic[i]['10th'] = np.percentile(swe_vic_binned[i], 10) # 10th
		bin_swe_stat_vic[i]['90th'] = np.percentile(swe_vic_binned[i], 90) # 90th
		bin_swe_stat_vic[i]['median'] = np.median(swe_vic_binned[i]) # median

fig = plt.figure()
for i in range(nbin):
	# vic 
	if bin_swe_stat_vic[i]['n']>=1:
		plt.plot([bin_swe_stat_vic[i]['min'], bin_swe_stat_vic[i]['max']], [T_min_bin+i,T_min_bin+i], 'r--')
		if bin_swe_stat_vic[i]['n']>=10:
			plt.plot([bin_swe_stat_vic[i]['10th'], bin_swe_stat_vic[i]['90th']], [T_min_bin+i,T_min_bin+i], 'r-')
			plt.plot([bin_swe_stat_vic[i]['10th'], bin_swe_stat_vic[i]['90th']], [T_min_bin+i,T_min_bin+i], 'r+')
		plt.plot(bin_swe_stat_vic[i]['median'], T_min_bin+i, 'ro')
	if i>=1 and bin_swe_stat_vic[i]['n']>=1 and bin_swe_stat_vic[i-1]['n']>=1: # connect median values
		plt.plot([bin_swe_stat_vic[i-1]['median'], bin_swe_stat_vic[i]['median']], [T_min_bin+i-1,T_min_bin+i], 'r-')
	# obs
	if bin_swe_stat_obs[i]['n']>=1:
		plt.plot([bin_swe_stat_obs[i]['min'], bin_swe_stat_obs[i]['max']], [T_min_bin+i-0.2,T_min_bin+i-0.2], 'b--')
		if bin_swe_stat_obs[i]['n']>=10:
			plt.plot([bin_swe_stat_obs[i]['10th'], bin_swe_stat_obs[i]['90th']], [T_min_bin+i-0.2,T_min_bin+i-0.2], 'b-')
			plt.plot([bin_swe_stat_obs[i]['10th'], bin_swe_stat_obs[i]['90th']], [T_min_bin+i-0.2,T_min_bin+i-0.2], 'b+')
		plt.plot(bin_swe_stat_obs[i]['median'], T_min_bin+i-0.2, 'bo')
	if i>=1 and bin_swe_stat_obs[i]['n']>=1 and bin_swe_stat_obs[i-1]['n']>=1: # connect median values
		plt.plot([bin_swe_stat_obs[i-1]['median'], bin_swe_stat_obs[i]['median']], [T_min_bin+i-0.2-1,T_min_bin-0.2+i], 'b-')

plt.plot(bin_swe_stat_vic[0]['median'], T_min_bin, 'ro', label='VIC simulated')
plt.plot(bin_swe_stat_obs[0]['median'], T_min_bin-0.2, 'bo', label='Observed')

plt.legend(loc='upper right', fontsize=16)
plt.ylim(T_min_bin-1, T_max_bin+1)
plt.xlabel('Apr 1 SWE (mm)', fontsize=16)
plt.ylabel('Nov-Mar winter temperature ($^0$C)', fontsize=16)
fig.savefig('%s/swe_obs_binned_snowband.png' %output_dir, format='png')





