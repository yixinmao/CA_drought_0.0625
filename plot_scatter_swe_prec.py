#!/usr/local/bin/python

###########################################
################ HEADER ###################
###########################################
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import pdb

##########################################
# user defined
##########################################
in_path = './output/prec_swe'   # year prec swe
out_path = './output/scatter_swe_prec.png'

#########################################
# load data and plot
#########################################
data = np.loadtxt(in_path)

fig = plt.figure(figsize=(8,8))
ax = plt.axes()
plt.plot(data[:,1], data[:,2], 'o')
plt.xlabel('Accumulative precipitation (Nov1-Apr1), km$^3$', fontsize=24)
plt.ylabel('Total SWE on Apr 1, km$^3$', fontsize=24)
for item in ([ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
	item.set_fontsize(20)
plt.xlim([0, 80])

fig.savefig(out_path, format='png')
plt.show()
