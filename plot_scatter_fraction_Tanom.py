#!/usr/local/bin/python

###########################################
################ HEADER ###################
###########################################
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import pdb

##########################################
# user defined
##########################################
in_path = './output/Tanom.largest5_fraction'   # year Tanom swe
out_path = './output/fraction_Tanom.largest5.png'

#########################################
# load data and plot
#########################################
data = np.loadtxt(in_path)

fig = plt.figure(figsize=(8,8))
ax = plt.axes()
plt.plot(data[:,1], data[:,2], 'o')
plt.xlabel('Average T anomaly (Nov1-Mar23), $^o$C', fontsize=24)
plt.ylabel('Total SWE on Mar 23, km$^3$', fontsize=24)
for item in ([ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
	item.set_fontsize(20)
plt.xlim([-4,4])

fig.savefig(out_path, format='png')
plt.show()
