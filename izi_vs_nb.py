#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 11:04:24 2018

@author: mugpol
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import pickle
import pandas as pd

if sys.platform == 'linux2':
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results')

else:
    os.chdir('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/')
from astropy.table import Table
inputfile = 'RESOLVE_SDSS_dext.fits'
dat = Table.read(inputfile, format='fits')
infile = dat.to_pandas()

#filter out bad data and merge gooddata with blue E/S0's
gooddata = ((infile['oiii_5007_flux_ext'] > 0) & 
(infile['h_alpha_flux_ext'] > 0 ) & (infile['h_beta_flux_ext'] > 0 ) & (infile['oi_6300_flux_ext'] > 0 ) & (infile['sii_6717_flux_ext'] + infile['sii_6731_flux_ext'] >0 ) & (infile['Flux_HeII_4685_ext'] > 0))

goodnames = list(infile["NAME"][gooddata])

if sys.platform == 'linux2':
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/results')

else:
    os.chdir('C:/Users/mugdhapolimera/github/izi/results/')

def intersection(interval_1, interval_2):
    start = max(interval_1[0], interval_2[0])
    end = min(interval_1[1], interval_2[1])
    if start < end:
        return end - start
    return None
idl = np.genfromtxt('IZI_Z_idl.txt', dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])
py  = np.genfromtxt("IZI_Z.txt", dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])
idl_good = [i for i in range(len(idl)) if idl["name"][i] in goodnames]
idl = idl[idl_good]
py_good = [i for i in range(len(py)) if py["name"][i] in goodnames]
py = py[py_good]
#os.chdir('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/')
#gpy = np.genfromtxt('IZI_Z2_gpy.txt', dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])

nb = pd.read_csv(r"C:\Anaconda2\Lib\site-packages\NebulaBayes\docs\results\RESOLVE_param_estimates.csv")

nb_Z_ind = nb['Parameter'] == 'LOGZ'
nb_Z = nb[nb_Z_ind]
nb_Z = nb_Z.set_index(np.arange(len(nb_Z)))
nb_good = [i for i in range(len(nb_Z)) if nb_Z["Galaxy Name"][i] in goodnames]
nb_Z = nb_Z.loc[nb_good]
nb_Z = nb_Z.set_index(np.arange(len(nb_Z)))

xarr = np.arange(7.3, 9.3, 0.1)
def line(p, xarr):
    return p[0]*xarr + p[1]
#print list(nb_Z['Galaxy Name'])
unmatched_py = []
unmatched_nb = []
matched_nb = []
matched_py = []

plot_flag = 0
if plot_flag:    
    fig1, ax1 = plt.subplots()
    ax1.plot(np.arange(7.3,9.3), np.arange(7.3,9.3), 'b-')
    ax1.set_xlabel("NebulaBayes Z estimates")
    ax1.set_ylabel("Python Z estimates (with scipy interpolation)")

    fig2, ax2 = plt.subplots()
    ax2.plot(np.arange(7.3,9.3), np.arange(7.3,9.3), 'b-')
    ax2.set_xlabel("NebulaBayes Z estimates")
    ax2.set_ylabel("IDL Z estimates")

    for i in range (len(nb_Z)):
        if nb_Z['Galaxy Name'][i] in py['name'] :
                j  = np.where(py['name'] == nb_Z["Galaxy Name"][i])[0][0]
                x_down = nb_Z['Estimate'][i] - nb_Z['CI68_low'][i]
                x_up = nb_Z['CI68_high'][i] - nb_Z['Estimate'][i]
                ax1.errorbar(nb_Z['Estimate'][i], py['Z'][j], 
                             xerr = [[x_down], [x_up]], yerr = [[py['err_down'][j]], [py['err_up'][j]]], fmt = 'bo')
                if (abs(py['Z'][j]+ - nb_Z['Estimate'][i]) >= intersection((-py['err_down'][i], abs(py['err_up'][i])), (-x_down, x_up))):
                    unmatched_py.append(j)
                    unmatched_nb.append(i)
                else:
                    matched_py.append(j)
                    matched_nb.append(i)

                ax2.errorbar(nb_Z['Estimate'][i], idl['Z'][j], 
                             xerr = [[x_down], [x_up]], yerr = [[idl['err_down'][j]], [idl['err_up'][j]]], fmt = 'bo')

else:
    for i in range (len(nb_Z)):
        if nb_Z['Galaxy Name'][i] in py['name'] :       
                j  = np.where(py['name'] == nb_Z["Galaxy Name"][i])[0][0]
                x_down = nb_Z['Estimate'][i] - nb_Z['CI68_low'][i]
                x_up = nb_Z['CI68_high'][i] - nb_Z['Estimate'][i]
                if (abs(py['Z'][j]+ - nb_Z['Estimate'][i]) >= intersection((-py['err_down'][i], abs(py['err_up'][i])), (-x_down, x_up))):
                    unmatched_py.append(j)
                    unmatched_nb.append(i)
                else:
                    matched_py.append(j)
                    matched_nb.append(i)
         
    #p = np.polyfit(idl['Z'],py['Z'],1)
    #plt.plot(xarr, line(p,xarr))    
    #plt.errorbar(idl['Z'], gpy['Z'], xerr = [idl['err_down'], idl['err_up']], yerr = [gpy['err_down'], gpy['err_up']], fmt = 'ro')
count = 0
#Some galaxies have error_up as negative!
for i,j in zip(unmatched_nb, unmatched_py):
    if py['err_up'][j] >0:
        count += 1        
        #print py['name'][j], nb_Z['Estimate'][i], py['Z'][i], nb_Z['CI68_high'][i] - nb_Z['Estimate'][i], py['err_up'][j]
print count

#Galaxies with either negative or 0 error estimates
weird = [i for i in range(len(py)) if py['err_up'][i]<=0.] #idl['err_up'][i]<=0. or 


#python_data = pickle.load(open("IZI_results.pkl", "r"))
os.chdir('C:/Users/mugdhapolimera/github/BPT/')
comp = pd.read_csv("Composite.csv", header = None, names = ["Galaxy Name", "Composite Flag"])
comp_weird = [comp["Galaxy Name"][i] for i in range(len(comp)) if comp["Composite Flag"][i] and comp["Galaxy Name"][i] in py["name"][weird]]
print comp_weird
print len(np.where(comp["Composite Flag"])[0]), len(comp_weird)

if sys.platform == 'linux2':
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results')

else:
    os.chdir('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/')

for i in range(len(nb_Z)):
    if nb_Z["Galaxy Name"][i] in comp_weird:
        x_down = nb_Z['Estimate'][i] - nb_Z['CI68_low'][i]
        x_up = nb_Z['CI68_high'][i] - nb_Z['Estimate'][i]
        print nb_Z["Galaxy Name"][i], nb_Z["Estimate"][i], x_down, x_up
for i in range(len(py)):
    if py["name"][i] in comp_weird:
        print py[i]

f = open("IZI_results2 - Copy.pkl", "r")
count = 0
while (f is not EOFError):
    d = pickle.load(f)
    if not count:
        print d.id[np.where(d.error!= -666.0)]        
    count += 1
    if (d['name'] in comp_weird):
        print d['name'] #, d.flux[np.where(d.id == 'hbeta')[0]]
        print d.flux[np.where(d.flux!= -666.0)]                        
        print d.error[np.where(d.error!= -666.0)]        



