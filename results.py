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
#idl = np.genfromtxt('RESOLVE_izioutv1.txt', dtype = None)
py = np.genfromtxt('IZI_Z_idl_Amandalines.txt', dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])
idl  = np.genfromtxt('IZI_Z_idl_Amanda.txt', dtype = None, names = ['name', 'Z', 'err_up', 'err_down'])
#os.chdir('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/')
gpy = np.genfromtxt('IZI_Z_Amanda.txt', dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])

#np.polyfit(idl['Z'],py['Z'],1)#,w=1/np.sqrt(y_err**2 + x_err**2),full=False,cov=True)
xarr = np.arange(7.3, 9.3, 0.1)
def line(p, xarr):
    return p[0]*xarr + p[1]

plot_flag = 1
if plot_flag:    
    plt.figure()
    plt.plot(np.arange(7.3,9.3), np.arange(7.3,9.3), 'b-')
    plt.xlabel("IDL Z estimates")
    plt.ylabel("Python Z estimates (with scipy interpolation)")
    plt.errorbar(idl['Z'], py['Z'], xerr = [idl['err_down'], idl['err_up']], yerr = [py['err_down'], py['err_up']], fmt = 'o')
    p = np.polyfit(idl['Z'],py['Z'],1)
    plt.plot(xarr, line(p,xarr))    
    #plt.errorbar(idl['Z'], gpy['Z'], xerr = [idl['err_down'], idl['err_up']], yerr = [gpy['err_down'], gpy['err_up']], fmt = 'ro')
    
    
    plt.figure()
    plt.plot(np.arange(7.3,9.3), np.arange(7.3,9.3), 'b-')
    plt.xlabel("IDL Z estimates")
    plt.ylabel("Python Z estimates (with gpy interpolation)")
    plt.errorbar(idl['Z'], gpy['Z'], xerr = [idl['err_down'], idl['err_up']], yerr = [gpy['err_down'], gpy['err_up']], fmt = 'ro')    
    p = np.polyfit(idl['Z'],gpy['Z'],1)
    plt.plot(xarr, line(p,xarr))    
    #plt.errorbar(py['Z'], gpy['Z'], xerr = [py['err_down'], py['err_up']], yerr = [gpy['err_down'], gpy['err_up']], fmt = 'ro')
    
    plt.figure()
    plt.plot(np.arange(7.3,9.3), np.zeros(len(np.arange(7.3,9.3))), 'b-')
    plt.errorbar(py['Z'], py['Z']-idl['Z'],fmt = 'ro')
    plt.xlabel('Python + scipy estimate')
    plt.ylabel('Residuals ( python[Z] - idl[Z])')
    #p = np.polyfit(idl['Z'],py['Z'],1)
    #plt.plot(xarr, line(p,xarr))    
    '''    
    
    py_matched = [i for i in range(len(idl)) if (abs(py['Z'][i]+ - idl['Z'][i]) 
         <= intersection((-py['err_down'][i], abs(py['err_up'][i])), (-idl['err_down'][i], idl['err_up'][i])))]
    correct = len(py_matched) * 100.0 /len(py['Z'])
    plt.title("%4.2f%% same predictions" %correct)

    plt.figure()
    plt.plot(np.arange(7.3,9.3), np.zeros(len(np.arange(7.3,9.3))), 'b-')
    plt.errorbar(gpy['Z'], gpy['Z']-idl['Z'],fmt = 'ro')
    plt.xlabel('Python + gpy estimate')
    plt.ylabel('Residuals ( python[Z] - idl[Z])')
    gpy_matched = [i for i in range(len(idl)) if (abs(gpy['Z'][i]+ - idl['Z'][i]) 
         <= intersection((-gpy['err_down'][i], abs(gpy['err_up'][i])), (-idl['err_down'][i], idl['err_up'][i])))]
    correct = len(gpy_matched) * 100.0 /len(gpy['Z'])
    plt.title("%4.2f%% same predictions" %correct)
'''
#unmatched = np.where(abs(py['Z'] - idl['Z']) >= 10**-4)[0]
#Some galaxies have error_up as negative!
unmatched = [i for i in range(len(idl)) if (abs(py['Z'][i]+ - idl['Z'][i]) 
         >= intersection((-py['err_down'][i], abs(py['err_up'][i])), (-idl['err_down'][i], idl['err_up'][i])))]

#for i in unmatched:
#    print idl['name'][i], idl['Z'][i], py['Z'][i], idl['err_up'][i], py['err_up'][i]

#Galaxies with either negative or 0 error estimates
weird = [i for i in range(len(idl)) if idl['err_up'][i]<=0. or py['err_up'][i]<=0.]

if sys.platform == 'linux2':
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results')

else:
    os.chdir('C:/Users/mugdhapolimera/github/izi/results/')

#python_data = pickle.load(open("IZI_results.pkl", "r"))
'''f = open("IZI_results.pkl", "r")
for i in range(11):#len(idl)):
    d = pickle.load(f)
    if (i ==2) or (i in weird):
        print d['name']
        print d.flux[np.where(d.flux != -666.0)]        
'''