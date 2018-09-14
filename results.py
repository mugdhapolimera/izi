#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 11:04:24 2018

@author: mugpol
"""

import numpy as np
import os
import matplotlib.pyplot as plt
#os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results')
os.chdir('C:/Users/mugdhapolimera/github/izi/results/')
#idl = np.genfromtxt('RESOLVE_izioutv1.txt', dtype = None)
idl = np.genfromtxt('IZI_Z_idl.txt', dtype = None, names = ['name', 'Z', 'err_up', 'err_down'])
py  = np.genfromtxt('IZI_Z.txt', dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])
os.chdir('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/')
gpy = np.genfromtxt('IZI_Z2_gpy.txt', dtype = None, names = ['name', 'Z', 'err_down', 'err_up'])

print py['Z'].shape
np.polyfit(idl['Z'],py['Z'],1)#,w=1/np.sqrt(y_err**2 + x_err**2),full=False,cov=True)

plt.figure()
plt.plot(np.arange(7,10), np.arange(7,10), 'b-')
plt.xlabel("IDL Z estimates")
plt.ylabel("Python Z estimates (with scipy interpolation)")
#for i in range(len(py)):
#    plt.errorbar(idl[i][1], py[i][1], xerr = [[idl[i][3]], [idl[i][2]]], yerr = [[py[i][2]], [py[i][3]]], fmt = 'ro')
plt.errorbar(idl['Z'], py['Z'], xerr = [idl['err_down'], idl['err_up']], yerr = [py['err_down'], py['err_up']], fmt = 'o')

plt.errorbar(idl['Z'], gpy['Z'], xerr = [idl['err_down'], idl['err_up']], yerr = [gpy['err_down'], gpy['err_up']], fmt = 'ro')

plt.figure()
plt.plot(np.arange(7,10), np.arange(7,10), 'b-')
plt.errorbar(py['Z'], gpy['Z'], xerr = [py['err_down'], py['err_up']], yerr = [gpy['err_down'], gpy['err_up']], fmt = 'ro')
