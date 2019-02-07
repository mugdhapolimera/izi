# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 15:34:10 2019

@author: mugdhapolimera
"""
import pandas as pd 
from astropy.io import fits
from astropy.table import Table  # Used in converting to pandas DataFrame 

gridfile1 = r'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_STB99secular_csf_n1e2_5.0Myr_newgrid.fits'
#gridfile = 'C:/Users/mugdhapolimera/github/izi/richardson_agnfrac-0-1_STB99secular_csf_n1e2_5.0Myr.fits'
gridfile2 = r'C:/Users/mugdhapolimera/github/izi/l09_high_csf_n1e2_6.0Myr_new.fits'

gridfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\grids\l09_high_csf_n1e2_6.0Myr.fits'
grid = Table.read(gridfile, format='fits')

#grid = pd.read_csv(gridfile)
#gridfile = r'C:/Users/mugdhapolimera/github/izi/richardson_bpass_secular_csf_n1e2_10_0Myr.fits'
grid0 = Table.read(gridfile1, format='fits')
grid1 = grid0.to_pandas()
print grid1
grid1 = grid1[(grid1["AGNFRAC"] == 0) & (grid1["LOGQ"] < 8.0)]
grid1.index = range(len(grid1))
grid00 = Table.read(gridfile2, format='fits')
grid2 = grid00.to_pandas()

#for i in range(len(grid1)):
#    print grid1.LOGZ[i]+8.66, grid1.LOGQ[i], grid1.oii3726[i]

for i in range(len(grid2)):
    print grid2.LOGZ[i], grid2.LOGQ[i], grid2.oii3726[i]

for i in range(len(grid)):
    print grid['LOGZ'][i]+8.66, grid['LOGQ'][i], grid['FLUX'][i][6]

