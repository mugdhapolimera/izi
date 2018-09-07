import numpy as np
import os
from astropy.table import Table
import pandas as pd

# Setting System Path
import sys
sys.path.append('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/izi_utils/')
#sys.path.append('C:\Users\mugdhapolimera\github\izi')

#Importing Custom Utility Files
from izi import izi


#Set Path to Data Directory
os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/')
#os.chdir('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes')

#Load Resolve Catalog 
#TODO : Load from the server
inputfile = 'RESOLVE_SDSS_dext.fits'
dat = Table.read(inputfile, format='fits')
infile = dat.to_pandas()

cols = infile.keys()


#Define Names of Flux lines according to Catalog
#TODO: Make automated script   
fluxnames  = ['oii_3726_flux_ext', 'oii_3729_flux_ext', 'h_gamma_flux_ext', 'h_beta_flux_ext', 'oiii_4959_flux_ext', 
             'oiii_5007_flux_ext', 'oi_6300_flux_ext', 'nii_6548_flux_ext', 'h_alpha_flux_ext', 'nii_6584_flux_ext', 
             'sii_6717_flux_ext', 'sii_6731_flux_ext'] 
errornames = ['oii_3726_flux_ext_err', 'oii_3729_flux_ext_err', 'h_gamma_flux_ext_err', 'h_beta_flux_ext_err', 
              'oiii_4959_flux_ext_err', 'oiii_5007_flux_ext_err', 'oi_6300_flux_ext_err', 'nii_6548_flux_ext_err', 
              'h_alpha_flux_ext_err', 'nii_6584_flux_ext_err', 'sii_6717_flux_ext_err', 'sii_6731_flux_ext_err']
   

idin = ['oii3726', 'oii3729', 'hgamma', 'hbeta', 'oiii4959', 'oiii5007', 'oi6300', 'nii6548', 'halpha', 'nii6584', 
        'sii6717', 'sii6731']


'''
#For the old data table
fluxnames = ['F_OIII_5007_ERR_BROAD' ,
'F_OIII_5007_BROAD' ,
'F_OIII_4960_ERR_BROAD' ,
'F_OIII_4960_BROAD' ,
'F_OIII_4363_ERR_BROAD' ,
'F_OIII_4363_BROAD' ,
'F_NII_6586_ERR_BROAD' ,
'F_NII_6586_BROAD' ,
'F_NII_6548_ERR_BROAD' ,
'F_NII_6548_BROAD' ,
'F_HB_ERR_BROAD' ,
'F_HB_BROAD' ,
'F_HA_ERR_BROAD' ,
'F_HA_BROAD']
'''
#Making Flux and Error Arrays
fluxin = []
errorin = []
for i in range(0,len(fluxnames)):
        fluxin.append(infile[fluxnames[i]][0])
        errorin.append(infile[errornames[i]][0])

#Defining IDs for different Flux Lines

#idin = ['oiii5007', 'oiii4959', 'oiii4363' , 'nii6584',  'nii6548', 'hbeta', 'halpha' ]
#kwargs =  logOHsun, intergridfile, logzlimits, logqlimits ,logzprior, logqprior (logz/q prior have no application in the original code)
izi(fluxin, errorin, idin)
