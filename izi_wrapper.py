import numpy as np
import os
from astropy.table import Table
import pandas as pd
import pickle
import time

# Setting System Path
import sys
#sys.path.append('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/izi_utils/')
sys.path.append('C:\Users\mugdhapolimera\github\izi')

#Importing Custom Utility Files
from izi import izi


#Set Path to Data Directory
#os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/')
os.chdir('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes')

#Load Resolve Catalog 
#TODO : Load from the server
inputfile = 'RESOLVE_SDSS_dext.fits'
dat = Table.read(inputfile, format='fits')
infile = dat#.to_pandas()
#print infile
cols = infile.keys()


#Define Names of Flux lines according to Catalog
#TODO: Make automated script   
fluxnames = ['oii_3726_flux_ext', 'oii_3729_flux_ext', 'neiii_3869_flux_ext',
              'h_gamma_flux_ext', 'oiii_4363_flux_ext', 
             'h_beta_flux_ext', 'oiii_4959_flux_ext', 'oiii_5007_flux_ext',
             'hei_5876_flux_ext', 'oi_6300_flux_ext', 'nii_6548_flux_ext', 
             'h_alpha_flux_ext', 'nii_6584_flux_ext', 'sii_6717_flux_ext', 
             'sii_6731_flux_ext', 'ariii_7135_flux_ext'] 

errornames = ['oii_3726_flux_ext_err', 'oii_3729_flux_ext_err', 
              'neiii_3869_flux_ext_err',  
              'h_gamma_flux_ext_err', 'oiii_4363_flux_ext_err', 
              'h_beta_flux_ext_err', 'oiii_4959_flux_ext_err', 
              'oiii_5007_flux_ext_err', 'hei_5876_flux_ext_err', 
              'oi_6300_flux_ext_err', 'nii_6548_flux_ext_err', 
              'h_alpha_flux_ext_err', 'nii_6584_flux_ext_err', 
              'sii_6717_flux_ext_err', 'sii_6731_flux_ext_err', 
              'ariii_7135_flux_ext_err']
   
#Defining IDs for different Flux Lines
#idin = ['oiii5007', 'oiii4959', 'oiii4363' , 'nii6584',  'nii6548', 'hbeta', 'halpha' ]

idin = ['oii3726', 'oii3729', 'neiii3869', 'oiii4363', 'hgamma', 
        'hbeta', 'oiii4959', 'oiii5007', 'hei5875', 'oi6300', 'nii6548', 
        'halpha', 'nii6584', 'sii6717', 'sii6731', 'ariii7136']

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
#f1 = open('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results/IZI_Z2.txt', 'w')
#f1 = open('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\results\IZI_Z2.txt', 'w+')
f1 = open('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/IZI_Z2.txt', 'a+')

f1.write("#Name \t Z_Estimate \t Err_down \t Err_up\r\n")
#f2 = open("/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results/IZI_q2.txt", "w")  
#f2 = open('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\results\IZI_q2.txt', 'w+')
f2 = open('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/IZI_q2.txt', 'a+')


f2.write("#Name \t q_Estimate \t Err_down \t Err_up\r\n")
#f3 = open("/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results/IZI_results2.pkl", "w")  
#f3 = open('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\results\IZI_results2.pkl', 'w+')
f3 = open('C:/Users/mugdhapolimera/Desktop/UNC/Courses/Research/Codes/results/IZI_results2.pkl', 'w+')

print len(infile['NAME'])
t1 = time.time()         
for gal in range(648,len(infile['NAME'])):
    fluxin = []
    errorin = []
    print gal, infile['NAME'][gal]
    for i in range(0,len(fluxnames)):
            fluxin.append(infile[fluxnames[i]][gal])
            errorin.append(infile[errornames[i]][gal])
    #kwargs =  logOHsun, intergridfile, logzlimits, logqlimits ,logzprior, logqprior (logz/q prior have no application in the original code)
    #d = izi(fluxin, errorin, idin, plot_flag = 0, name = str(infile['NAME'][gal]), intergridfile = '/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/outputgrid_linear.fits', interpolate_flag = False)
    d = izi(fluxin, errorin, idin, name = str(infile['NAME'][gal]), plot_flag = 0, print_flag = 0, method = 'gpy')
    #print d
    f1.write(" %s \t %f \t %f \t %f \r\n" %(d['name'], d.Zgrid, d.edownZgrid, d.eupZgrid))
    f2.write(" %s \t %f \t %f \t %f \r\n" %(d['name'], d.qgrid, d.edownqgrid, d.eupqgrid))
    pickle.dump(d, f3)    
    
f1.close()
f2.close()
f3.close()
t2 = time.time()
print t2-t1