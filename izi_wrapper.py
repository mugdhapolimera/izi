import numpy as np
import os
from astropy.table import Table
import pandas as pd
import pickle
import time
import matplotlib.pyplot as plt
from scipy.stats import norm

# Setting System Path
import sys
#Importing Custom Utility Files
#Set Path to Data Directory
if sys.platform == 'linux2':
    print 'Linux'
    sys.path.append('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/izi_utils/')
    from izi import izi
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/')

elif sys.platform == 'win32':
    print 'Windows'
    sys.path.append('C:\Users\mugdhapolimera\github\izi')
    from izi import izi
    #os.chdir('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes')

#Load Resolve Catalog 
#TODO : Load from the server
#inputfile = 'RESOLVE_SDSS_dext.fits'
#inputfile = 'RESOLVE_all2.fits'
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_filtered.pkl'
infile = pd.read_pickle(inputfile)

#dat = Table.read(inputfile, format='fits')
#infile = dat#.to_pandas()
#print infile
cols = infile.keys()
#infile.columns[0] = 'NAME'
#print infile.columns

#Define Names of Flux lines according to Catalog
#TODO: Make automated script   
'''fluxnames = ['oii_3726_flux_port_ext', 'oii_3729_flux_port_ext', 'neiii_3869_flux_ext',
              'h_gamma_flux_ext', 'oiii_4363_flux_ext', 
             'h_beta_flux_ext', 'oiii_4959_flux_ext', 'oiii_5007_flux_ext',
             'hei_5876_flux_ext', 'oi_6300_flux_ext', 'nii_6548_flux_ext', 
             'h_alpha_flux_ext', 'nii_6584_flux_ext', 'sii_6717_flux_ext', 
             'ariii_7135_flux_ext'] 
#'sii_6731_flux_ext', 
errornames = ['oii_3726_flux_port_ext_err', 'oii_3729_flux_port_ext_err', 
              'neiii_3869_flux_ext_err',  
              'h_gamma_flux_ext_err', 'oiii_4363_flux_ext_err', 
              'h_beta_flux_ext_err', 'oiii_4959_flux_ext_err', 
              'oiii_5007_flux_ext_err', 'hei_5876_flux_ext_err', 
              'oi_6300_flux_ext_err', 'nii_6548_flux_ext_err', 
              'h_alpha_flux_ext_err', 'nii_6584_flux_ext_err', 
              'sii_6731_flux_ext_err', 
              'ariii_7135_flux_ext_err']
#'sii_6717_flux_ext_err', 
'''   
fluxnames = ['oiii_4363_flux_ext', 
             'oiii_4959_flux_ext','h_beta_flux_ext', 'oiii_5007_flux_ext',
             'nii_6548_flux_ext', 'h_alpha_flux_ext', 'nii_6584_flux_ext', 
             'sii_6717_flux_ext'] 

errornames = ['oiii_4363_flux_ext_err', 
              'oiii_4959_flux_ext_err', 'h_beta_flux_ext_err',
              'oiii_5007_flux_ext_err', 'nii_6548_flux_ext_err', 
              'h_alpha_flux_ext_err', 'nii_6584_flux_ext_err', 
              'sii_6717_flux_ext_err']

#Defining IDs for different Flux Lines
linelist_full = ['oiii4363' , 'oiii4959', 'hbeta','oiii5007',  'nii6584', 'halpha', 'nii6548', 'sii6731;sii6731' ]

#linelist_full = ['oii3726', 'oii3729', 'neiii3869', 'oiii4363', 'hgamma', 
#        'hbeta', 'oiii4959', 'oiii5007', 'hei5875', 'oi6300', 'nii6548', 
#        'halpha', 'nii6584', 'sii6717;sii6731', 'ariii7136']

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
if sys.platform == 'linux2':
    f1 = open('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results/IZI_Z_3.txt', 'w')
    f2 = open("/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results/IZI_q_3.txt", "w")  
    f3 = open("/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/results/IZI_results_3.pkl", "w")  
else:
    f1 = open('C:/Users/mugdhapolimera/github/izi/results/IZI_Z_prior3.txt', 'a+')
    f2 = open('C:/Users/mugdhapolimera/github/izi/results/IZI_q_prior3.txt', 'a+')
    f3 = open('C:/Users/mugdhapolimera/github/izi/results/IZI_results_prior3.pkl', 'a+')

#f1.write("#Name \t Z_Estimate \t Err_down \t Err_up\r\n")
#f2.write("#Name \t q_Estimate \t Err_down \t Err_up\r\n")
#infile['NAME'] = infile['col0']
print len(infile['NAME'])
t1 = time.time()         

def mz(m):
    if (m < 8.8):
        p = np.poly1d([0.53518733, 3.66817274])
        z = p(m)
    else:
        m = m-10    
        z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
    #p = np.poly1d(np.polyfit(x, z[np.where((m<9.6) & (m>9.0))[0]], 1))    
    return z

Z = [mz(m) for m in infile.logmstar]
Z_hist = np.histogram(Z, bins = 'fd')
logzprior = np.column_stack((Z_hist[1][:-1],Z_hist[0]))
plt.hist(Z,bins = 'fd')
plt.plot(logzprior[:,0], logzprior[:,1])
zprior = np.poly1d(np.polyfit(logzprior[:,0], logzprior[:,1],4))
zlim = np.linspace(min(logzprior[:,0])-0.1, max(logzprior[:,0])+0.1, 100)
plt.plot(zlim,zprior(zlim))
for gal in range(len(infile['NAME'])):
    tgal1 = time.time()             
    fluxin = []
    errorin = []
    idin = []
    print gal, infile['NAME'][gal]    
    for i in range(len(fluxnames)):
        if (infile[fluxnames[i]][gal] > 0) & (infile[errornames[i]][gal] > 0) & (infile[fluxnames[i]][gal]/infile[errornames[i]][gal] > 2)  : 
            if fluxnames[i] == 'sii_6717_flux_ext':
                fluxin.append((infile[fluxnames[i]][gal]+infile['sii_6731_flux_ext'][gal])/2)
                errorin.append((infile[errornames[i]][gal]**2 + infile['sii_6731_flux_ext_err'][gal]**2)**0.5)
            else:   
                fluxin.append(infile[fluxnames[i]][gal])
                errorin.append(infile[errornames[i]][gal])
            idin.append(linelist_full[i])
    #print fluxin
    #print errorin
    #kwargs =  logOHsun, intergridfile, logzlimits, logqlimits ,logzprior, logqprior (logz/q prior have no application in the original code)
    #d = izi(fluxin, errorin, idin, plot_flag = 0, name = str(infile['NAME'][gal]), intergridfile = '/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/outputgrid_linear.fits', interpolate_flag = False)
    
    #zlim = np.linspace(min(infile.logmstar)-8.66, max(infile.logmstar)-8.66, 100)
    #zprior = norm.pdf(zlim, loc = mz(infile.logmstar[gal])-8.66, scale = 10**0.08)
    #logzprior = np.column_stack((zlim,zprior))
    d = izi(fluxin, errorin, idin, name = str(infile['NAME'][gal]), 
            plot_flag = 0, print_flag = 0, method = 'gpy', logzprior = logzprior)

    #print d
#    f1.write(" %s \t %f \t %f \t %f \r\n" %(d['name'], d.Zgrid, d.edownZgrid, d.eupZgrid))
#    f2.write(" %s \t %f \t %f \t %f \r\n" %(d['name'], d.qgrid, d.edownqgrid, d.eupqgrid))
#    pickle.dump(d, f3)    
    tgal2= time.time()         
    print 'Time for ', infile['NAME'][gal], ': ',tgal2-tgal1
f1.close()
f2.close()
f3.close()
t2 = time.time()
print t2-t1