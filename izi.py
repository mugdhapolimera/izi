#IZI Gas-Phase Metallicity Estimator -version Python
#---------------------------------------------------

import numpy as np
import math
import matplotlib.pyplot as plt
#import PyQt4
import os
from astropy.table import Table
import pandas as pd
import scipy
import scipy.interpolate

# Setting System Path
import sys
sys.path.append('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/izi_utils/')
#sys.path.append('C:\Users\mugdhapolimera\github\izi\izi_utils')

#Importing Custom Utility Files
#from izi_utils.tabulate import idl_tabulate
from tabulate import idl_tabulate
from interpolate import grid_interpolate
import izi_plots
#---------------------------------------------------


def uprior(xaxis):
    return 1./(xaxis[1]-xaxis[0])


#def izi(fluxin, errorin, idin, gridfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\l09_high_csf_n1e2_6.0Myr.fits', 
def izi(fluxin, errorin, idin, name, gridfile = '/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/grids/l09_high_csf_n1e2_6.0Myr.fits', 
        plot_flag = 1, print_flag = True, epsilon = 0.15, nz1 = 50, nq1 = 50, interpolate_flag = True, outgridfile = True, nonorm = False, **kwargs ):
    
#kwargs =  logOHsun, intergridfile, logzlimits, logqlimits ,logzprior, logqprior (logz/q prior have no application in the original code)

    # RENAME INPUT ARRAYS
    flux = fluxin
    error = errorin
    idno = idin 

    #TODO Make IZI a function
    #CHECK INPUT FOR CONSISTENCY
    nlines=len(flux)
    if (len(error) != nlines | len(idno) != nlines): 
        print 'ERROR: Flux, Error, and ID arrays do not have the same number of elements'

    #DEFAULT GRID: 
    #Levesque 2010, HIGH MASS LOSS, CSF 6Myr, n=100 cm^-3

    #READ GRID
    if 'intergridfile' in kwargs.keys():
        gridfile = kwargs['intergridfile']
    else:
        gridfile = gridfile

    grid0 = Table.read(gridfile, format='fits')
    
    id0=grid0['ID'][0]
    nlines0 = len(id0)
    for i in range (nlines0-1):
        'in'+id0[i]+"=where(id0 eq '"+id0[i]+"')"

    grid0['ID'] = [np.char.strip(x) for x in grid0['ID']]
    ngrid=len(grid0['LOGZ'])


    if 'logOHsun' in kwargs.keys():
        logOHsun = kwargs['logOHsun']
    else:
        logOHsun = grid0['LOGOHSUN'][0] 
    
    
    #CUT GRID TO LOGZLIMITS AND LOGQLIMITS 
    try:
        logzlimits
    except NameError: 
        logzlimits = [min(grid0['LOGZ']+logOHsun), max(grid0['LOGZ']+logOHsun)]

    try:
        logqlimits
    except NameError:    
        logqlimits = [min(grid0['LOGQ']), max(grid0['LOGQ'])]

    grid0=grid0[np.where((grid0['LOGZ']+logOHsun >= logzlimits[0]) & (grid0['LOGZ']+logOHsun <= logzlimits[1]) & 
                         (grid0['LOGQ'] >= logqlimits[0]) & (grid0['LOGQ'] <= logqlimits[1]))]

    #CHANGE LOGZPRIOR TO SOLAR UNITS
    #TODO: Why is this a comment?
    #logZprior[:,0]=logZprior-logOHsun

    #INCLUDE SYSTEMATIC UNCERTAINTY IN THE PHOTO-IONIZATION MODELS
    # default is 0.15 dex systematic uncertainty
    epsilon2 = epsilon*math.log(10) # convert to scaling factor

    #Default number of interpolation steps
    nz1 = 50
    nq1 = 50
    zarr=np.linspace(min(grid0['LOGZ']), max(grid0['LOGZ']), nz1)
    qarr=np.linspace(min(grid0['LOGQ']), max(grid0['LOGQ']), nq1)
    
    if (interpolate_flag):    
        grid, ngrid, zarr, qarr, dlogz, dlogq = grid_interpolate(grid0, method = 'scipy', nz1 = nz1, nq1 = nq1)
    else:
        grid  = grid0
        ngrid = len(grid['LOGZ'])
        zarr  = np.unique(grid0['LOGZ'])
        qarr  = np.unique(grid0['LOGQ'])
        nz1   = len(zarr)
        nq1   = len(qarr)
        dlogz = (zarr[nz1-1]-zarr[0])/(nz1-1)
        dlogq = (qarr[nq1-1]-qarr[0])/(nq1-1)
        
    #CREATE DATA STRUCTURE CONTAINING LINE FLUXES AND ESTIMATED PARAMETERS
    d = pd.Series({      'name'             : str(name),
                         'id'               : id0,
                         'flux'             : np.zeros(nlines0) -666,
                         'error'            : np.zeros(nlines0) -666,
                         'chi2'             : 0., 
                         'Zgrid'            : 0., 
                         'eupZgrid'         : 0.,
                         'edownZgrid'       : 0.,
                         'qgrid'            : 0., 
                         'eupqgrid'         : 0.,
                         'edownqgrid'       : 0.,
                         'Zgridmarmod'      : 0.,
                         'eupZgridmarmod'   : 0.,
                         'edownZgridmarmod' : 0.,
                         'qgridmarmod'      : 0.,
                         'eupqgridmarmod'   : 0.,
                         'edownqgridmarmod' : 0.,
                         'Zgridmarmean'     : 0.,
                         'eupZgridmarmean'  : 0.,
                         'edownZgridmarmean': 0.,
                         'qgridmarmean'     : 0.,
                         'eupqgridmarmean'  : 0.,
                         'edownqgridmarmean': 0.,
                         'zarr'             : np.zeros(nz1),
                         'zpdfmar'          : np.zeros(nz1),
                         'qarr'             : np.zeros(nq1),
                         'qpdfmar'          : np.zeros(nq1),
                         'flag'            : [0,0,1,1], 
                         'pdfjoint'         : np.zeros([nz1,nq1]) 
                    })

    grid['ID'] = [np.char.strip(x) for x in grid['ID']]
    #FILL STRUCTURE WITH LINE FLUXES
    for i in range(nlines):
        #auxind=np.where(d.id == idno[i])[0]
        auxind = [x for x,item in enumerate(d.id) if idno[i] in item]
        if (auxind == 0):
            print 'ERROR: ===== Line ID '+idno[i]+'not recognized ====='
        d.flux[auxind]=flux[i]
        d.error[auxind]=error[i]

    # INDEX OF LINES WITH MEASUREMENTS
    good = np.where(d.error != -666)[0]
    ngood = len(good)
    measured = np.where(d.flux != -666)[0]
    nmeasured = len(measured)
    upperlim = np.where((d.error != -666) & (d.flux == -666))[0]
    flag0=np.zeros(nlines0, dtype = float)
    if (measured != []):
        flag0[measured] = 1      #measured flux
    if (upperlim == []):
        flag0[upperlim] = 2      #upper limit on flux
    flag=flag0[good]

    # NORMALIZE LINE FLUXES TO H-BETA OR
    # IF ABSENT NORMALIZE TO BRIGHTEST LINE
    if not (nonorm):
        if print_flag:
            print 'Normalizing Fluxes'

        idnorm = 'hbeta'
        if (d.flux[[x for x,item in enumerate(d.id) if idnorm in item][0]] == -666):
            idnorm = (d.id[measured])[np.argsort(d.flux[measured])][::-1][0] 

        #normalize data
        norm = d.flux[[x for x,item in enumerate(d.id) if idnorm in item][0]]
        d.flux[measured] = d.flux[measured]/norm
        d.error[good] = d.error[good]/norm

        #normalize grid
        for i in range(ngrid):
            norm = grid[i]['FLUX'][[x for x,item in enumerate(grid[i]['ID']) if idnorm in item][0]]
            grid[i]['FLUX'] = grid[i]['FLUX']/norm

    like=np.ones(ngrid, dtype = np.double)
    post=np.ones(ngrid, dtype = np.double)
    zrange=[min(grid['LOGZ']), max(grid['LOGZ'])]
    qrange=[min(grid['LOGQ']), max(grid['LOGQ'])]
    
    for i in range(ngrid):
        for j in range(ngood):
                            #CALCULATE LIKELIHOOD
            if (flag[j] == 1): # If measured          
                normalization = np.sqrt(d.error[good][j]**2+(epsilon2*grid[i]['FLUX'][good][j])**2)
                flux_diff = (d.flux[good][j] - grid[i]['FLUX'][good][j])**2
                error_quad = (d.error[good][j]**2 + (epsilon2*grid[i]['FLUX'][good][j])**2)
                exponent = np.exp(-1.0*flux_diff/(2.0*error_quad))
                like[i] = like[i]*1.0/np.sqrt(2.0*3.14)*exponent/normalization
                
            if (flag[j] == 2): # if upper limit
                like[i] = like[i]*0.5*( 1 + scipy.special.erf((d.error[good][j] - grid[i]['FLUX'][good][j])/(np.sqrt(d.error[good][j]**2+(epsilon2*grid[i]['FLUX'][good][j])**2)*np.sqrt(2))))
                            #CALCULATE POSTERIOR BY INCLUDING PRIORS AND NORMALIZING
        if (('logzprior' in locals()) == 0) & (('logqprior' in locals()) == 0):
            post[i] = uprior(zrange)*uprior(qrange)*like[i]
        if (('logzprior' in locals()) == 1) & (('logqprior' in locals()) == 0):
            post[i] = userprior(grid['LOGZ'][i], logzprior[:,0], logzprior[:,1])*uprior(qrange)*like[i]
        if (('logzprior' in locals()) == 0) & (('logqprior' in locals()) == 1):
            post[i] = uprior(zrange)*userprior(grid[i].logq, logqprior[:,0], logqprior[:,1])*like[i]
        if (('logzprior' in locals()) == 1) & (('logqprior' in locals()) == 1):
            post[i] = userprior(grid[i].logz, logzprior[:,0], logzprior[:,1])*userprior(grid['LOGQ'][i], logqprior[:,0], logqprior[:,1])*like[i]
    like[np.where(np.isfinite(like) == 0)]=0
    post[np.where(np.isfinite(post) == 0)]=0

    goodlike = np.where(np.isfinite(like))[0]
    sortlike = like[goodlike][np.argsort(like[goodlike])[::-1]]
    sortz = grid['LOGZ'][goodlike][np.argsort(like[goodlike])[::-1]]
    sortq = grid['LOGQ'][goodlike][np.argsort(like[goodlike])[::-1]]
    sumlike=np.zeros(len(sortlike))
    for i in range (len(sortlike)):
        sumlike[i]=np.sum(sortlike[:i])/np.sum(sortlike)                                 

    goodpost = np.where(np.isfinite(post))[0]  
    sortpost = (post[goodpost])[np.argsort(post[goodpost])[::-1]]
    sortz = np.array(grid['LOGZ'][goodpost][np.argsort(post[goodpost])[::-1]])
    sortq = np.array(grid['LOGQ'][goodpost][np.argsort(post[goodpost])[::-1]])
    sumpost = np.zeros(len(sortpost))
    for i in range(len(sortpost)):
        sumpost[i] = np.sum(sortpost[0:i])/np.sum(sortpost) 

    # CALCULATE BEST FIT METALLICITY, IONIZATION PARAMETER AND ERRORS

    post1sig=(sortpost[np.where(sumpost >= 0.683)])#[0]
    post2sig=(sortpost[np.where(sumpost >= 0.955)])#[0]
    post3sig=(sortpost[np.where(sumpost >= 0.997)])#[0]

    like1sig=(sortlike[np.where(sumlike >= 0.683)])#[0]
    like2sig=(sortlike[np.where(sumlike >= 0.955)])#[0]
    like3sig=(sortlike[np.where(sumlike >= 0.997)])#[0]

    d.Zgrid = sortz[0]+logOHsun
    d.edownZgrid = sortz[0]-min( list(sortz[np.where(sumpost <= 0.683)]) or [0])
    d.eupZgrid = max(list(sortz[np.where(sumpost <= 0.683)])  or [0])-sortz[0]

    d.qgrid = sortq[0]
    d.edownqgrid = sortq[0]-min(list(sortq[np.where(sumpost <= 0.683)])  or [0])
    d.eupqgrid = max(list(sortq[np.where(sumpost <= 0.683)])  or [0])-sortq[0]


    # COMPUTE chi2

    bestgrid = np.where((grid['LOGZ'] == sortz[0]) & (grid['LOGQ'] == sortq[0]))[0][0]
    fobs = d.flux[np.where(d.flux != -666)[0]]  
    eobs = d.error[np.where(d.flux != -666)[0]]  
    fmod = grid['FLUX'][bestgrid][np.where(d.flux != -666)[0]]
    emod = epsilon2*fmod
    d.chi2 = np.sum((fobs-fmod)**2/(eobs**2+emod**2))/len(fobs)

    # posterior for Z, marginalizing over q
    postz = np.zeros(nz1, dtype = np.float64)
    for j in range(nz1):
        postz[j] = idl_tabulate(sortq[np.where(sortz==zarr[j])[0]],sortpost[np.where(sortz==zarr[j])[0]])
    postz = postz/np.sum(postz)

    sumpostz = np.zeros(len(postz))
    sumpz = 0
    for i in range(nz1):
        sumpz += postz[i]
        sumpostz[i] += sumpz 

    d.Zgridmarmod = zarr[np.where(postz == max(postz))[0]] + logOHsun # max of PDF
    d.Zgridmarmean = np.sum(zarr*postz)/np.sum(postz) + logOHsun # first moment of PDF

    if len(np.where(sumpostz >= (1.0-0.683)/2.0)[0]) != 0:
        d.edownZgridmarmod = d.Zgrid - logOHsun - zarr[np.where(sumpostz >= (1.0-0.683)/2.0)[0][0]]
        d.eupZgridmarmod = zarr[np.where(sumpostz >= (1.0-(1.0-0.683)/2.0))[0][0]] - d.Zgrid + logOHsun
        d.edownZgridmarmean = d.Zgrid - logOHsun - zarr[np.where(sumpostz >= (1.0-0.683)/2.0)[0][0]]
        d.eupZgridmarmean = zarr[np.where(sumpostz >= (1.0-(1.0-0.683)/2.0))[0][0]] - d.Zgrid + logOHsun
    else:
        d.edownZgridmarmod = d.eupZgridmarmod = d.edownZgridmarmean = d.eupZgridmarmean = 0
    #posterior for q, marginalizing over Z
    postq = np.zeros(nq1)
    for j in range(nq1):
        postq[j] = idl_tabulate(sortz[np.where(sortq == qarr[j])[0]], sortpost[np.where(sortq == qarr[j])[0]])
    postq = postq/np.sum(postq)

    sumpostq = np.zeros(len(postq))
    sumpq = 0
    for i in range(nq1):
        sumpq += postq[i]
        sumpostq[i] = sumpq 

    d.qgridmarmod = qarr[np.where(postq == max(postq))[0]] #MAx of PDF
    d.qgridmarmean = np.sum(qarr*postq)/np.sum(postq) # first moment of PDF

    if len(np.where(sumpostz >= (1.0-0.683)/2.0)[0]) != 0:
        d.edownqgridmarmod = d.qgrid - qarr[np.where(sumpostq >= ((1.0-0.683)/2.0))[0][0]]
        d.eupqgridmarmod = qarr[np.where(sumpostq >= (1.0-(1.0-0.683)/2.0))[0][0]] - d.qgrid
        d.edownqgridmarmean = d.qgrid - qarr[np.where(sumpostq >= ((1.0-0.683)/2.0))[0][0]]
        d.eupqgridmarmean = qarr[np.where(sumpostq >= (1.0-(1.0-0.683)/2.0))[0][0]] - d.qgrid
    else:
        d.edownqgridmarmod = d.eupqgridmarmod = d.edownqgridmarmean = d.eupqgridmarmean = 0

    # WRITE MARGINALIZED PDFS
    d.zarr=zarr+logOHsun
    d.zpdfmar=postz
    d.qarr=qarr
    d.qpdfmar=postq

    # Set FLAGS to warn bout multiple peaks and lower/upper limits
    dzpdf=np.gradient(postz)
    ddzpdf=np.gradient(dzpdf)
    auxpeak = np.zeros(nz1, dtype=np.int)
    for i in range (nz1-1):
        if ((dzpdf[i] > 0) & (dzpdf[i+1] < 0) & (ddzpdf[i] < 0)):
            auxpeak[i] = 1
    zpeaks=np.where(auxpeak == 1)[0]
    d.flag[0]=len(zpeaks)

    dqpdf=np.gradient(postq)
    ddqpdf=np.gradient(dqpdf)
    auxpeak = np.zeros(nq1, dtype=np.int)
    for i in range (nq1-1):
        if ((dqpdf[i] > 0) & (dqpdf[i+1] < 0) & (ddqpdf[i] < 0)):
            auxpeak[i] = 1
    qpeaks=np.where(auxpeak == 1)[0]
    d.flag[1]=len(qpeaks)

    if (max(postz[0:1]) > 0.5*max(postz)): 
        d.flag[2]=2
    if (max(postz[nz1-2:nz1-1]) > 0.5*max(postz)):
        d.flag[2]=3
    if ((max(postz[0:1]) > 0.5*max(postz)) & (max(postz[nz1-2:nz1-1]) > 0.5*max(postz))):
        d.flag[2]=0
    if (max(postq[0:1]) > 0.5*max(postq)): 
        d.flag[3]=2
    if (max(postq[nq1-2:nq1-1]) > 0.5*max(postq)):
        d.flag[3]=3
    if ((max(postq[0:1]) > 0.5*max(postq)) & (max(postq[nq1-2:nq1-1]) > 0.5*max(postq))): 
        d.flag[3]=0
    
    if print_flag:
        print '===== BEST FIT FROM JOINT PDF MODE ====='
        print '===== Z =====', d.Zgrid, d.edownZgrid, d.eupZgrid, d.Zgrid-d.edownZgrid, d.Zgrid+d.eupZgrid 
        print '===== q =====', d.qgrid, d.edownqgrid, d.eupqgrid, d.qgrid-d.edownqgrid, d.qgrid+d.eupqgrid
        c=2.99792458e10 # cm/s 
        print '=== U=q/c ===', d.qgrid-np.log10(c), d.edownqgrid, d.eupqgrid, d.qgrid-np.log10(c)-d.edownqgrid, d.qgrid-np.log10(c)+d.eupqgrid
        print '================ FLAGS ====================='
        print d.flag
        print '============================================'

    # PLOT RESULTS
    directory = '/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/izi_plots/'+str(d['name'])
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)
    if (np.sum(post) == 0. or np.sum(like) == 0.):
        return d
    else:
        izi_plots.izi_pdf(d = d, grid = grid, postz = postz, postq = postq, post = post, like = like, plot_flag = plot_flag)
        izi_plots.zratios_plots(grid = grid, grid0 = grid0, d = d, flag0 = flag0, plot_flag = plot_flag)
        izi_plots.qratios_plots(grid = grid, grid0 = grid0, d = d, flag0 = flag0, plot_flag = plot_flag)
    
    return d