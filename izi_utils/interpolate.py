import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
import GPy as gpy
from astropy.table import Table

def grid_interpolate(grid0, plot = False, method = 'scipy', nz1 = 50, nq1 = 50):
    
    #GET LINE IDs IN THE GRID AND INDEX OF EACH LINE
    id0=grid0['ID'][0]
    nlines0=len(id0)
    for i in range (nlines0-1):
        'in'+id0[i]+"=where(id0 eq '"+id0[i]+"')"
    
    if (method == 'gpy'):
        
        Z = np.array(grid0['LOGZ'])
        q = np.array(grid0['LOGQ'])

        # Fit a GP
        
        #Defining a Matern 3/2 + Constant Kernel
        km = gpy.kern.Matern32(input_dim=2, ARD = True)
        kb = gpy.kern.Bias(input_dim=2)
        k = km + kb

        train_input = np.zeros([len(Z), 2])
        train_input[:,0] = Z
        train_input[:,1] = q
        
        # Compute the model prediction on training input
        for no in range(nlines0):
            train_output = np.array(grid0['FLUX'][:,no])[:,None]
            m = gpy.models.GPRegression(train_input, train_output, k, normalizer=True)
            #print m
            #m.rbf.variance.constrain_bounded(1e-3, 1e5)
            #m.bias.variance.constrain_bounded(1e-3,1e5)
            #m.rbf.lengthscale.constrain_bounded(.1,200.)
            m.Gaussian_noise.variance.constrain_fixed(1e-5)
            m.randomize()
            
            if not(np.isnan(m.Y_normalized).all()):
                m.optimize()
            
            test_input = np.zeros([nz1*nq1,2])
            for i in range(nz1):
                for j in range(nq1):
                    ngrid = i*nq1 + j
                    test_input[ngrid,0] = zarr[i]
                    test_input[ngrid,1] = qarr[j]
            test_output = m.predict(test_input)[0]
            if (np.isnan(test_output).all()):
                test_output = np.ones(len(test_input))
        
            if (plot):
                #Plot GP Model and 2-D Interpolation Results
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                #ax.plot_surface(gridx1, gridy1, flux_idl[k], color = 'b')
                ax.plot_surface(gridx1, gridy1, flux_idl[no].transpose(), color = 'b')
                ax.scatter(test_input[:,0],test_input[:,1],test_output,c = 'r')
                ax.plot_surface(gridx, gridy, fluxarr[no], color = 'g')
                ax.set_xlabel('Log(Z)')
                ax.set_ylabel('Log(q)')
                ax.set_zlabel('Relative'+grid0[0]['ID'][no]+'Intensity')

                #Calculate and Plot Residuals
                residual = np.zeros(np.shape(test_output))
                for i in range(len(zarr)):
                    for j in range(len(qarr)):
                        ind = np.where((test_input[:,0] == zarr[i]) & (test_input[:,1] == qarr[j]))
                        residual[ind] = test_output[ind] - flux_idl[no][i][j]
                
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(test_input[:,0],test_input[:,1],residual,c = 'r')
                ax.set_xlabel('Log(Z)')
                ax.set_ylabel('Log(q)')
                ax.set_zlabel('Residuals of Relative'+grid0[0]['ID'][no]+'Intensity')

    elif (method == 'scipy'):
  
        Z = np.array(grid0['LOGZ'])
        q = np.array(grid0['LOGQ'])

        x = np.unique(Z)
        y = np.unique(q)
        X, Y = np.meshgrid(x,y)
        nlines0 = 35
        z = np.zeros([len(y), len(x),nlines0])

        fluxarr = np.zeros([nlines0,nz1,nq1])
        zarr=np.linspace(min(grid0['LOGZ']), max(grid0['LOGZ']), nz1)
        qarr=np.linspace(min(grid0['LOGQ']), max(grid0['LOGQ']), nq1)
        gridx, gridy = np.meshgrid(zarr, qarr)

        for k in range (nlines0):
            for i in range(len(x)):
                for j in range(len(y)):
                    n = np.where((q == y[j]) & (Z == x[i]))[0]
                    z[j,i,k] = grid0['FLUX'][n,k]

            f = scipy.interpolate.interp2d(x,y,z[:,:,k],'linear')
            fluxarr[k,:,:] = f(zarr,qarr)

            if (plot):
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_surface(gridx, gridy, fluxarr[k])
                ax.scatter(X,Y,z[:,:,k],'r')

    #CREATE AN EMPTY GRID 
    id_str = []
    for i in range(nz1*nq1):
        id_str.append(id0)
    id_str = np.array(id_str)

    names = []
    for i in range(nz1*nq1):
        names.append(grid0[0]['NAME'])    
    names = np.array(names)

    grid_flux = []
    for i in range(nz1*nq1):
        grid_flux.append(np.zeros(nlines0))    
    grid_flux = np.array(grid_flux)
    
    logOHsun = grid0['LOGOHSUN'][0]
    
    data = { 'NAME'     : names,
             'LOGOHSUN' : logOHsun*np.ones(nz1*nq1),
             'LOGZ'     : np.zeros(nz1*nq1),
             'LOGQ'     : np.zeros(nz1*nq1),
             'ID'       : id_str,
             'FLUX'     : grid_flux}

    grid = Table(data, names = ('NAME', 'LOGOHSUN', 'LOGZ', 'LOGQ', 'ID', 'FLUX'))

    for i in range(nz1):
        for j in range(nq1):
            ngrid = i*nq1 + j
            grid[ngrid]['LOGZ'] = zarr[i]
            grid[ngrid]['LOGQ'] = qarr[j]
            for k in range(nlines0):
                grid[ngrid]['FLUX'][k] = fluxarr[k,i,j]
            #print 'yes'

    ngrid=ngrid+1
    return grid, ngrid    
