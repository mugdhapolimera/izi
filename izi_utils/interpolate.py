import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate
import GPy as gpy
from astropy.table import Table

def grid_interpolate(grid0, plot = False, method = 'scipy', nz1 = 50, nq1 = 50):
    
    #CREATE AN EMPTY GRID 
    nlines0=len(grid0['ID'][0])
    
    id_str = []
    for i in range(nz1*nq1):
        id_str.append(grid0['ID'][0])
    id_str = np.array(id_str)

    names = []
    for i in range(nz1*nq1):
        names.append(grid0[0]['NAME'])    
    names = np.array(names)

    fluxlines = []
    for i in range(nz1*nq1):
        fluxlines.append(np.zeros(nlines0))    
    fluxlines = np.array(fluxlines)

    data = { 'NAME'    : names,
             'LOGOHSUN': np.ones(nz1*nq1)*grid0['LOGOHSUN'][0],
             'LOGZ'    : np.zeros(nz1*nq1),
             'LOGQ'    : np.zeros(nz1*nq1),
             'ID'      : id_str,
             'FLUX'    : fluxlines}

    grid = Table(data, names = ('NAME','LOGOHSUN', 'LOGZ', 'LOGQ', 'ID', 'FLUX'))

    
    #GET LINE IDs IN THE GRID AND INDEX OF EACH LINE
    if (method == 'gpy'):
        Z = np.array(grid0['LOGZ'])
        q = np.array(grid0['LOGQ'])
        zarr=np.linspace(min(grid0['LOGZ']), max(grid0['LOGZ']), nz1)
        qarr=np.linspace(min(grid0['LOGQ']), max(grid0['LOGQ']), nq1)
        

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
            
            dlogz = (zarr[nz1-1]-zarr[0])/(nz1-1)
            dlogq = (qarr[nq1-1]-qarr[0])/(nq1-1)
            
            ngrid += 1
            test_output = m.predict(test_input)[0]
            if (np.isnan(test_output).all()):
                test_output = np.ones(len(test_input))
                
            grid['FLUX'][:,no] = test_output.reshape(ngrid)

        
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
        grid['LOGZ'] = test_input[:,0]
        grid['LOGQ'] = test_input[:,1]


    elif (method == 'scipy'):
        
        Z = np.array(grid0['LOGZ'])
        q = np.array(grid0['LOGQ'])
        x = np.unique(Z)
        y = np.unique(q)

        nlines0 = len(grid0['ID'][0])
        
        z = np.zeros([nlines0, len(y), len(x)])
        fluxarr = np.zeros([nlines0,nz1,nq1])
        zarr=np.linspace(min(grid0['LOGZ']), max(grid0['LOGZ']), nz1)
        qarr=np.linspace(min(grid0['LOGQ']), max(grid0['LOGQ']), nq1)
        gridx, gridy = np.meshgrid(zarr, qarr)

        X, Y = np.meshgrid(x,y)

        for k in range (nlines0):
            for i in range(len(x)): 
                for j in range(len(y)):
                    ind = np.where((grid0['LOGZ'] == x[i]) & (grid0['LOGQ'] == y[j]))[0]

                    z[k][j][i] = grid0['FLUX'][:,k][ind]

            f = scipy.interpolate.interp2d(x,y,z[k,:,:],'linear')
            fluxarr[k,:,:] = f(zarr,qarr)
            fluxarr[k] = fluxarr[k].transpose()

            if plot:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_surface(gridx, gridy, fluxarr[k])
                ax.scatter(X,Y,z[:,:,k],'r')
        
        for i in range(nz1):
            for j in range(nq1):
                ngrid = i*nq1 + j
                grid[ngrid]['LOGZ'] = zarr[i]
                grid[ngrid]['LOGQ'] = qarr[j]
                for k in range(nlines0):
                    grid[ngrid]['FLUX'][k] = fluxarr[k,i,j]
                    #print 'yes'
        dlogz = (zarr[nz1-1]-zarr[0])/(nz1-1)
        dlogq = (qarr[nq1-1]-qarr[0])/(nq1-1)
        ngrid=ngrid+1

    
    return grid, ngrid, zarr, qarr, dlogz, dlogq    
