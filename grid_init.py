from astropy.table import Table
import numpy as np
gridfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\l09_high_csf_n1e2_6.0Myr.fits'
grid0 = Table.read(gridfile, format='fits')

new_gridfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\l09_high_csf_n1e2_6.0Myr_new.fits'
cols = []
cols.append(list(grid0['LOGZ']+grid0['LOGOHSUN']))
cols.append(list(grid0['LOGQ']))
for i in range(len(grid0['FLUX'][0])):
    cols.append(list(grid0['FLUX'][:,i]))#.transpose()
col_names = ['LOGZ', 'LOGQ'] + list( grid0['ID'][0])

t = Table(cols, names=(col_names))
print t.keys()
#t.write(new_gridfile, format='fits')