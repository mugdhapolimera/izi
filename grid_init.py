from astropy.table import Table
import pandas as pd
import numpy as np
gridfile = 'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\grids\l09_high_csf_n1e2_6.0Myr.fits'
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

'''
gridfile = 'C:/Users/mugdhapolimera/github/izi/Richardson_bpass_binary_csf_n1e2_40.0Myr.csv'
df = pd.read_csv(gridfile)
lines = "cii1335  civ1548  civ1551  ciii1909 ciii1911 cii2324  oii3726  oii3729  neiii3869 sii4069  hgamma   oiii4363 hei4471  hbeta    oiii4959 oiii5007 hei5016  ariii5192 ni5198   nii5755  hei5875  oi6300   siii6312 nii6548  halpha   nii6584  hei6678  sii6717  sii6731  ariii7136 oii7318  oii7320  ariii7751 siii9068 siii9532"
line_id = lines.split()
for i in range(len(df)):
    df['FLUX'][i] = float(df['FLUX'][i][3:])
    df['Unnamed: 39'][i] = float(df['Unnamed: 39'][i][:-2])
  
for i in range(5,len(df.columns)):
    df.rename(columns={df.columns[i]:line_id[i-5]}, inplace=True)    

df['LOGZ'] = df['LOGZ'] +df['LOGOHSUN']    
matched_Z = [df['LOGZ'][x] - 8.76 in np.around(np.unique(grid0['LOGZ']),5) for x in range(len(df['LOGZ']))]
matched_q = [(df['LOGQ'][x] >= 6.9) & (df['LOGQ'][x] <= 8.8)for x in range(len(df['LOGQ']))]
df = df[np.array(matched_Z)]# & np.array(matched_Z)]
df = df.drop(columns = ['NAME','LOGOHSUN','ID'])
df = df.apply(pd.to_numeric)
#for i in range(len(df)):
    #if (df['LOGZ'][i] < min(grid0['LOGZ'])+grid0['LOGOHSUN'][0]) | (df['LOGZ'][i] > max(grid0['LOGZ'])+grid0['LOGOHSUN'][0]) : 
#        df = df.drop([i])
df.index = range(len(df))
new_gridfile = r'C:\Users\mugdhapolimera\github\izi\Richardson_bpass_binary_csf_n1e2_40.0Myr_limitedZ.fits'
t = Table.from_pandas(df)
print t
print t.keys()
t.write(new_gridfile, format='fits')
'''