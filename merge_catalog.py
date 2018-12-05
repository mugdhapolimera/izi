import numpy as np
from astropy.io import fits
import pandas as pd
import os
#data1 = np.genfromtxt("RESOLVE_liveOctober2018.csv", delimiter=",",dtype=None,names=True)
os.chdir('C:\Users\mugdhapolimera\github\SDSS_Spectra')
df1 = pd.read_csv("RESOLVE_liveOctober2018.csv",index_col='name')
print df1

#data2 = np.genfromtxt("RESOLVE_izioutv1.txt", delimiter="",dtype=None,names=True)
#df2 = pd.DataFrame(data=data2,index=data2['name'])

#dftemp = pd.merge(df1,df2,how="left",left_index=True,right_index=True)

#hdulist = fits.open("ECO_SDSS_full_dext.fits")
#name3 = hdulist[1].data['NAME']
#sort3 = np.argsort(name3)
#data3 = hdulist[1].data[sort3]
##newdata3 = data3.byteswap().newbyteorder()
#df3 = pd.DataFrame(data=newdata3,index=newdata3['NAME'])
df3 = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_Spectra\RESOLVE_SDSS_full_raw_flux.csv')
df3.index = np.array(df3['name'])
dfres = pd.merge(df1,df3,how="left",left_index=True,right_index=True)
#print dfres.h_alpha_flux
#dfres.to_pickle('RESOLVE_full_raw.pkl')


