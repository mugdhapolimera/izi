

input = getenv('izi') + 'resolvecatalog_str.fits'

file=mrdfits(input, 1)
idin = ['oiii5007', 'oiii4959', 'oiii4363' , 'nii6584',  'nii6548', 'hbeta', 'halpha' ]


for i=0, 0 do begin
  fluxin = [file[i].F_OIII_5007_BROAD, file[i].F_OIII_4960_BROAD, file[i].F_OIII_4363_BROAD,  file[i].F_NII_6586_BROAD,  file[i].F_NII_6548_BROAD,  file[i].F_HB_BROAD, file[i].F_HA_BROAD]
  errorin = [file[i].F_OIII_5007_ERR_BROAD, file[i].F_OIII_4960_ERR_BROAD, file[i].F_OIII_4363_ERR_BROAD, file[i].F_NII_6586_ERR_BROAD, file[i].F_NII_6548_ERR_BROAD, file[i].F_HB_ERR_BROAD, file[i].F_HA_ERR_BROAD]
  outgridfile = 'outputgrid_linear.fits'
  result = izi(fluxin, errorin, idin, /plot, /outgridfile, /linear)
  print, result
endfor 

end