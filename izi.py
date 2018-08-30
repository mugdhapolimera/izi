#IZI Gas-Phase Metallicity Estimator -version Python
#---------------------------------------------------

#Import Packages

import numpy as np
import math
import matplotlib as plt
#---------------------------------------------------

def uprior(xrange):
    return 1d/(xrange[1]-xrange[0]);

def jprior(x, xrange):
    return 1d/(x*math.log10(xrange[1]/xrange[0]));

def userprior(x,xarr,yarr):
    if (x < min(xarr)) | (x > max(xarr)):
        return 0;
    return interpol(yarr,xarr,x,nan); 

def  plotratioq (ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, range=range, title=title)
# plot a line ratio vs the ionization parameter

  if !(keyword_set(range)):
      range=[-2.0, 2.0]
  if !(keyword_set(title)):
      title='log(ratio)'
  
  auxzarr=grid[uniq(grid.logz, sort(grid.logz))].logz
  auxzarr0=grid0[uniq(grid0.logz, sort(grid0.logz))].logz

  sel1=np.where(abs(grid.logz+logOHsun-d.zgrid) ==  min(abs(grid.logz+logOHsun-d.zgrid))) 
  sel2=np.where(abs(grid.logz+logOHsun-d.zgrid-d.eupzgrid) == min(abs(grid.logz+logOHsun-d.zgrid-d.eupzgrid)))
  sel3=np.where(abs(grid.logz+logOHsun-d.zgrid+d.edownzgrid) == min(abs(grid.logz+logOHsun-d.zgrid+d.edownzgrid)))
  sel1=sel1[np.sort(grid[sel1].logq)]
  sel2=sel2[np.sort(grid[sel2].logq)]
  sel3=sel3[np.sort(grid[sel3].logq)]

  gridratio=ga/gb
  gridratio0=ga0/gb0
 
  plt.figure()
  plt.plot(grid[sel1].logq, alog10(gridratio[sel1]), yrange=range, xrange=[min(grid.logq)-0.1, max(grid.logq)+0.1])
  # xthick=3, ythick=3, charthick=3, charsize=3, 
  plt.xlabel('log(q)')
  plt.ylabel(title) 
  #thick=3, linestyle=2, xstyle=1, ystyle=1
  
  #plotsym, 0, 1.2, /fill; 
  for i in range (n_elements(auxzarr)-1):
      sel=np.where(grid.logz ==  auxzarr[i]) 
      sel=sel[np.sort(grid[sel].logq)]
      plot.plot(grid[sel].logq, math.log10(gridratio[sel]), color='gray', psym='.')

  for i in range(n_elements(auxzarr0)-1):
     sel=np.where(grid0.logz eq auxzarr0[i]) 
     sel=sel[np.sort(grid0[sel].logq)]
     plt.plot(grid0[sel].logq, math.log10(gridratio0[sel]), color='gray', thick=3;, psym='.')
















  oplot, grid[sel1].logq, alog10(gridratio[sel1]), thick=4, color=cgcolor('red')
  oplot, grid[sel2].logq, alog10(gridratio[sel2]), thick=4, color=cgcolor('red'), linestyle=2
  oplot, grid[sel3].logq, alog10(gridratio[sel3]), thick=4, color=cgcolor('red'), linestyle=2
  if (flaga eq 1 and flagb eq 1) then begin
     dratio=da/db
     edratio=sqrt((eda/db)^2+(da/db^2*edb)^2)/dratio/alog(10.)
     oploterror, [d.qgrid], [alog10(dratio)], [d.eupqgrid], [edratio], psym=3, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.qgrid], [alog10(dratio)], [d.edownqgrid], [edratio], psym=3, thick=4, /lobar , color=cgcolor('black')
  endif else if (flaga eq 2 and flagb eq 1) then begin
     dratio=eda/db
     edratio=eda/(db-edb)-eda/db
     oploterror, [d.qgrid], [alog10(dratio)], [d.eupqgrid], [edratio], psym=8, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.qgrid], [alog10(dratio)], [d.edownqgrid], [0], psym=8, thick=4, /lobar , color=cgcolor('black')
     plotsym, 1, 10, thick=4
     oplot, [d.qgrid], [alog10(dratio)], psym=8, thick=4 , color=cgcolor('black')
   endif else if (flaga eq 1 and flagb eq 2) then begin
     dratio=da/edb
     edratio=da/edb-(da-eda)/edb
     oploterror, [d.qgrid], [alog10(dratio)], [d.eupqgrid], [0], psym=8, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.qgrid], [alog10(dratio)], [d.edownqgrid], [edratio], psym=8, thick=4, /lobar , color=cgcolor('black')
     plotsym, 2, 10, thick=4
     oplot, [d.qgrid], [alog10(dratio)], psym=8, thick=4 , color=cgcolor('black')
  endif
END
 
