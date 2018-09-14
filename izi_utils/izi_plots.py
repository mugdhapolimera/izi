import numpy as np
import matplotlib.pyplot as plt
import os
import sys
#sys.path.append('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/izi_utils/')
sys.path.append('C:\Users\mugdhapolimera\github\izi\izi_utils')
from plotratio import plotratioz 
from plotratio import plotratioq

#os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/izi_plots/')
os.chdir('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\Codes\izi_plots')
def izi_pdf(d, grid, postz, postq, post, like, plot_flag):
    
    #CALCULATE 1-,2-,3- SIGMA VALUES
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

    post1sig=(sortpost[np.where(sumpost >= 0.683)])
    post2sig=(sortpost[np.where(sumpost >= 0.955)])
    post3sig=(sortpost[np.where(sumpost >= 0.997)])

    if post1sig.shape == 0:
        post1sig = post1sig[0]
    else:
        post1sig = 0
    
    if post3sig.shape == 0:
        post3sig = post3sig[0]
    else:
        post3sig = 0
    if post2sig.shape == 0:
        post2sig = post2sig[0]
    else:
        post2sig = 0
    
    like1sig=(sortlike[np.where(sumlike >= 0.683)])
    like2sig=(sortlike[np.where(sumlike >= 0.955)])
    like3sig=(sortlike[np.where(sumlike >= 0.997)])
    
    if like1sig.shape == 0:
        like1sig = like1sig[0]
    else:
        like1sig = 0
    
    if like2sig.shape == 0:
        like2sig = like2sig[0]
    else:
        like2sig = 0

    if like3sig.shape == 0:
        like3sig = like3sig[0]
    else:
        like3sig = 0
    


    logOHsun = grid['LOGOHSUN'][0]
    zarr = np.unique(grid['LOGZ'])
    qarr = np.unique(grid['LOGQ'])
    
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(50,10))
    col_range = [0,1]
    col_range[0] = 0*np.max(post[np.where(np.isfinite(np.log10(post)))[0]])
    col_range[1] = 1*np.max(post[np.where(np.isfinite(np.log10(post)))[0]])

    ncolors=256
    levels = col_range[0] + (col_range[1]-col_range[0])/(ncolors-1.0)*np.arange(ncolors)
    if (np.sum(levels) == 0.):
        return
    cmap = 'gist_earth'
    ax1.tricontourf(grid['LOGZ']+logOHsun, grid['LOGQ'], post,levels = levels)#, cmap = cmap)
    ax1.tricontour(grid[np.where(np.isfinite(np.log10(post)))[0]]['LOGZ']+logOHsun, 
                   grid[np.where(np.isfinite(np.log10(post)))[0]]['LOGQ'], 
                   post[np.where(np.isfinite(np.log10(post)))[0]], levels=[post1sig], colors = 'cyan')
    ax1.plot([logOHsun, logOHsun], [0,1e2], linestyle='--', linewidth=5, color='orange')
    ax1.plot([logOHsun-1.0, logOHsun-1.0], [0,1e2], linestyle='--', linewidth=5, color='orange')
    ax1.text(logOHsun+0.05, min(grid['LOGQ'])+0.1, r'$Z_{\odot}$', color = 'orange', size = 30)
    ax1.text(logOHsun-1+0.05, min(grid['LOGQ'])+0.1, r'$0.1Z_{\odot}$', color = 'orange', size = 30)
    ax1.plot(grid['LOGZ']+logOHsun, grid['LOGQ'], '.', color = 'white')
    ax1.plot(d.Zgrid, d.qgrid, marker = 'o', color='orange', markersize = 10)

    ax1.set_xlabel('12+log(O/H)')
    ax1.set_ylabel('log(q)')
    ax1.set_xlim([np.min(grid['LOGZ'])+logOHsun, np.max(grid['LOGZ'])+logOHsun])
    ax1.set_ylim([np.min(grid['LOGQ']), np.max(grid['LOGQ'])])
    ax1.set_title('P(Z,q)d(log Z)d(log q)')
    
    sel=np.where(np.isfinite(postz))
    ax2.plot(zarr[sel]+logOHsun, postz[sel], linewidth = '4',color = 'black')
    ax2.plot([d.Zgrid, d.Zgrid], [0,1.1*np.max(postz)], linestyle='-', linewidth=2, color='red')
    ax2.plot([d.Zgridmarmod, d.Zgridmarmod], [0,1.1*np.max(postz)], linestyle='-', linewidth=2, color='blue')
    ax2.plot([d.Zgridmarmean, d.Zgridmarmean], [0,1.1*np.max(postz)], linestyle='-', linewidth=2, color='green')
    ax2.plot([d.Zgrid-d.edownZgrid, d.Zgrid-d.edownZgrid], [0,1.1*np.max(postz)],linestyle=':', linewidth=1, color='red')
    ax2.plot([d.Zgrid+d.eupZgrid, d.Zgrid+d.eupZgrid], [0,1.1*np.max(postz)], linestyle=':', linewidth=1, color='red')
    ax2.plot([logOHsun, logOHsun], [0,1.1*np.max(postz)], linestyle='--', linewidth=2)
    ax2.plot([logOHsun-1, logOHsun-1], [0,1.1*np.max(postz)], linestyle='--', linewidth=2)

    ax2.set_xlabel('12+log(O/H)')
    ax2.set_ylabel('P(Z)d(log Z)')
    ax2.set_xlim([np.min(zarr)+logOHsun, np.max(zarr)+logOHsun])
    ax2.set_ylim([0, 1.1*np.max(postz)])

    sel=np.where(np.isfinite(postq))
    ax3.plot(qarr[sel], postq[sel],linewidth=4, color = 'black')
    ax3.plot([d.qgrid, d.qgrid], [0,1e2], linestyle='-', linewidth=4, color='red')
    ax3.plot([d.qgridmarmod, d.qgridmarmod], [0,1e2], linestyle='-', linewidth=2, color='blue')
    ax3.plot([d.qgridmarmean, d.qgridmarmean], [0,1e2], linestyle='-', linewidth=2, color='green')
    ax3.plot([d.qgrid-d.edownqgrid, d.qgrid-d.edownqgrid], [0,1e2], linestyle='--', linewidth=2, color='red')
    ax3.plot([d.qgrid+d.eupqgrid, d.qgrid+d.eupqgrid], [0,1e2], linestyle='--', linewidth=2, color='red')
    ax3.text(min(grid['LOGQ'])+1.1, 1.0*max(postq), 'Joint Mode', color='red', size=30)
    ax3.text(min(grid['LOGQ'])+1.1, 0.925*max(postq), 'Marg Mode', color='blue', size = 30)
    ax3.text(min(grid['LOGQ'])+1.1, 0.85*max(postq), 'Marg. Mean', color='green', size = 30)

    ax3.set_xlabel('log(q)')
    ax3.set_ylabel('P(q)d(log q)')
    ax3.set_xlim([np.min(qarr), np.max(qarr)])
    ax3.set_ylim([0, 1.1*max(postq)])
    plt.savefig(str(d['name'])+'_pdf.png')

    if plot_flag:
        plt.show(block=True)
    #plt.hold(True)
    else:
        plt.close()

    
def zratios_plots (grid, grid0, d, flag0, plot_flag):
    logOHsun = grid['LOGOHSUN'][0]
    
    #======== PLOTS VS METALLICITY ============

    fig2, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4,figsize=(50,10))



    #======= R23 ========
    #define line ratio

    inoiii5007 = np.where('oiii5007' == grid['ID'][0])[0]
    inoii3726 = np.where('oii3726' == grid['ID'][0])[0]
    inhbeta = np.where('hbeta' == grid['ID'][0])[0]

    ga = np.array(grid['FLUX'])[:,inoiii5007] + np.array(grid['FLUX'])[:,inoii3726]
    gb = np.array(grid['FLUX'])[:,inhbeta]
    ga0 = np.array(grid0['FLUX'])[:,inoiii5007] + np.array(grid0['FLUX'])[:,inoii3726]
    gb0 = np.array(grid0['FLUX'])[:,inhbeta]
    da = d.flux[inoiii5007] + d.flux[inoii3726]
    db = d.flux[inhbeta]
    eda = np.sqrt(d.error[inoiii5007]**2 + d.error[inoii3726]**2)
    edb = d.error[inhbeta]     
    flaga = 1
    if (flag0[inoiii5007] == 2 or flag0[inoii3726] == 2):
        flaga=2
        da=-666
        eda=np.sum(([d.flux[inoiii5007], d.flux[inoii3726]])[np.where([flag0[inoiii5007],flag0[inoii3726]] == 1)[0]]) + np.sum(([d.error[inoiii5007],d.error[inoii3726]])[np.where([flag0[inoiii5007],flag0[inoii3726]] == 2)[0]])

    if (flag0[inoiii5007] == 0 or flag0[inoii3726] == 0):
        flaga = 0
        
    flagb = flag0[inhbeta]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax1, title='log(R23)', yrange=[-1.5,1.5])

    # ======= N2O2 ========
    # define line ratio

    innii6584 = np.where('nii6584' == grid['ID'][0])[0]
    ga=np.array(grid['FLUX'])[:,innii6584]
    gb=np.array(grid['FLUX'])[:,inoii3726]
    ga0=np.array(grid0['FLUX'])[:,innii6584]
    gb0=np.array(grid0['FLUX'])[:,inoii3726]
    da=d.flux[innii6584]
    db=d.flux[inoii3726]
    eda=d.error[innii6584]
    edb=d.error[inoii3726]
    flaga=flag0[innii6584]
    flagb=flag0[inoii3726]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax2, title='log(N2O2)', yrange=[-2,1])

    #======= N2 ========
    #define line ratio
    inhalpha = np.where('halpha' == grid['ID'][0])[0]

    ga=np.array(grid['FLUX'])[:,innii6584]
    gb=np.array(grid['FLUX'])[:,inhalpha]
    ga0=np.array(grid0['FLUX'])[:,innii6584]
    gb0=np.array(grid0['FLUX'])[:,inhalpha]
    da=d.flux[innii6584]
    db=d.flux[inhalpha]
    eda=d.error[innii6584]
    edb=d.error[inhalpha]
    flaga=flag0[innii6584]
    flagb=flag0[inhalpha]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax3, title='log(N2)', yrange=[-3.0,0.0])


    # ======= O3N2 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,inoiii5007]
    gb=np.array(grid['FLUX'])[:,innii6584]
    ga0=np.array(grid0['FLUX'])[:,inoiii5007]
    gb0=np.array(grid0['FLUX'])[:,innii6584]
    da=d.flux[inoiii5007]
    db=d.flux[innii6584]
    eda=d.error[inoiii5007]
    edb=d.error[innii6584]
    flaga=flag0[inoiii5007]
    flagb=flag0[innii6584]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax4, title='log(O3N2)', yrange=[-0.5,2.5])

    # ======= O3O2 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,inoiii5007]
    gb=np.array(grid['FLUX'])[:,inoii3726]
    ga0=np.array(grid0['FLUX'])[:,inoiii5007]
    gb0=np.array(grid0['FLUX'])[:,inoii3726]
    da=d.flux[inoiii5007]
    db=d.flux[inoii3726]
    eda=d.error[inoiii5007]
    edb=d.error[inoii3726]
    flaga=flag0[inoiii5007]
    flagb=flag0[inoii3726]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax5, title='log(O3O2)', yrange=[-1.5,1.5])

    # ======= R3 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,inoiii5007]
    gb=np.array(grid['FLUX'])[:,inhbeta]
    ga0=np.array(grid0['FLUX'])[:,inoiii5007]
    gb0=np.array(grid0['FLUX'])[:,inhbeta]
    da=d.flux[inoiii5007]
    db=d.flux[inhbeta]
    eda=d.error[inoiii5007]
    edb=d.error[inhbeta]
    flaga=flag0[inoiii5007]
    flagb=flag0[inhbeta]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax6, title='log(R3)', yrange=[-2,1])

    #======= N2S2 ========
    # define line ratio
    insii6717 = np.where('sii6717' == grid['ID'][0])[0]
    insii6731 = np.where('sii6731' == grid['ID'][0])[0]

    ga=np.array(grid['FLUX'])[:,innii6584]
    gb=np.array(grid['FLUX'])[:,insii6717]+np.array(grid['FLUX'])[:,insii6731]
    ga0=np.array(grid0['FLUX'])[:,innii6584]
    gb0=np.array(grid0['FLUX'])[:,insii6717]+np.array(grid0['FLUX'])[:,insii6731]
    da=d.flux[innii6584]
    db=d.flux[insii6717]
    eda=d.error[innii6584]
    edb=d.error[insii6717]
    flaga=flag0[innii6584]
    flagb=flag0[insii6717]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax7, title='log(N2S2)', yrange=[-1.5,1.5])

    # ======= S2 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,insii6717]+np.array(grid['FLUX'])[:,insii6731]
    gb=np.array(grid['FLUX'])[:,inhalpha]
    ga0=np.array(grid0['FLUX'])[:,insii6717]+np.array(grid0['FLUX'])[:,insii6731]
    gb0=np.array(grid0['FLUX'])[:,inhalpha]
    da=d.flux[insii6717]
    db=d.flux[inhalpha]
    eda=d.error[insii6717]
    edb=d.error[inhalpha]
    flaga=flag0[insii6717]
    flagb=flag0[inhalpha]
    plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, ax8, title='log(S2)', yrange=[-3.0,0.0])

    plt.savefig(str(d['name'])+'_zratios.png')

    if plot_flag:
        plt.show(block=True)
    #plt.hold(True)
        
    else:
        plt.close()

    
    
def qratios_plots(grid, grid0, d, flag0, plot_flag):
    logOHsun = grid['LOGOHSUN'][0]
    
    #============= PLOTS VERSUS IONIZATION PARAMETER ===========


    fig3, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, figsize=(70,10))

    #======= R23 ========
    #define line ratio

    inoiii5007 = np.where('oiii5007' == grid['ID'][0])[0]
    inoii3726 = np.where('oii3726' == grid['ID'][0])[0]
    inhbeta = np.where('hbeta' == grid['ID'][0])[0]

    ga = np.array(grid['FLUX'])[:,inoiii5007] + np.array(grid['FLUX'])[:,inoii3726]
    gb = np.array(grid['FLUX'])[:,inhbeta]
    ga0 = np.array(grid0['FLUX'])[:,inoiii5007] + np.array(grid0['FLUX'])[:,inoii3726]
    gb0 = np.array(grid0['FLUX'])[:,inhbeta]
    da = d.flux[inoiii5007] + d.flux[inoii3726]
    db = d.flux[inhbeta]
    eda = np.sqrt(d.error[inoiii5007]**2 + d.error[inoii3726]**2)
    edb = d.error[inhbeta]     
    flaga = 1
    if (flag0[inoiii5007] == 2 or flag0[inoii3726] == 2):
        flaga=2
        da=-666
        eda=total(([d.flux[inoiii5007], d.flux[inoii3726]])[np.where([flag0[inoiii5007],flag0[inoii3726]] == 1)[0]]) + np.sum(([d.error[inoiii5007],d.error[inoii3726]])[np.where([flag0[inoiii5007],flag0[inoii3726]] == 2)[0]])

    if (flag0[inoiii5007] == 0 or flag0[inoii3726] == 0): 
        flaga = 0
    flagb = flag0[inhbeta]

    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax1, title='log(R23)', yrange=[-1.5,1.5])

    # ======= N2O2 ========
    # define line ratio

    innii6584 = np.where('nii6584' == grid['ID'][0])[0]
    ga=np.array(grid['FLUX'])[:,innii6584]
    gb=np.array(grid['FLUX'])[:,inoii3726]
    ga0=np.array(grid0['FLUX'])[:,innii6584]
    gb0=np.array(grid0['FLUX'])[:,inoii3726]
    da=d.flux[innii6584]
    db=d.flux[inoii3726]
    eda=d.error[innii6584]
    edb=d.error[inoii3726]
    flaga=flag0[innii6584]
    flagb=flag0[inoii3726]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax2, title='log(N2O2)', yrange=[-2,1])

    #======= N2 ========
    #define line ratio
    inhalpha = np.where('halpha' == grid['ID'][0])[0]

    ga=np.array(grid['FLUX'])[:,innii6584]
    gb=np.array(grid['FLUX'])[:,inhalpha]
    ga0=np.array(grid0['FLUX'])[:,innii6584]
    gb0=np.array(grid0['FLUX'])[:,inhalpha]
    da=d.flux[innii6584]
    db=d.flux[inhalpha]
    eda=d.error[innii6584]
    edb=d.error[inhalpha]
    flaga=flag0[innii6584]
    flagb=flag0[inhalpha]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax3, title='log(N2)', yrange=[-3.0,0.0])


    # ======= O3N2 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,inoiii5007]
    gb=np.array(grid['FLUX'])[:,innii6584]
    ga0=np.array(grid0['FLUX'])[:,inoiii5007]
    gb0=np.array(grid0['FLUX'])[:,innii6584]
    da=d.flux[inoiii5007]
    db=d.flux[innii6584]
    eda=d.error[inoiii5007]
    edb=d.error[innii6584]
    flaga=flag0[inoiii5007]
    flagb=flag0[innii6584]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax4, title='log(O3N2)', yrange=[-2.0,2.0])

    # ======= O3O2 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,inoiii5007]
    gb=np.array(grid['FLUX'])[:,inoii3726]
    ga0=np.array(grid0['FLUX'])[:,inoiii5007]
    gb0=np.array(grid0['FLUX'])[:,inoii3726]
    da=d.flux[inoiii5007]
    db=d.flux[inoii3726]
    eda=d.error[inoiii5007]
    edb=d.error[inoii3726]
    flaga=flag0[inoiii5007]
    flagb=flag0[inoii3726]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax5, title='log(O3O2)', yrange=[-1.5,1.5])

    # ======= R3 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,inoiii5007]
    gb=np.array(grid['FLUX'])[:,inhbeta]
    ga0=np.array(grid0['FLUX'])[:,inoiii5007]
    gb0=np.array(grid0['FLUX'])[:,inhbeta]
    da=d.flux[inoiii5007]
    db=d.flux[inhbeta]
    eda=d.error[inoiii5007]
    edb=d.error[inhbeta]
    flaga=flag0[inoiii5007]
    flagb=flag0[inhbeta]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax6, title='log(R3)', yrange=[-2,1])

    #======= N2S2 ========
    # define line ratio
    insii6717 = np.where('sii6717' == grid['ID'][0])[0]
    insii6731 = np.where('sii6731' == grid['ID'][0])[0]

    ga=np.array(grid['FLUX'])[:,innii6584]
    gb=np.array(grid['FLUX'])[:,insii6717]+np.array(grid['FLUX'])[:,insii6731]
    ga0=np.array(grid0['FLUX'])[:,innii6584]
    gb0=np.array(grid0['FLUX'])[:,insii6717]+np.array(grid0['FLUX'])[:,insii6731]
    da=d.flux[innii6584]
    db=d.flux[insii6717]
    eda=d.error[innii6584]
    edb=d.error[insii6717]
    flaga=flag0[innii6584]
    flagb=flag0[insii6717]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, ax7, title='log(N2S2)', yrange=[-1.5,1.5])

    # ======= S2 ========
    # define line ratio
    ga=np.array(grid['FLUX'])[:,insii6717]+np.array(grid['FLUX'])[:,insii6731]
    gb=np.array(grid['FLUX'])[:,inhalpha]
    ga0=np.array(grid0['FLUX'])[:,insii6717]+np.array(grid0['FLUX'])[:,insii6731]
    gb0=np.array(grid0['FLUX'])[:,inhalpha]
    da=d.flux[insii6717]
    db=d.flux[inhalpha]
    eda=d.error[insii6717]
    edb=d.error[inhalpha]
    flaga=flag0[insii6717]
    flagb=flag0[inhalpha]
    plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, ax8, title='log(S2)', yrange=[-3.0,0.0])

    plt.savefig(str(d['name'])+'_qratios.png')

    if plot_flag:
        plt.show(block=True)
    #plt.hold(True)
    
    else:
        plt.close()

