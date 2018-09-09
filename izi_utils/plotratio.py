import numpy as np

def plotratioz( ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, ax, yrange=[-2.0, 2.0], title='log(ratio)'):
    #plot line ratio vs metallicity
        #print flaga, flagb
        auxqarr = np.unique(grid['LOGQ'])
        auxqarr0 = np.unique(grid0['LOGQ'])

        sel1=np.where(abs(grid['LOGQ']-d.qgrid) == min(abs(grid['LOGQ']-d.qgrid))) 
        sel2=np.where(abs(grid['LOGQ']-d.qgrid-d.eupqgrid) == min(abs(grid['LOGQ']-d.qgrid-d.eupqgrid)))
        sel3=np.where(abs(grid['LOGQ']-d.qgrid+d.edownqgrid) == min(abs(grid['LOGQ']-d.qgrid+d.edownqgrid)))

        gridratio=ga/gb
        gridratio = gridratio.flatten()
        gridratio0=ga0/gb0
        gridratio0 = gridratio0.flatten()
        dratio=da/db
        edratio=np.sqrt((eda/db)**2 + (da/db**2*edb)**2)/dratio/np.log(10.)

        sel=np.where(abs(grid['LOGQ']-d.qgrid) == np.min(abs(grid['LOGQ']-d.qgrid)))[0] 
        
        ax.plot(grid['LOGZ'][sel]+logOHsun, np.log10(gridratio[sel]),'--')
        ax.set_ylim(yrange)
        ax.set_xlim(min(grid['LOGZ'])+logOHsun-0.1, max(grid['LOGZ'])+logOHsun+0.1)
        ax.set_xlabel('12+log(O/H)')
        ax.set_ylabel(title)
        
        for i in range(len(auxqarr0)):
            sel=np.where(grid0['LOGQ'] == auxqarr0[i]) 
            ax.plot(grid0[sel]['LOGZ']+logOHsun, np.log10(gridratio0[sel]),color='gray',linewidth = 1)

        ax.plot(grid[sel1]['LOGZ']+logOHsun, np.log10(gridratio[sel1]), linewidth=2, color='red')
        ax.plot(grid[sel2]['LOGZ']+logOHsun, np.log10(gridratio[sel2]),'--', linewidth=2, color='red')
        ax.plot(grid[sel3]['LOGZ']+logOHsun, np.log10(gridratio[sel3]),'--', linewidth=2, color='red')

        if ((flaga == 1) and (flagb == 1)):
            dratio=da/db
            edratio1=np.sqrt((eda/db)**2+(da/db**2*edb)**2)/dratio/np.log10(10.)
            ax.errorbar([d.Zgrid], np.log10(dratio), yerr = edratio1, fmt = 'o', xerr = [[d.edownZgrid]*len([d.Zgrid]), [d.eupZgrid]*len([d.Zgrid])], linewidth=2, color='k')

        elif (flaga == 2 and flagb == 1):
            dratio=eda/db
            edratio=eda/(db-edb)-eda/db
            ax.errorbar([d.Zgrid], np.log10(dratio), yerr = edratio1, fmt = 'o', xerr = [[d.edownZgrid]*len([d.Zgrid]), [d.eupZgrid]*len([d.Zgrid])], linewidth=2, color='k')
            symbols = [u'\u2193']
            ax.plot( d.Zgrid, np.log10(dratio), symbols,size = 10, linewidth=2 , color='k')

        elif (flaga == 1 and flagb == 2):
            dratio=da/edb
            edratio=da/edb-(da-eda)/edb
            ax.errorbar([d.Zgrid], np.log10(dratio), yerr = edratio1, fmt = 'o', xerr = [[d.edownZgrid]*len([d.Zgrid]), [d.eupZgrid]*len([d.Zgrid])], linewidth=2, color='k')
            symbols = [u'\u2191']
            ax.plot( d.Zgrid, np.log10(dratio), symbols,size = 10, linewidth=4 , color='k')
            
        

def  plotratioq (ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, ax, yrange=[-2.0,2.0], title='log(ratio)'):
# plot a line ratio vs the ionization parameter
        auxzarr = np.unique(grid['LOGZ'])
        auxzarr0 = np.unique(grid0['LOGZ'])

        sel1=np.where(abs(grid['LOGZ']+logOHsun-d.Zgrid) == min(abs(grid['LOGZ']+logOHsun-d.Zgrid)))[0] 
        sel2=np.where(abs(grid['LOGZ']+logOHsun-d.Zgrid-d.eupZgrid) == min(abs(grid['LOGZ']+logOHsun-d.Zgrid-d.eupZgrid)))[0]
        sel3=np.where(abs(grid['LOGZ']+logOHsun-d.Zgrid+d.edownZgrid) == min(abs(grid['LOGZ']+logOHsun-d.Zgrid+d.edownZgrid)))[0]
        sel1=sel1[np.argsort(grid[sel1]['LOGQ'])]
        sel2=sel2[np.argsort(grid[sel2]['LOGQ'])]
        sel3=sel3[np.argsort(grid[sel3]['LOGQ'])]

        gridratio=ga/gb
        gridratio = gridratio.flatten()
        gridratio0=ga0/gb0
        gridratio0 = gridratio0.flatten()

        ax.plot(grid['LOGQ'][sel1], np.log10(gridratio[sel1]), color = 'red', linewidth = 2)
        ax.set_ylim(yrange)
        ax.set_xlim(min(grid['LOGQ'])-0.1, max(grid['LOGQ'])+0.1)
        ax.set_xlabel('log(q)')
        ax.set_ylabel(title)

        for i in range(len(auxzarr0)):
            sel=np.where(grid0['LOGZ'] == auxzarr0[i])[0]
            sel=sel[np.argsort(grid0[sel]['LOGQ'])]
            ax.plot(grid0[sel]['LOGQ'], np.log10(gridratio0[sel]), color='gray', linewidth = 1)

        ax.plot(grid[sel2]['LOGQ'], np.log10(gridratio[sel2]),'--', linewidth=2, color='red')
        ax.plot(grid[sel3]['LOGQ'], np.log10(gridratio[sel3]),'--', linewidth=2, color='red')

        if ((flaga == 1) and (flagb == 1)):
            dratio=da/db
            edratio1=np.sqrt((eda/db)**2+(da/db**2*edb)**2)/dratio/np.log10(10.)
            ax.errorbar([d.qgrid], np.log10(dratio), yerr = edratio1, fmt = 'o', xerr = [[d.edownqgrid]*len([d.qgrid]), [d.eupqgrid]*len([d.qgrid])], linewidth=2, color='k')

        elif (flaga == 2 and flagb == 1):
            dratio=eda/db
            edratio=eda/(db-edb)-eda/db
            ax.errorbar([d.qgrid], np.log10(dratio), yerr = edratio1, fmt = 'o', xerr = [[d.edownqgrid]*len([d.qgrid]), [d.eupqgrid]*len([d.qgrid])], linewidth=2, color='k')
            symbols = [u'\u2193']
            ax.plot( d.qgrid, np.log10(dratio), symbols,size = 10, linewidth=2 , color='k')

        elif (flaga == 1 and flagb == 2):
            dratio=da/edb
            edratio=da/edb-(da-eda)/edb
            ax.errorbar([d.qgrid], np.log10(dratio), yerr = edratio1, fmt = 'o', xerr = [[d.edownqgrid]*len([d.qgrid]), [d.eupqgrid]*len([d.qgrid])], linewidth=2, color='k')
            symbols = [u'\u2191']
            ax.plot( d.qgrid, np.log10(dratio), symbols,size = 10, linewidth=4 , color='k')
