import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats


def BkgSigHistos (background, signals, data, variable_binning,x_label,savepath=None) :
    
    bkg = background
    sigModels = signals
    binning = variable_binning
    x_name, x_unit = x_label
    data_hist = data
    
    binw = binning[1]-binning[0]
    #binw = np.array([10,10,10,5,5,5,5,5,5,10,20])
    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True,figsize=(8,10))

    m_H = [85,90,95]
    
    for i in xrange(3) :

        ax = axs[i]

        ax.bar(binning[:-1]+binw/2.,bkg[i],width=binw,label='background',color='yellow',edgecolor='k',linewidth=0.1)

        ax.bar(binning[:-1]+binw/2.,sigModels[i],width=binw,linewidth=0.1,
               label=r'H signal ($m_\mathrm{H}= $'+str(m_H[i])+' GeV)',color='red',bottom=bkg[i],edgecolor='k')
        
        ax.errorbar(x=binning[:-1]+binw/2.,y=data_hist[i], xerr=binw/2., yerr=np.sqrt(data_hist[i]),
                    fmt='o', color='k',label='data',linewidth=1)
        ax.legend()
        ax.set_ylabel('Events / '+ str(round(binw,1))+' '+x_unit)
        if (i == 2) :
            ax.set_xlabel(x_name+' ['+x_unit+']')
        plt.tight_layout()

    fig.subplots_adjust(hspace=0)
    
    if (savepath != None) :
        #plt.savefig('plots/test')
        plt.savefig(savepath)
    else :
        plt.show()

        
        
        
#----------------------------------------------------------------------------------------------
def LogLikRatioPlots(arrays,obs,Nbins=30,savepath=None) :
 
    
    fig, axs = plt.subplots(nrows=3, ncols=1,figsize=(8,10))
    m_H = [85,90,95]

    CLlist = []
    
    QuantileList_b = []
    QuantileList_sPlusb = []
    
    for i in xrange(3) :
        ax = axs[i]

        llr_b, llr_sPlusb = arrays[i]
        norm = len(llr_b)
        binning = np.linspace(np.minimum(llr_b,llr_sPlusb).min(),np.maximum(llr_b,llr_sPlusb).max(),Nbins)
        
        llr_b_hist = 1.*np.histogram(llr_b,bins=binning)[0]/norm
        QuantileList_b.append(GetQuantiles(llr_b_hist,binning)) 
        
        
        pos =  np.where(binning <= obs[i])[0][-1]
        print pos
        OneMinusCLb =  sum(llr_b_hist[:pos]) 
        
        llr_sPlusb_hist = 1.*np.histogram(llr_sPlusb,bins=binning)[0]/norm
        QuantileList_sPlusb.append(GetQuantiles(llr_sPlusb_hist,binning)) 
        
        CLsPlusb =  sum(llr_sPlusb_hist[pos:])
        CLlist.append([OneMinusCLb, CLsPlusb])
       
        ax.step(x=binning[:-1],y=llr_b_hist,color='blue',label='bkg-like')
        ax.step(x=binning[:-1],y=llr_sPlusb_hist,color='red',label='sig+bkg-like')

        x1 = binning[binning<=obs[i]]
        x2 = binning[binning>obs[i]]
        
        #ax.fill_between(x1,llr_b_hist[:len(x1)],color='red',alpha=0.5,interpolate=True)
        #ax.fill_between(x2,llr_sPlusb_hist[-len(x2):],color='blue',alpha=0.5,interpolate=True)
        
        ax.set_xlabel(r'$-2 \ln (Q)$')
        ax.set_ylabel('p.d.f.')
        ax.set_title('signal model ' + r'($m_\mathrm{H} = $'+str(m_H[i])+' GeV)')
        ax.axvline(obs[i],label='observed',color='k')
        ax.legend()
    
    plt.tight_layout()
    if (savepath != None) :
        #plt.savefig('plots/test')
        plt.savefig(savepath)
    else :
        plt.show()      
    
    
    
    return CLlist, QuantileList_b, QuantileList_sPlusb



def GetQuantiles (hist,binning) :
    # note that the histogram must be normalized
    cumulative = np.cumsum(hist)
    
    twoSigmaLeft = 0.023
    oneSigmaLeft = 0.16
    median = 0.5
    oneSigmaRight = 1.-0.16  
    twoSigmaRight = 1.-0.023
    
    
    TwoSigmaLeft = binning[np.where(cumulative <= twoSigmaLeft)[0][-1]] 
    OneSigmaLeft = binning[np.where(cumulative <= oneSigmaLeft)[0][-1]] 
    Median = binning[np.where(cumulative <= median)[0][-1]]
    OneSigmaRight = binning[np.where(cumulative < oneSigmaRight)[0][-1]]  
    TwoSigmaRight = binning[np.where(cumulative < twoSigmaRight)[0][-1]] 
    
        
    return [Median,[OneSigmaLeft,OneSigmaRight],[TwoSigmaLeft,TwoSigmaRight]]


"""
OLD BACK-UP STUFF
"""
def BkgSigHistos_OLD (background, signals, data, variable_binning,x_label,savepath=None) :
    
    bkg = background
    sigModels = signals
    binning = variable_binning
    x_name, x_unit = x_label
    data_hist = data
    
    binw = binning[1]-binning[0]
    #binw = np.array([10,10,10,5,5,5,5,5,5,10,20])
    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True,figsize=(8,10))

    m_H = [85,90,95]
    
    for i in xrange(3) :

        ax = axs[i]

        ax.bar(binning[:-1]+binw/2.,bkg,width=binw,label='background',color='yellow',edgecolor='k',linewidth=0.1)

        ax.bar(binning[:-1]+binw/2.,sigModels[i],width=binw,linewidth=0.1,
               label=r'H signal ($m_\mathrm{H}= $'+str(m_H[i])+' GeV)',color='red',bottom=bkg,edgecolor='k')
        
        ax.errorbar(x=binning[:-1]+binw/2.,y=data_hist, xerr=binw/2., yerr=np.sqrt(data_hist),
                    fmt='o', color='k',label='data',linewidth=1)
        ax.legend()
        ax.set_ylabel('Events / '+ str(round(binw,1))+' '+x_unit)
        if (i == 2) :
            ax.set_xlabel(x_name+' ['+x_unit+']')
        plt.tight_layout()

    fig.subplots_adjust(hspace=0)
    
    if (savepath != None) :
        #plt.savefig('plots/test')
        plt.savefig(savepath)
    else :
        plt.show()
