import numpy as np


def LogLikRatio (background, signal, N_experiments=10000) :
    
           
    b = background
    s = signal
    s_tot = s.sum()
    
    llr_b_like = []
    llr_sPlusb_like = []

    for k in xrange(N_experiments) :
        N_b = np.random.poisson(lam=b)
        N_sPlusb = np.random.poisson(lam=(s+b))

        llr_b_like.append(2*s_tot - 2*np.dot(N_b,np.log(1+s/b)))
        llr_sPlusb_like.append(2*s_tot - 2*np.dot(N_sPlusb,np.log(1+s/b)))
    
    return llr_b_like, llr_sPlusb_like

def LogLikRatioObserved (background, signals, data) :
    
           
    b = background
    N = data
    llr_data_is_b_like = []
    for s in signals :
        s_tot = s.sum()
    
        llr_data_is_b_like.append(2*s_tot - 2*np.dot(N,np.log(1+s/b))) 
   
    
    return llr_data_is_b_like