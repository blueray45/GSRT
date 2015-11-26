# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:25:08 2015

@author: irnakat
"""
def sh_tf_fix(hl,vs,dn,qs,freq):
    """
    Computes the theoretical SH transfer function of a given 1D model
    model vertical (fixed) wave propagation
    Translated from sh_tf_fix.m written by Poggi Valerio 2009
    
    Version 26 November 2015
    
    """
    from numpy import zeros,exp
    from numpy.linalg import solve
    import numpy as np
    
    # number of layer and frequencies
    nlayer = len(hl)
    fnum = len(freq)
    
    # angular frequencies conversion
    angf = 2.*np.pi*freq
    
    # elastic parameters
    
    # attenuation using complex velocities
    vs = vs*((2.*qs*1j)/(2.*qs*1j-1.))
    
    # modulus of rigidity
    mu = dn*(vs**2)
    
    # slowness parameter    
    ns = 1./vs
    
    # arays initialization
    tf = zeros((fnum,1))
    
    # data vectors
    Ds = zeros((nlayer*2,1))
    Ds[-1] = 1.
    
    # building core matrix
    CORE = zeros((nlayer*2,nlayer*2))
    
    # free surface constraints
    
    CORE[0,0] = 1
    CORE[0,1] = 1
    
    # interfaces constraints
    for nl in range(nlayer-1):
        row = ((nl)*2)+1
        col = ((nl)*2)
        CORE[row+1,col+3] = 1
        CORE[row+1,col+4] = -1
        CORE[row+2,col+3] = -mu[nl+1]*ns[nl+1]
        CORE[row+2,col+4] = -mu[nl+1]*ns[nl+1]
        
    # input constraints
    CORE[-1,-1] = 1
        
    # loop over frequencies
    for nf in range(fnum):
        #----------------------------------------------
        # Interfaces Constraints
        
        for nl in range(nlayer-1):
            row = (nl*2)+1
            col = (nl*2)
            
            expDSA = exp(1j*angf[nf]*ns[nl]*hl[nl])
            expUSA = exp(-1j*angf[nf]*ns[nl]*hl[nl])
            CORE[row+1,col+1] = -expDSA
            CORE[row+1,col+2] = expUSA
            CORE[row+2,col+1] = mu[nl]*ns[nl]*expDSA
            CORE[row+2,col+2] = mu[nl]*ns[nl]*expUSA
            
        # solving linear system
        As = solve(CORE,Ds)
        
        # transfer function
        tf[nf] = (As[1]-As[0])/(2.*As[-1])
        
        return tf