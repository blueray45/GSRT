# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 20:04:33 2015
Program to calculate transfer function using several different methods

Theodosius Marwan Irnaka
marwan.irnaka@gmail.com
www.github.com/GSRT
"""
import commondata as cd
import numpy as np
from copy import deepcopy

class TFCalculator:
    def __init__(self,data,freq=None):
        self.mode = data['mode']                  # calculation mode
        self.modeID = cd.mode.index(self.mode)
        self.ntf = data['ntf']                   # number of transfer function pairs
        self.tfpair = np.array(data['tfPair'])         # list of transfer function pairs
        self.nlayer = data['nlayer']                # number of soil layer
        self.hl = np.array(data['hl'])          # thickness of each layer
        self.vs = np.array(data['vs'])          # shear wave velocity
        self.dn = np.array(data['dn'])          # density
        self.qs = np.array(data['qs'])          # quality factor
        self.slayer = data['sourceloc']
        self.iang = data['iang']

        if self.modeID>4:                
            self.vp = np.array(data['vp'])
            self.qp = np.array(data['qp'])
            self.comp = np.array(data['comp'])
        if freq==None:
            self.freq = np.linspace(1.,50.,1024)
        else:
            self.freq = freq
        
    def tf_kramer286_sh(self):
        """
        Calculate transfer function using direct and simple implementation from Kramer Book page 286
        
        freq   : frequency array
        vs  : Shear wave velocity
        dn : density
        hl   : thickness of layers
        qs  : Quality factor of Shear wave velocity
        lt  : index of top layer
        lb  : index of bottom layer
        """
        # checking calculation validity
        #if self.mode != cd.mode[0] and self.mode != cd.mode[1]:
        #    raise ValueError("Calculation mode is not supported! Use another moethod!")
        Gc = np.zeros((len(self.hl),len(self.freq)),dtype='complex128')
        A = np.zeros_like(Gc)
        B = np.zeros_like(Gc)
        eta = np.zeros_like(self.qs)
        amp = np.zeros((len(self.freq)),dtype='complex128')
        
        Yi = np.complex(0.,1.)
        A[0,:] = np.complex(1.,0.)
        B[0,:] = np.complex(1.,0.)
        ai = 0.0
        # calculating complex modulus
        for i in range(len(self.hl)):
            eta[i] = 1./(2.*self.qs[i])
            for j in range(len(self.freq)):
                Gc[i,j] = self.dn[i]*self.vs[i]**2 * np.complex(1.,2.*eta[i])
        # calculating propagator terms (A and B)
        for j in range(len(self.freq)):
            for i in range(len(self.hl)-1):
                fc = np.complex(self.freq[j],-ai)
                alphaz = np.sqrt((self.dn[i]*Gc[i,j])/(self.dn[i+1]*Gc[i+1,j]))
                alpha = 2.*np.pi*np.sqrt(self.dn[i]/Gc[i,j])*fc
                A[i+1,j] = 0.5*A[i,j]*(1.0+alphaz) * np.exp(Yi*alpha*self.hl[i]) + \
                    0.5*B[i,j]*(1.-alphaz)*np.exp(-Yi*alpha*self.hl[i])
                B[i+1,j] = 0.5*A[i,j]*(1.0-alphaz) * np.exp(Yi*alpha*self.hl[i]) + \
                    0.5*B[i,j]*(1.+alphaz)*np.exp(-Yi*alpha*self.hl[i])
        # calculating transfer function
        self.tf = []
        for j in range(self.ntf):
            amp = np.zeros((len(self.freq)),dtype='complex128')
            vtf = np.zeros((len(self.freq)),dtype='complex128')
            for i in range(len(self.freq)):
                amp[i] = (A[self.tfpair[j][0],i]+B[self.tfpair[j][0],i])/(2.*A[self.tfpair[j][1],i])
                #amp[i] = (A[self.tfpair[j][0],i]+B[self.tfpair[j][0],i])/(A[self.tfpair[j][1],i]+B[self.tfpair[j][1],i])
            self.tf.append(amp)
            self.tf.append(vtf)
        return self.tf
        
    def tf_knopoff_sh(self):
        """
        Calculate Transfer function for SH-wave vertical incident using knopoff method (1964)
        modified from matlab code created by Poggi Valerio (2009)
        This method tries to solve huge block of matrix. Not very efficient.
        """
        
        from numpy import zeros,exp
        from numpy.linalg import solve
        import numpy as np
        
        # uniforming variable
        mode = self.mode
        ntf = self.ntf
        tfpair = self.tfpair
        nlayer = self.nlayer
        hl = self.hl
        vs = self.vs
        dn = self.dn
        qs = self.qs
        freq = self.freq
        
        # checking calculation validity
        #if mode != cd.mode[0] and mode!=cd.mode[1]:
        #    raise ValueError("Calculation mode is not supported! Use another method!")
        
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
        CORE = zeros((nlayer*2,nlayer*2),dtype='complex128')
        
        # free surface constraints
        
        CORE[0,0] = 1
        CORE[0,1] = 1
        
        # interfaces constraints
        for nl in range(nlayer-1):
            row = ((nl)*2)+1
            col = ((nl)*2)
            CORE[row+0,col+2] = 1
            CORE[row+0,col+3] = -1
            CORE[row+1,col+2] = -mu[nl+1]*ns[nl+1]
            CORE[row+1,col+3] = -mu[nl+1]*ns[nl+1]
            
        # input constraints
        CORE[-1,-1] = 1
            
        # loop over frequencies and number of tfpair
        self.tf = []
        for tfp in range(ntf):
            hft = np.zeros((fnum),dtype='complex128')
            vft = np.zeros((fnum),dtype='complex128')
            for nf in range(fnum):
                #----------------------------------------------
                # Interfaces Constraints
                
                for nl in range(nlayer-1):
                    row = (nl*2)+1
                    col = (nl*2)
                    
                    expDSA = exp(1j*angf[nf]*ns[nl]*hl[nl])
                    expUSA = exp(-1j*angf[nf]*ns[nl]*hl[nl])
                    CORE[row+0,col+0] = -expDSA
                    CORE[row+0,col+1] = expUSA
                    CORE[row+1,col+0] = mu[nl]*ns[nl]*expDSA
                    CORE[row+1,col+1] = mu[nl]*ns[nl]*expUSA
                    
                # solving linear system
                As = solve(CORE,Ds)
                
                # transfer function
                hft[nf] = (As[tfpair[tfp][0]*2+1][0]-As[tfpair[tfp][0]*2][0])/ \
                    (2.*As[tfpair[tfp][1]*2+1][0])
                    
            self.tf.append(hft)
            self.tf.append(vft)
        return self.tf
        
    def tf_knopoff_sh_adv(self):
        """
        
        """
        from numpy.linalg import solve
        # uniforming variable
        mode = self.mode
        ntf = self.ntf
        tfpair = self.tfpair
        nlayer = self.nlayer
        hl = self.hl
        vs = self.vs
        dn = self.dn
        qs = self.qs
        freq = self.freq
        # checking calculation validity
        #if mode != cd.mode[2] and mode!=cd.mode[3]:
        #    raise ValueError("Calculation mode is not supported! Use another method!")
        
        # angular frequency conversion
        angf = 2.*np.pi*freq
        
        # attenuation using complex velocities
        vs = vs*((2.*qs*1j)/(2.*qs*1j-1.))
        
        # angle of propagation within layers
        slayer = self.slayer
        iang = self.iang
        
        iD = np.zeros((nlayer,1))
        iCORE = np.zeros((nlayer,nlayer),dtype='complex128')
        
        iD[0] = np.sin(iang)
        iCORE[0,slayer] = 1.
        
        for nl in range(nlayer-1):
            iCORE[nl+1,nl] = 1./vs[nl]
            iCORE[nl+1,nl+1] = -1./vs[nl+1]
            
        iA = solve(iCORE,iD)
        
        iS = np.arcsin(iA)
        
        # Lame Parameter(s)
        mu = np.zeros((nlayer,1),dtype='complex128')
        
        for nl in range(nlayer):
            mu[nl]=dn[nl]*(vs[nl]**2)
        
        # horizontal and vertical slowness
        ns = np.zeros((nlayer,1),dtype='complex128')
        
        for nl in range(nlayer):
            ns[nl]=np.cos(iS[nl])/vs[nl]
            
        # building data vector
        A = np.zeros((nlayer*2,1))
        D = np.zeros((nlayer*2,1))
        D[-1] = 1.
        
        # Dispacement and transfer function initialization
        
        fnum = len(freq)
        
        # loop over frequencies and tfpair
        self.tf = []
        for tfp in range(ntf):
            hft = np.zeros((fnum),dtype='complex128')
            vft = np.zeros((fnum),dtype='complex128')
            for nf in range(fnum):
                # building core matrix
                CORE = np.zeros((nlayer*2,nlayer*2),dtype='complex128')
                
                # free surface constraints
                CORE[0,0] = 1.
                CORE[0,1] = -1.
                
                # Interfaces constraints
                for nl in range(nlayer-1):
                    row = (nl*2)+1
                    col = nl*2
                    
                    expDSA = np.exp(1j*angf[nf]*ns[nl]*hl[nl])
                    expUSA = np.exp(-1j*angf[nf]*ns[nl]*hl[nl])
                    
                    CORE[row,col+0] = expDSA[0]
                    CORE[row,col+1] = expUSA[0]
                    CORE[row,col+2] = -1.
                    CORE[row,col+3] = -1.
                    
                    CORE[row+1,col+0] =  mu[nl][0]*ns[nl][0]*expDSA[0]
                    CORE[row+1,col+1] = -mu[nl][0]*ns[nl][0]*expUSA[0]
                    CORE[row+1,col+2] = -mu[nl+1][0]*ns[nl+1][0]
                    CORE[row+1,col+3] =  mu[nl+1][0]*ns[nl+1][0]
                    
                # input constraints
                CORE[-1,-1]=1.
                
                # solving linear system
                try:
                    A=solve(CORE,D)
                except:
                    A[:] = np.nan
                    
                # transfer function
                hft[nf] = (A[tfpair[tfp][0]*2+1][0]+A[tfpair[tfp][0]*2][0])/ \
                    (2.*A[tfpair[tfp][1]*2+1][0])
                    
            self.tf.append(hft)
            self.tf.append(vft)
        return self.tf
        
    def tf_knopoff_psv_adv(self):
        """
        
        """
        from numpy.linalg import solve
        # uniforming variable
        mode = self.mode
        ntf = self.ntf
        tfpair = self.tfpair
        nlayer = self.nlayer
        hl = self.hl
        vs = self.vs
        dn = self.dn
        qs = self.qs
        freq = self.freq
        vp = self.vp
        qp = self.qp
        comp = self.comp
        iang = self.iang        
        
        # checking calculation validity
        #if mode != cd.mode[2] and mode!=cd.mode[3]:
        #    raise ValueError("Calculation mode is not supported! Use another method!")
        
        # angular frequency conversion
        angf = 2.*np.pi*freq
        
        # attenuation using complex velocities
        vp = vp*((2.*qp*1j)/(2.*qp*1j-1.))
        vs = vs*((2.*qs*1j)/(2.*qs*1j-1.))
        
        # angle of propagation within layers
        slayer = self.slayer
        iang = self.iang
        
        iD = np.zeros((nlayer*2),dtype='complex128')
        iCORE = np.zeros((nlayer*2,nlayer*2),dtype='complex128')
        
        if comp=='p':
            iD[0] = np.sin(iang)
            iD[1] = np.sin(iang)*vs[slayer]/vp[slayer]
        elif comp=='s':
            iD[0] = np.sin(iang)*vp[slayer]/vs[slayer]
            iD[1] = np.sin(iang)
        
        iCORE[0,2*slayer] = 1. # check the index
        iCORE[1,2*slayer+1] = 1.
        
        for nl in range(nlayer-1):
            
            row = nl*2
            col = nl*2
            
            iCORE[row+2,col+0] = 1./vp[nl]
            iCORE[row+2,col+2] = -1./vp[nl+1]
            
            iCORE[row+3,col+1] = 1/vs[nl]
            iCORE[row+3,col+3] = -1/vs[nl+1]
            
        iA = np.linalg.solve(iCORE,iD)
        
        iP = np.zeros((nlayer),dtype='complex128')
        iS = np.zeros((nlayer),dtype='complex128')

        for nl in range(nlayer):
            row = nl*2
            iP[nl] = np.arcsin(iA[row]) if np.real(iA[row]) <= 1.0 else np.arcsin(np.complex(np.real(iA[row]),-1.*np.abs(np.imag(iA[row]))))
            iS[nl] = np.arcsin(iA[row+1]) if np.real(iA[row+1]) <= 1.0 else np.arcsin(np.complex(np.real(iA[row+1]),-1.*np.abs(np.imag(iA[row+1]))))
        #print iA
        #print iP
        #print iS
        # Lame Parameter(s)
        mu = np.zeros((nlayer),dtype='complex128')
        la = np.zeros((nlayer),dtype='complex128')
        
        for nl in range(nlayer):
            mu[nl]=dn[nl]*(vs[nl]**2)
            la[nl]=dn[nl]*((vp[nl]**2)-(2.*(vs[nl]**2)))            

        # horizontal and vertical slowness
        rp = np.zeros((nlayer),dtype='complex128')
        rs = np.zeros((nlayer),dtype='complex128')
        np1 = np.zeros((nlayer),dtype='complex128')
        ns = np.zeros((nlayer),dtype='complex128')
        
        for nl in range(nlayer):
            rp[nl] = np.sin(iP[nl])/vp[nl]
            rs[nl] = np.sin(iS[nl])/vs[nl]
            np1[nl] = np.cos(iP[nl])/vp[nl]
            ns[nl] = np.cos(iS[nl])/vs[nl]

        # building data vector
        A = np.zeros((nlayer*4),dtype='complex128')
        D = np.zeros((nlayer*4),dtype='complex128')
        
        if comp=='p':
            D[-2] = vp[-1]
        elif comp=='s':
            D[-1] = vs[-1]
        
        # Dispacement and transfer function initialization
        
        fnum = len(freq)
        
        # loop over frequencies and tfpair
        self.tf = []
        for tfp in range(ntf):
            htft = np.zeros((fnum),dtype='complex128')
            vtft = np.zeros((fnum),dtype='complex128')
            for nf in range(fnum):
                # building core matrix
                CORE = np.zeros((nlayer*4,nlayer*4),dtype='complex128')
                
                # free surface constraints
                CORE[0,0] = la[0]*rp[0]**2+(la[0]+2.*mu[0])*np1[0]**2
                CORE[0,1] = 2.*mu[0]*rs[0]*ns[0]
                CORE[0,2] = la[0]*rp[0]**2+(la[0]+2.*mu[0])*np1[0]**2
                CORE[0,3] = -2.*mu[0]*rs[0]*ns[0]
                
                CORE[1,0] = -2.*mu[0]*rp[0]*np1[0]
                CORE[1,1] = -mu[0]*(rs[0]**2-ns[0]**2)
                CORE[1,2] = 2.*mu[0]*rp[0]*np1[0]
                CORE[1,3] = -mu[0]*(rs[0]**2-ns[0]**2)
                #print la[0],rp[0],mu[0],np1[0]
                # Interfaces constraints
                for nl in range(nlayer-1):
                    row = (nl*4)+2
                    col = nl*4
                    
                    expDPA = np.exp(1j*angf[nf]*np1[nl]*hl[nl])
                    expDSA = np.exp(1j*angf[nf]*ns[nl]*hl[nl])
                    expUPA = np.exp(-1j*angf[nf]*np1[nl]*hl[nl])
                    expUSA = np.exp(-1j*angf[nf]*ns[nl]*hl[nl])
                    
                    CORE[row+0,col+0] = rp[nl]*expDPA
                    CORE[row+0,col+1] = -ns[nl]*expDSA
                    CORE[row+0,col+2] = rp[nl]*expUPA
                    CORE[row+0,col+3] = ns[nl]*expUSA
                    CORE[row+0,col+4] = -rp[nl+1]
                    CORE[row+0,col+5] = ns[nl+1]
                    CORE[row+0,col+6] = -rp[nl+1]
                    CORE[row+0,col+7] = -ns[nl+1]
                    
                    CORE[row+1,col+0] = np1[nl]*expDPA
                    CORE[row+1,col+1] = rs[nl]*expDSA
                    CORE[row+1,col+2] = -np1[nl]*expUPA
                    CORE[row+1,col+3] = rs[nl]*expUSA
                    CORE[row+1,col+4] = -np1[nl+1]
                    CORE[row+1,col+5] = -rs[nl+1]
                    CORE[row+1,col+6] = np1[nl+1]
                    CORE[row+1,col+7] = -rs[nl+1]
                    
                    CORE[row+2,col+0] = (la[nl]*rp[nl]**2+(la[nl]+2.*mu[nl])*np1[nl]**2)*expDPA
                    CORE[row+2,col+1] = 2.*mu[nl]*rs[nl]*ns[nl]*expDSA
                    CORE[row+2,col+2] = (la[nl]*rp[nl]**2+(la[nl]+2.*mu[nl])*np1[nl]**2)*expUPA
                    CORE[row+2,col+3] = -2.*mu[nl]*rs[nl]*ns[nl]*expUSA
                    CORE[row+2,col+4] = -(la[nl+1]*rp[nl+1]**2+(la[nl+1]+2.*mu[nl+1])*np1[nl+1]**2)
                    CORE[row+2,col+5] = -2.*mu[nl+1]*rs[nl+1]*ns[nl+1]
                    CORE[row+2,col+6] = -(la[nl+1]*rp[nl+1]**2+(la[nl+1]+2.*mu[nl+1])*np1[nl+1]**2)
                    CORE[row+2,col+7] = 2.*mu[nl+1]*rs[nl+1]*ns[nl+1]
                    
                    CORE[row+3,col+0] = -2.*mu[nl]*rp[nl]*np1[nl]*expDPA
                    CORE[row+3,col+1] = -mu[nl]*(rs[nl]**2-ns[nl]**2)*expDSA
                    CORE[row+3,col+2] = 2.*mu[nl]*rp[nl]*np1[nl]*expUPA
                    CORE[row+3,col+3] = -mu[nl]*(rs[nl]**2-ns[nl]**2)*expUSA
                    CORE[row+3,col+4] = 2.*mu[nl+1]*rp[nl+1]*np1[nl+1]
                    CORE[row+3,col+5] = mu[nl+1]*(rs[nl+1]**2-ns[nl+1]**2)
                    CORE[row+3,col+6] = -2.*mu[nl+1]*rp[nl+1]*np1[nl+1]
                    CORE[row+3,col+7] = mu[nl+1]*(rs[nl+1]**2-ns[nl+1]**2)
                    
                # input constraints
                CORE[-2,-2]=1.  # <-- P
                CORE[-1,-1]=1.  # <-- S
                
                # solving linear system
                try:
                    #A=solve(CORE,D)
                    #CORE = np.matrix(CORE)
                    #D = np.matrix(D)
                    A = np.linalg.solve(CORE,D)
                    #A=np.dot(np.linalg.pinv(CORE,rcond=1e-40),D)
                    #import sympy as sp
                    #X = sp.Matrix(CORE)
                    #X=np.dot(np.linalg.inv(CORE),D)
                except ValueError:
                    print 'ValueError TFCalculater.py'
                    A[:] = np.nan
                    
                
                # transfer function
                htft[nf] = (rp[tfpair[tfp][0]]*(A[tfpair[tfp][0]+2]+A[tfpair[tfp][0]+0])+ \
                            ns[tfpair[tfp][0]]*(A[tfpair[tfp][0]+3]-A[tfpair[tfp][0]+1]))/ \
                            (2.*(rp[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+2])+ \
                            ns[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+3])))
                #print rp[tfpair[tfp][0]],(A[tfpair[tfp][0]+2]+A[tfpair[tfp][0]+0]), \
                #            ns[tfpair[tfp][0]],(A[tfpair[tfp][0]+3]-A[tfpair[tfp][0]+1])
                vtft[nf] = (np1[tfpair[tfp][0]]*(-A[tfpair[tfp][0]+2]+A[tfpair[tfp][0]+0])+ \
                            rs[tfpair[tfp][0]]*(A[tfpair[tfp][0]+3]+A[tfpair[tfp][0]+1]))/ \
                            (2.*(-np1[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+2])+ \
                            rs[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+3])))
                #vtft[nf] = (np1[tfpair[tfp][0]]*(-A[tfpair[tfp][0]+2]+A[tfpair[tfp][0]+0])+ \
                #            rs[tfpair[tfp][0]]*(A[tfpair[tfp][0]+3]+A[tfpair[tfp][0]+1]))/ \
                #            (2.*(-np1[-1]*(A[-2])+ \
                #            rs[-1]*(A[-1])))
                #if nf==0:
                #    print A[4*tfpair[tfp][1]+2],A[4*tfpair[tfp][1]+3]
                #    print rs
                #    print -np1[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+2]),rs[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+3])
            if iang==0.:
                if comp == 'p':
                    htft[:] = 0.
                elif comp == 's':
                    vtft[:] = 0.
            self.tf.append(htft)
            self.tf.append(vtft)
        return self.tf
    
    def tf_kennett_psv(self):
        """
        porting from geopsy based on kennet formalism
        Pierre-Yves Bard (LGIT, Grenoble, France)
        Cecile Cornou (LGIT, Grenoble, France), minor modifications
        Marc Wathelet (LGIT, Grenoble, France), port to subroutine
        Theodosius Marwan Irnaka (GEM, Pavia, Italy), port to python
        """
        
        def kenpsv(jcas,cwx,omega,nlayer,th,dn,cvp,cvs,verbose=False):
            """
            kennet formalism
            """
            
            # if th[0]=0, given the depth of the interfaces, else given the thickness
            hc = np.zeros(nlayer)            
            if th[0]==0.:
                for ic in range(nlayer):
                    hc[ic]=th[ic+1]-th[ic]
            else:
                hc = th

            cwx2 = cwx**2
            
            cka2 = (omega/cvp)**2
            ckb2= (omega/cvs)**2
            cnu = np.sqrt(cka2-cwx2)
            cgam = np.sqrt(ckb2-cwx2)
            for ic in range(nlayer):
                cnu[ic] =-cnu[ic] if np.imag(cnu[ic])>0. else cnu[ic]
                cgam[ic]=-cgam[ic] if np.imag(cgam[ic])>0. else cgam[ic]

            # calculation of reflection/transmission coefficients matrix and phase shift
            
            # coefficient for the convention of PSI (coef) and for TF
            aki = -1.
            
            ru = np.zeros((nlayer,2,2),dtype='complex128')
            tu = np.zeros((nlayer,2,2),dtype='complex128')
            rd = np.zeros((nlayer,2,2),dtype='complex128')
            td = np.zeros((nlayer,2,2),dtype='complex128')
            me1 = np.zeros((nlayer),dtype='complex128')
            me2 = np.zeros((nlayer),dtype='complex128')
            if jcas == 0:
                # coefficient for free surface
                cf1 = ckb2[0]-2.*cwx2
                cf2 = cf1**2
                cf3 = 4.*cnu[0]*cwx2*cgam[0]
                cdd = cf2+cf3
                
                ru[0,0,0] = (-cf2+cf3)/cdd
                ru[0,1,0] = 4.*cwx*cnu[0]*cf1/cdd*aki
                ru[0,1,1] = (cf2-cf3)/cdd*aki
                ru[0,0,1] = 4.*cwx*cgam[0]*cf1/cdd
                #tu[0,0,0] = 0.
                #tu[0,0,1] = 0.
                #tu[0,1,0] = 0.
                #tu[0,1,1] = 0.
            else:
                # coefficient for infinite space
                #ru[0,0,0] = 0.
                #ru[0,1,0] = 0.
                #ru[0,1,1] = 0.
                #ru[0,0,1] = 0.
                tu[0,0,0] = 1.
                #tu[0,0,1] = 0.
                #tu[0,1,0] = 0.
                tu[0,1,1] = 1.

            # coefficients at the interfaces between layers
            for ic in range(1,nlayer):
                cb1 = cwx2/ckb2[ic-1]
                cb2 = cwx2/ckb2[ic]
                ca1d= dn[ic-1]*(1.-2.*cb1)
                ca2d= dn[ic]*(1.-2.*cb2)
                ca  = ca2d-ca1d
                cb  = ca2d+2.*dn[ic-1]*cb1
                cc  = ca1d*2.*dn[ic]*cb2
                cd  = 2.*(dn[ic]/ckb2[ic]-dn[ic-1]/ckb2[ic-1])
                ce  = cb*cnu[ic-1]+cc*cnu[ic]
                cf  = cb*cgam[ic-1]+cc*cgam[ic]
                cg  = ca-cd*cnu[ic-1]*cgam[ic]
                ch  = ca-cd*cnu[ic]*cgam[ic-1]
                cdd = ce*cf+cg*ch*cwx2
                
                rd[ic,0,0] = (cf*(cb*cnu[ic-1]-cc*cnu[ic])- \
                    ch*cwx2*(ca+cd*cnu[ic-1]*cgam[ic]))/cdd
                rd[ic,0,1] =-2.*cwx*cgam[ic-1]* \
                    (ca*cb+cc*cd*cnu[ic]*cgam[ic])/cdd*aki
                rd[ic,1,1] =-(ce*(cb*cgam[ic-1]-cc*cgam[ic])- \
                    cg*cwx2*(ca+cd*cnu[ic]*cgam[ic-1]))/cdd*aki
                rd[ic,1,0] =-2.*cwx*cnu[ic-1]* \
                    (ca*cb+cc*cd*cnu[ic]*cgam[ic])/cdd
                    
                td[ic,0,0] = 2.*dn[ic-1]*cnu[ic-1]*cf/cdd
                td[ic,0,1] =-2.*dn[ic-1]*cgam[ic-1]*cg*cwx/cdd*aki
                td[ic,1,1] = 2.*dn[ic-1]*cgam[ic-1]*ce/cdd
                td[ic,1,0] = 2.*dn[ic-1]*cnu[ic-1]*ch*cwx/cdd*aki
                
                ru[ic,0,0] =-(cf*(cb*cnu[ic-1]-cc*cnu[ic])+ \
                    cg*cwx2*(ca+cd*cnu[ic]*cgam[ic-1]))/cdd
                ru[ic,0,1] = 2.*cwx*cgam[ic]* \
                    (ca*cc+cb*cd*cd*cnu[ic-1]*cgam[ic-1])/cdd
                ru[ic,1,1] = (ce*(cb*cgam[ic-1]-cc*cgam[ic])+ \
                    ch*cwx2*(ca+cd*cnu[ic-1]*cgam[ic]))/cdd*aki
                ru[ic,1,0] = 2.*cwx*cnu[ic]* \
                    (ca*cc+cb*cd*cnu[ic-1]*cgam[ic-1])/cdd*aki
                tu[ic,0,0] = 2.*dn[ic]*cnu[ic]*cf/cdd
                tu[ic,0,1] = 2.*dn[ic]*cgam[ic]*ch*cwx/cdd
                tu[ic,1,1] = 2.*dn[ic]*cgam[ic]*ce/cdd
                tu[ic,1,0] =-2.*dn[ic]*cnu[ic]*cg*cwx/cdd
                
                me1[ic-1] = np.exp(-1j*cnu[ic-1]*hc[ic-1])
                me2[ic-1] = np.exp(-1j*cgam[ic-1]*hc[ic-1])
                
            # calculation of reflectivity matrix : mt(), mb(), nt() nb()
            
            # calculation for the layers above the source
            nc = deepcopy(nlayer)
            nt   = np.zeros((nlayer,2,2),dtype='complex128')
            mt   = np.zeros((nlayer,2,2),dtype='complex128')
            fdo  = np.zeros((nlayer,2,2),dtype='complex128')
            fup  = np.zeros((nlayer,2,2),dtype='complex128')
            nb = np.zeros((2,2),dtype='complex128')
            mb = np.zeros((2,2),dtype='complex128')
            
            nt[0,0,0] = ru[0,0,0]
            nt[0,0,1] = ru[0,0,1]
            nt[0,1,0] = ru[0,1,0]
            nt[0,1,1] = ru[0,1,1]

            for ic in range(nc-1):
                nb[0,0] = me1[ic]*me1[ic]*nt[ic,0,0]
                nb[0,1] = me1[ic]*me2[ic]*nt[ic,0,1]
                nb[1,0] = me2[ic]*me1[ic]*nt[ic,1,0]
                nb[1,1] = me2[ic]*me2[ic]*nt[ic,1,1]
                
                ca1 = 1.-(rd[ic+1,0,0]*nb[0,0]+rd[ic+1,0,1]*nb[1,0])
                ca2 =   -(rd[ic+1,0,0]*nb[0,1]+rd[ic+1,0,1]*nb[1,1])
                ca3 =   -(rd[ic+1,1,0]*nb[0,0]+rd[ic+1,1,1]*nb[1,0])
                ca4 = 1.-(rd[ic+1,1,0]*nb[0,1]+rd[ic+1,1,1]*nb[1,1])
                cadet = ca1*ca4-ca2*ca3
                
                cb1 = td[ic+1,0,0]*nb[0,0]+td[ic+1,0,1]*nb[1,0]
                cb2 = td[ic+1,0,0]*nb[0,1]+td[ic+1,0,1]*nb[1,1]
                cb3 = td[ic+1,1,0]*nb[0,0]+td[ic+1,1,1]*nb[1,0]
                cb4 = td[ic+1,1,0]*nb[0,1]+td[ic+1,1,1]*nb[1,1]
                
                cc1 = (ca4*tu[ic+1,0,0]-ca2*tu[ic+1,1,0])/cadet
                cc2 = (ca4*tu[ic+1,0,1]-ca2*tu[ic+1,1,1])/cadet
                cc3 = (-ca3*tu[ic+1,0,0]+ca1*tu[ic+1,1,0])/cadet
                cc4 = (-ca3*tu[ic+1,0,1]+ca1*tu[ic+1,1,1])/cadet
                
                nt[ic+1,0,0] = ru[ic+1,0,0]+cb1*cc1+cb2*cc3
                nt[ic+1,0,1] = ru[ic+1,0,1]+cb1*cc2+cb2*cc4
                nt[ic+1,1,0] = ru[ic+1,1,0]+cb3*cc1+cb4*cc3
                nt[ic+1,1,1] = ru[ic+1,1,1]+cb3*cc2+cb4*cc4
                
                fup[ic,0,0] = cc1*me1[ic]
                fup[ic,0,1] = cc2*me1[ic]
                fup[ic,1,0] = cc3*me2[ic]
                fup[ic,1,1] = cc4*me2[ic]
            
            # calculation for the laters below the source
            
            #mt[nc-1,0,0] = 0.
            #mt[nc-1,0,1] = 0.
            #mt[nc-1,1,0] = 0.
            #mt[nc-1,1,1] = 0.

            for ic in range(nc-2,-1,-1):
                ca1 = 1.-(ru[ic+1,0,0]*mt[ic+1,0,0]+ru[ic+1,0,1]*mt[ic+1,1,0])
                ca2 =   -(ru[ic+1,0,0]*mt[ic+1,0,1]+ru[ic+1,0,1]*mt[ic+1,1,1])
                ca3 =   -(ru[ic+1,1,0]*mt[ic+1,0,0]+ru[ic+1,1,1]*mt[ic+1,1,0])
                ca4 = 1.-(ru[ic+1,1,0]*mt[ic+1,0,1]+ru[ic+1,1,1]*mt[ic+1,1,1])
                cadet = ca1*ca4-ca2*ca3
                
                cb1 = tu[ic+1,0,0]*mt[ic+1,0,0]+tu[ic+1,0,1]*mt[ic+1,1,0]
                cb2 = tu[ic+1,0,0]*mt[ic+1,0,1]+tu[ic+1,0,1]*mt[ic+1,1,1]
                cb3 = tu[ic+1,1,0]*mt[ic+1,0,0]+tu[ic+1,1,1]*mt[ic+1,1,0]
                cb4 = tu[ic+1,1,0]*mt[ic+1,0,1]+tu[ic+1,1,1]*mt[ic+1,1,1]
                
                cc1 = ( ca4*td[ic+1,0,0]-ca2*td[ic+1,1,0])/cadet
                cc2 = ( ca4*td[ic+1,0,1]-ca2*td[ic+1,1,1])/cadet
                cc3 = (-ca3*td[ic+1,0,0]+ca1*td[ic+1,1,0])/cadet
                cc4 = (-ca3*td[ic+1,0,1]+ca1*td[ic+1,1,1])/cadet
                
                mb[0,0] = rd[ic+1,0,0]+cb1*cc1+cb2*cc3
                mb[0,1] = rd[ic+1,0,1]+cb1*cc2+cb2*cc4
                mb[1,0] = rd[ic+1,1,0]+cb3*cc1+cb4*cc3
                mb[1,1] = rd[ic+1,1,1]+cb3*cc2+cb4*cc4
                
                mt[ic,0,0] = me1[ic]*me1[ic]*mb[0,0]
                mt[ic,0,1] = me1[ic]*me2[ic]*mb[0,1]
                mt[ic,1,0] = me2[ic]*me1[ic]*mb[1,0]
                mt[ic,1,1] = me2[ic]*me2[ic]*mb[1,1]

                fdo[ic+1,0,0] = cc1*me1[ic]
                fdo[ic+1,0,1] = cc2*me2[ic]
                fdo[ic+1,1,0] = cc3*me1[ic]
                fdo[ic+1,1,1] = cc4*me2[ic]
            
            # calculation upgoing and downgoing P and S wave for each layer"
                # - upgoing P wave (jcas = 0) in the layer ln
                # - upgoing S wave (jcas = 0) in the layer ln
                # - downgoing P wave (jcas = 1) in the layer ln
                # - downgoing S wave (jcas = 1) in the layer ln
            
            # reflect4(jcas)
            
            ftup = np.zeros((nlayer,2,2),dtype='complex128')
            pu = np.zeros((nlayer,2,2),dtype='complex128')
            pd = np.zeros((nlayer,2,2),dtype='complex128')
            ftdo = np.zeros((nlayer,2,2),dtype='complex128')
            cfwave = np.zeros((nlayer*4,3),dtype='complex128')
            if jcas==0:
                # case jcas=0 free surface reflection
                ftup[nc-1,0,0] = 1.
                #ftup[nc-1,0,1] = 0.
                #ftup[nc-1,1,0] = 0.
                ftup[nc-1,1,1] = 1.
                
                for ic in range(nc-2,-1,-1):
                    ftup[ic,0,0] = fup[ic,0,0]*ftup[ic+1,0,0]+fup[ic,0,1]*ftup[ic+1,1,0]
                    ftup[ic,0,1] = fup[ic,0,0]*ftup[ic+1,0,1]+fup[ic,0,1]*ftup[ic+1,1,1]
                    ftup[ic,1,0] = fup[ic,1,0]*ftup[ic+1,0,0]+fup[ic,1,1]*ftup[ic+1,1,0]
                    ftup[ic,1,1] = fup[ic,1,0]*ftup[ic+1,0,1]+fup[ic,1,1]*ftup[ic+1,1,1]
                    
                # ???? Vectors potential amount (pu) and down (pd) in each layer receiver for 6 elementary sources
                # Receivers above the source
                    
                for ic in range(nc):
                    ic1 = 4*ic
                    ic2 = ic1+1
                    ic3 = ic2+1
                    ic4 = ic3+1
                    pu[ic,0,0] = ftup[ic,0,0]
                    pu[ic,1,0] = ftup[ic,1,0]
                    pd[ic,0,0] = nt[ic,0,0]*pu[ic,0,0]+nt[ic,0,1]*pu[ic,1,0]
                    pd[ic,1,0] = nt[ic,1,0]*pu[ic,0,0]+nt[ic,1,1]*pu[ic,1,0]
                    cfwave[ic1,0] = pu[ic,0,0]
                    cfwave[ic2,0] = pd[ic,0,0]
                    cfwave[ic3,0] = pu[ic,1,0]
                    cfwave[ic4,0] = pd[ic,1,0]
                    
                    pu[ic,0,1] = ftup[ic,0,1]
                    pu[ic,1,1] = ftup[ic,1,1]
                    pd[ic,0,1] = nt[ic,0,0]*pu[ic,0,1]+nt[ic,0,1]*pu[ic,1,1]
                    pd[ic,1,1] = nt[ic,1,0]*pu[ic,0,1]+nt[ic,1,1]*pu[ic,1,1]
                    cfwave[ic1,1] = pu[ic,0,1]
                    cfwave[ic2,1] = pd[ic,0,1]
                    cfwave[ic3,1] = pu[ic,1,1]
                    cfwave[ic4,1] = pd[ic,1,1]
            else:
                # case jcas=1 : no upwave incident
                ftdo[0,0,0]=1.
                #ftdo[0,0,1]=0.
                #ftdo[0,1,0]=0.
                ftdo[0,1,1]=1.
                
                for ic in range(1,nc):
                    ftdo[ic,0,0] = fdo[ic,0,0]*ftdo[ic-1,0,0]+fdo[ic,0,1]*ftdo[ic-1,1,0]
                    ftdo[ic,0,1] = fdo[ic,0,0]*ftdo[ic-1,0,1]+fdo[ic,0,1]*ftdo[ic-1,1,1]
                    ftdo[ic,1,0] = fdo[ic,1,0]*ftdo[ic-1,0,0]+fdo[ic,1,1]*ftdo[ic-1,1,0]
                    ftdo[ic,1,1] = fdo[ic,1,0]*ftdo[ic-1,0,1]+fdo[ic,1,1]*ftdo[ic-1,1,1]
                    
                for ic in range(nc):
                    ic1 = 4*ic
                    ic2 = ic1+1
                    ic3 = ic2+1
                    ic4 = ic3+1
                    
                    pd[ic,0,0] = ftdo[ic,0,0]
                    pd[ic,1,0] = ftdo[ic,1,0]
                    pu[ic,0,0] = mt[ic,0,0]*pd[ic,0,0]+mt[ic,0,1]*pd[ic,1,0]
                    pu[ic,1,0] = mt[ic,1,0]*pd[ic,0,0]+mt[ic,1,1]*pd[ic,1,0]
                    cfwave[ic1,0] = pu[ic,0,0]
                    cfwave[ic2,0] = pd[ic,0,0]
                    cfwave[ic3,0] = pu[ic,1,0]
                    cfwave[ic4,0] = pd[ic,1,0]
                    
                    pd[ic,0,1] = ftdo[ic,0,1]
                    pd[ic,1,1] = ftdo[ic,1,1]
                    pu[ic,0,1] = mt[ic,0,0]*pd[ic,0,1]+mt[ic,0,1]*pd[ic,1,1]
                    pu[ic,1,1] = mt[ic,1,0]*pd[ic,0,1]+mt[ic,1,1]*pd[ic,1,1]
                    cfwave[ic1,1] = pu[ic,0,1]
                    cfwave[ic2,1] = pu[ic,0,1]
                    cfwave[ic3,1] = pu[ic,1,1]
                    cfwave[ic4,1] = pu[ic,1,1]
                    
            return cfwave
        
        # uniforming variable
        mode = self.mode
        modeID = self.modeID
        ntf = self.ntf
        tfpair = self.tfpair
        nlayer = self.nlayer
        hl = self.hl
        vs = self.vs
        dn = self.dn
        qs = self.qs
        freq = self.freq
        vp = self.vp
        qp = self.qp
        comp = self.comp
        
        #freq0 = 10.         # reference frequency for Q <-- WHY?
        q = 1.e+20          # ????
        iang = self.iang    # incidence angle (radians)
        if comp=='s':
            mode = 2
        else:
            mode = 1
        jcas = 0            # jcas <-- 0 means with free surface, 1 means infinite medium
        nr = len(tfpair)    # number of receiver given by the pair of input and output motion
        zr = np.array([tfpair[i][0]*hl for i in range(len(tfpair))])
                            # depth of receiver given as list (BEWARE OF INCOMPATIBILITY!!)
        
        # defining depth of the top layer (can be modified later)
        th = np.zeros(nlayer+1)
        for i in range(1,nlayer+1):
            th[i] = th[i-1]+hl[i-1]
            
        # portion of signal definition and calculationg of Fourier Transform
            # omitted
        
        nf = len(freq)
        aw = -np.pi/q           # I don't understand! it's basically 0!
        
        # index of receiver
        izr = np.array([tfpair[i][1] for i in range(len(tfpair))])
        
        # iterating over frequencies
        u = np.zeros((nf,nr),dtype='complex128')
        v = np.zeros((nf,nr),dtype='complex128')
        w = np.zeros((nf,nr),dtype='complex128')
       
        self.tf = []
        for i,fr in enumerate(freq):
            fr = 0.05*fr if fr==0. else fr  # correction for zero frequency
            rw = fr*np.pi*2.
            omega = np.complex(rw,aw)
            #df = fr - freq0
            #df = 0. if df<0. else df
            
            cvp = vp*((2.*qp*1j)/(2.*qp*1j-1.))
            cvs = vs*((2.*qs*1j)/(2.*qs*1j-1.)) 
            wb = omega/cvs      # omega for vs
            wb2 = wb**2
            wa = omega/cvp      # oemga for vp
            wa2 = wa**2
            
            #n = 1.
            if comp=='p':
                wx0 = wa[nlayer-1]*np.sin(iang)
                #mx0 = n*fr*np.sin(iang)/vp[nlayer-1]
                c1 = -1j/wa[nlayer-1]
            else:
                wx0 = wb[nlayer-1]*np.sin(iang)
                #mx0 = n*fr*np.sin(iang)/vs[nlayer-1]
                c1 = 1j/wb[nlayer-1]
            wx02 = wx0**2
            if i==0:
                f = kenpsv(jcas,wx0,omega,nlayer,th,dn,cvp,cvs,verbose=True)
            else:
                f = kenpsv(jcas,wx0,omega,nlayer,th,dn,cvp,cvs)
            # if mode is 1 or 2 :
            #   1 = upgoing p wave
            #   2 = downgoing p wave
            #   3 = upgoing sv wave
            #   4 = downgoing sv wave
            
            # if mode = 3 :
            #   1 = upgoing sh wave
            #   2 = downgoing sh wave
            wza = np.zeros_like(cvp)
            wzb = np.zeros_like(cvp)
            wza = np.sqrt(wa2-wx02)
            wzb = np.sqrt(wb2-wx02)
            for li in range(nlayer):
                wza[li] = -wza[li] if np.imag(wza[li])>0 else wza[li]
                wzb[li] = -wzb[li] if np.imag(wzb[li])>0 else wzb[li]
            if i==0:
                print np.shape(f)
            for ir in range(nr):
                li = izr[ir]
                li1=4*li
                li2=li1+1
                li3=li2+1
                li4=li3+1
                z=zr[0][ir]-th[li]
                phaspu=np.exp(1j*wza[li]*z)
                phassu=np.exp(1j*wzb[li]*z)
                phaspd=1./phaspu
                phassd=1./phassu
                u1 = -1j*wx0
                u4 = 1j*wzb[li]
                w1 = 1j*wza[li]
                modes = mode-1
                f1 = f[li1,modes]*phaspu
                f2 = f[li2,modes]*phaspd
                f3 = f[li3,modes]*phassu
                f4 = f[li4,modes]*phassd
                uh0 = u1*(f1+f2) + u4*(f4-f3)
                uv0 = w1*(f1-f2) + u1*(f3+f4)
                u[i,ir] = uh0*c1
                w[i,ir] = uv0*c1
                
                if i==0:
                    print np.abs(u[i,ir])
        if mode==1 or mode==2:
            htf = u/2.
            vtf = w/2.
            
        self.tf.append(htf[:,0])
        self.tf.append(vtf[:,0])
        return self.tf
        
    def RTcoefficientsSH(self,rho,beta,iangS,calctype='rssd'):        
        """
        Function to calculate Reflection and Transmission coefficient on SH case
        4 possible calculation types:
            - rssd
                     \    /
                      \  /
                       \/
                ================
                
            - rssu
                ================
                       /\
                      /  \
                     /    \
                     
            - tssd
                     \
                      \
                ================
                         \
                          \
                          
            - tssu
                          /
                         /
                ================
                      /
                     /
        """
        if calctype.lower() == 'tssd' or calctype.lower() == 'rssd':
            iang = [iangS, np.arcsin(np.sin(iangS)*beta[1]/beta[0])]
        elif calctype.lower() == 'tssu' or calctype.lower() == 'rssu':
            iang = [np.arcsin(np.sin(iangS)*beta[0]/beta[1]), iangS]
        Delta = rho[0]*beta[0]*np.cos(iang[0])+rho[1]*beta[1]*np.cos(iang[1])
        if calctype.lower() == 'rssd':
            return (rho[0]*beta[0]*np.cos(iang[0])-rho[1]*beta[1]*np.cos(iang[1]))/Delta
        elif calctype.lower() == 'rssu':
            return -(rho[0]*beta[0]*np.cos(iang[0])-rho[1]*beta[1]*np.cos(iang[1]))/Delta
        elif calctype.lower() == 'tssd':
            return (2.*rho[0]*beta[0]*np.cos(iang[0]))/Delta
        elif calctype.lower() == 'tssu':
            return (2.*rho[1]*beta[1]*np.cos(iang[1]))/Delta
        else:
            raise IOError('undefined calculation type')
            
    def RTcoefficientsPSV(self,rho,alpha,beta,iang,calctype='rppd',iangtype='p'):
        """
        Same analogy compared to RTcoefficientsSH
        
        iangtype = incident angle type --> 'p' or 's'
        """
        rho = np.asarray(rho)
        alpha = np.asarray(alpha)
        beta = np.asarray(beta)
        
        if iangtype.lower()=='p':
            if calctype[3].lower()=='d':
                iangP = [iang,np.arcsin(np.sin(iang)*alpha[1]/alpha[0])]
            elif calctype[3].lower()=='u':
                iangP = [np.arcsin(np.sin(iang)*alpha[0]/alpha[1]),iang]
            iangS = [np.arcsin(np.sin(iang)*beta[0]/alpha[0]), \
                     np.arcsin(np.sin(iang)*beta[1]/alpha[0])]
        elif iangtype.lower()=='s':
            if calctype[3].lower()=='d':
                iangS = [iang,np.arcsin(np.sin(iang)*beta[1]/beta[0])]
            elif calctype[3].lower()=='u':
                iangS = [np.arcsin(np.sin(iang)*beta[0]/beta[1]),iang]
            iangP = [np.arcsin(np.sin(iang)*alpha[0]/beta[0]), \
                     np.arcsin(np.sin(iang)*alpha[1]/beta[0])]
        else:
            raise IOError('incidence wave type is not supported!')
                     
        beta2 = beta**2
        p = np.sin(iangP[0])/alpha[0]      # ray parameter
        p2 = p**2
        
        a = rho[1]*(1.-2.*beta2[1]*p2)-rho[0]*(1.-2.*beta2[0]*p2)
        b = rho[1]*(1.-2.*beta2[1]*p2)+2.*rho[0]*beta2[0]*p2
        c = rho[0]*(1.-2.*beta2[0]*p2)+2.*rho[1]*beta2[1]*p2
        d = 2.*(rho[1]*beta2[1]-rho[0]*beta2[0])
         
        E = b*(np.cos(iangP[0])/alpha[0])+c*(np.cos(iangP[1])/alpha[1])
        F = b*(np.cos(iangS[0])/beta[0])+c*(np.cos(iangS[1])/beta[1])
        G = a - d*np.cos(iangP[0])*np.cos(iangS[1])/(alpha[0]*beta[1])
        H = a - d*np.cos(iangP[1])*np.cos(iangS[0])/(alpha[1]*beta[0])
        D = E*F+G*H*p2
        
        if calctype.lower() == 'rppd':
            return (((b*(np.cos(iangP[0])/alpha[0])-c*(np.cos(iangP[1])/alpha[1]))*F)- \
                    ((a + d*np.cos(iangP[0])*np.cos(iangS[1])/(alpha[0]*beta[1]))*H*p2))/D
        elif calctype.lower() == 'rpsd':
            return -2.*(np.cos(iangP[0])/alpha[0])*(a*b+c*d*np.cos(iangP[1])* \
                    np.cos(iangS[1])/(alpha[1]*beta[1]))*p*alpha[0]/(beta[0]*D)
        elif calctype.lower() == 'tppd':
            return 2.*rho[0]*(np.cos(iangP[0])/alpha[0])*F*alpha[0]/(alpha[1]*D)
        elif calctype.lower() == 'tpsd':
            return 2.*rho[0]*(np.cos(iangP[0])/alpha[0])*H*p*alpha[0]/(beta[1]*D)
        elif calctype.lower() == 'rspd':
            return -2.*(np.cos(iangS[0])/beta[0])*(a*b+c*d*np.cos(iangP[1])* \
                    np.cos(iangS[1])/(alpha[1]*beta[1]))*p*beta[0]/(alpha[0]*D)
        elif calctype.lower() == 'rssd':
            return -(((b*(np.cos(iangS[0])/beta[0])-c*(np.cos(iangS[1])/beta[1]))*E)- \
                    ((a + d*np.cos(iangP[1])*np.cos(iangS[0])/(alpha[1]*beta[0]))*G*p2))/D
        elif calctype.lower() == 'tspd':
            return -2.*rho[0]*(np.cos(iangS[0])/beta[0])*G*p*beta[0]/(alpha[1]*D)
        elif calctype.lower() == 'tssd':
            return 2.*rho[0]*(np.cos(iangS[0])/beta[0])*E*beta[0]/(beta[1]*D)
        elif calctype.lower() == 'tppu':
            return 2.*rho[1]*(np.cos(iangP[1])/alpha[1])*F*alpha[1]/(alpha[0]*D)
        elif calctype.lower() == 'tpsu':
            #return 2.*rho[1]*(np.cos(iangP[1])/alpha[1])*H*p*alpha[1]/(beta[0]*D)  # check it.. I think the formula on aki richard was wrong!
            return -2.*rho[1]*(np.cos(iangP[1])/alpha[1])*G*p*alpha[1]/(alpha[1]*D)
        elif calctype.lower() == 'rppu':
            return -(((b*(np.cos(iangP[0])/alpha[0])-c*(np.cos(iangP[1])/alpha[1]))*F)+ \
                    ((a + d*np.cos(iangP[1])*np.cos(iangS[0])/(alpha[1]*beta[0]))*G*p2))/D
        elif calctype.lower() == 'rpsu':
            return 2.*(np.cos(iangP[1])/alpha[1])*(a*c+b*d*np.cos(iangP[0])* \
                    np.cos(iangS[0])/(alpha[0]*beta[0]))*p*alpha[1]/(beta[1]*D)
        elif calctype.lower() == 'tspu':
            return 2.*rho[1]*(np.cos(iangS[1])/beta[1])*H*p*beta[1]/(alpha[0]*D)
        elif calctype.lower() == 'tssu':
            return  2.*rho[1]*(np.cos(iangS[1])/beta[1])*E*beta[1]/(beta[0]*D)
        elif calctype.lower() == 'rspu':
            return 2.*(np.cos(iangS[1])/beta[1])*(a*c+b*d*np.cos(iangP[0])* \
                    np.cos(iangS[0])/(alpha[0]*beta[0]))*p*beta[1]/(alpha[1]*D)
        elif calctype.lower() == 'rssu':
            return (((b*(np.cos(iangS[0])/beta[0])-c*(np.cos(iangS[1])/beta[1]))*E)+ \
                    ((a + d*np.cos(iangP[0])*np.cos(iangS[1])/(alpha[0]*beta[1]))*H*p2))/D
        else:
            raise IOError('Unsupported type of calculation!')
                         
    def tf_kennet_sh(self):
        """
        porting from geopsy based on kennet formalism
        Pierre-Yves Bard (LGIT, Grenoble, France)
        Cecile Cornou (LGIT, Grenoble, France), minor modifications
        Marc Wathelet (LGIT, Grenoble, France), port to subroutine
        Theodosius Marwan Irnaka (GEM, Pavia, Italy), port to python
        """
        
        def kensh(jcas,cwx,omega,nlayer,th,dn,cvs,verbose=False):
            """
            kennet formalism
            """
            
            # if th[0]=0, given the depth of the interfaces, else given the thickness
            hc = np.zeros(nlayer)            
            if th[0]==0.:
                for ic in range(nlayer):
                    hc[ic]=th[ic+1]-th[ic]
            else:
                hc = th
            
            cwx2 = cwx**2
            
            ckb2= (omega/cvs)**2
            cgam = np.sqrt(ckb2-cwx2)
            for ic in range(nlayer):
                cgam[ic]=-cgam[ic] if np.imag(cgam[ic])>0. else cgam[ic]

            # calculation of reflection/transmission coefficients matrix and phase shift
            
            # coefficient for the convention of PSI (coef) and for TF
            rush = np.zeros((nlayer),dtype='complex128')
            tush = np.zeros((nlayer),dtype='complex128')
            rdsh = np.zeros((nlayer),dtype='complex128')
            tdsh = np.zeros((nlayer),dtype='complex128')
            me = np.zeros((nlayer),dtype='complex128')
            if jcas == 0:
                # coefficient for free surface (perfect reflection)
                rush[0] = 1.
                #tush[0] = 0.
            else:
                # coefficient for infinite space (perfect transmission)
                #rush[0] = 0.
                tush[0] = 1.

            # coefficients at the interfaces between layers
            for ic in range(1,nlayer):
                me[ic-1] = np.exp(-1j*cgam[ic-1]*hc[ic-1])
                
                cs1 = dn[ic-1]/ckb2[ic-1]*cgam[ic-1]
                cs2 = dn[ic]/ckb2[ic]*cgam[ic]
                cdelt= cs1+cs2
                
                rush[ic] = (cs2-cs1)/cdelt
                rdsh[ic] =-rush[ic]
                tush[ic] = 2.*cs2/cdelt
                tdsh[ic] = 2.*cs1/cdelt
            
            # calculation of reflectivity matrix : mt(), mb(), nt() nb()
            
            # calculation for the layers above the source
            nc = deepcopy(nlayer)
            ntsh = np.zeros((nlayer),dtype='complex128')
            mtsh = np.zeros((nlayer),dtype='complex128')
            fupsh= np.zeros((nlayer),dtype='complex128')
            fdosh= np.zeros((nlayer),dtype='complex128')  
            
            ntsh[0] = rush[0]            
            for ic in range(nc-1):
                nbsh = me[ic]*me[ic]*ntsh[ic]
                cash = 1./(1.-rdsh[ic+1]*nbsh)
                cbsh = tdsh[ic+1]*nbsh
                ccsh = cash*tush[ic+1]
                ntsh[ic+1] = rush[ic+1]+cbsh*ccsh
                fupsh[ic] = ccsh*me[ic]
            
            # calculation for the laters below the source
            #mtsh[nc-1] = 0.
            for ic in range(nc-2,-1,-1):
                cash = 1./(1.-rush[ic+1]*mtsh[ic+1])
                cbsh = tush[ic+1]*mtsh[ic+1]
                ccsh = cash*tdsh[ic+1]
                mbsh = rdsh[ic+1]*cbsh*ccsh
                mtsh[ic] = me[ic]*me[ic]*mbsh
                fdosh[ic+1] = ccsh*me[ic]
            
            # reflect4(jcas)
            
            push = np.zeros((nlayer),dtype='complex128')
            pdsh = np.zeros((nlayer),dtype='complex128')
            ftupsh = np.zeros((nlayer),dtype='complex128')
            ftdosh = np.zeros((nlayer),dtype='complex128')
            cfwave = np.zeros((nlayer*2),dtype='complex128')
            if jcas==0:
                # case jcas=0 free surface reflection
                ftupsh[nc-1] = 1.
                
                for ic in range(nc-2,-1,-1):
                    ftupsh[ic] = fupsh[ic]*ftupsh[ic+1]
                    
                # ???? Vectors potential amount (pu) and down (pd) in each layer receiver for 6 elementary sources
                # Receivers above the source
                    
                for ic in range(nc):                    
                    icsh1 = 2*ic
                    icsh2 = icsh1+1
                    push[ic] = ftupsh[ic]
                    pdsh[ic] = ntsh[ic]*push[ic]
                    cfwave[icsh1] = push[ic]
                    cfwave[icsh2] = pdsh[ic]
            else:
                # case jcas=1 : no upwave incident
                ftdosh[0]=1.
                
                for ic in range(1,nc):
                    ftdosh[ic] = fdosh[ic]*ftdosh[ic-1]
                    
                for ic in range(nc):
                    icsh1 = 2*ic
                    icsh2 = icsh1+1
                    pdsh[ic] = ftdosh[ic]
                    push[ic] = mtsh[ic]*pdsh[ic]
                    cfwave[icsh1] = push[ic]
                    cfwave[icsh2] = pdsh[ic]
            return cfwave
        
        # uniforming variable
        mode = self.mode
        modeID = self.modeID
        ntf = self.ntf
        tfpair = self.tfpair
        nlayer = self.nlayer
        hl = self.hl
        vs = self.vs
        dn = self.dn
        qs = self.qs
        freq = self.freq
        
        q = 1.e+20          # ????
        iang = self.iang    # incidence angle (radians)
        jcas = 0            # jcas <-- 0 means with free surface, 1 means infinite medium
        nr = len(tfpair)    # number of receiver given by the pair of input and output motion
        zr = np.array([tfpair[i][0]*hl for i in range(len(tfpair))])
                            # depth of receiver given as list (BEWARE OF INCOMPATIBILITY!!)
        
        # defining depth of the top layer (can be modified later)
        th = np.zeros(nlayer+1)
        for i in range(1,nlayer+1):
            th[i] = th[i-1]+hl[i-1]
        
        nf = len(freq)
        aw = -np.pi/q           # I don't understand! it's basically 0!

        # index of receiver
        izr = np.array([tfpair[i][0] for i in range(len(tfpair))])
        
        # iterating over frequencies
        v = np.zeros((nf,nr),dtype='complex128')
        self.tf = []
        for i,fr in enumerate(freq):
            fr = 0.05*fr if fr==0. else fr  # correction for zero frequency
            rw = fr*np.pi*2.
            omega = np.complex(rw,aw)
            
            # calculating complex velocity
            cvs = vs*((2.*qs*1j)/(2.*qs*1j-1.)) 
            wb = omega/cvs      # omega for vs
            wb2 = wb**2
            
            wx0 = wb[nlayer-1]*np.sin(iang)
            wx02 = wx0**2
            if i==0:
                f = kensh(jcas,wx0,omega,nlayer,th,dn,cvs,verbose=True)
            else:
                f = kensh(jcas,wx0,omega,nlayer,th,dn,cvs)
            
            # if mode = 3 :
            #   1 = upgoing sh wave
            #   2 = downgoing sh wave
            wzb = np.zeros_like(cvs)
            wzb = np.sqrt(wb2-wx02)
            for li in range(nlayer):
                wzb[li] = -wzb[li] if np.imag(wzb[li])>0 else wzb[li]
                
            for ir in range(nr):
                li = izr[ir]
                lish1=2*li
                lish2=lish1+1
                z=zr[0][ir]-th[li]
                phassu=np.exp(1j*wzb[li]*z)
                phassd=1./phassu
                v[i,ir]= f[lish1]*phassu+f[lish2]*phassd
                
                if i==0:
                    print np.abs(v[i,ir])
            htft = v/2.
            vtft = np.zeros_like(htft)
        self.tf.append(htft[:,0])
        self.tf.append(vtft[:,0])
        return self.tf
        

# debugging
"""
import IOfile
fname2 = 'sampleinput_linear_elastic_6layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
theclass2 = TFCalculator(data2)
theclass2.tf_kennet_sh()
"""
"""
import IOfile
fname2 = 'sampleinput_psv_p_linear_elastic_1layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
theclass2 = TFCalculator(data2)
htf,vtf = theclass2.tf_kennett()
print np.shape(htf)
import pylab as plt
plt.plot(np.abs(htf))
plt.grid(True,which='both')
plt.xscale('log')
plt.yscale('log')
plt.axis('Tight')
print htf
"""