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

class TFCalculator:
    def __init__(self,data):
        self.mode = data['mode']                  # calculation mode
        self.ntf = data['ntf']                   # number of transfer function pairs
        self.tfpair = np.array(data['tfPair'])         # list of transfer function pairs
        self.nlayer = data['nlayer']                # number of soil layer
        self.hl = np.array(data['hl'])          # thickness of each layer
        self.vs = np.array(data['vs'])          # shear wave velocity
        self.dn = np.array(data['dn'])          # density
        self.qs = np.array(data['qs'])          # quality factor
        self.slayer = data['sourceloc']
        self.iang = data['iang']
        self.freq = np.linspace(1.,50.,1024)
        
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
        
        Gc = np.zeros((len(self.hl),len(self.freq)),dtype='complex64')
        A = np.zeros_like(Gc)
        B = np.zeros_like(Gc)
        eta = np.zeros_like(self.qs)
        amp = np.zeros((len(self.freq)),dtype='complex64')
        
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
            for i in range(len(self.freq)):
                amp[i] = (A[self.tfpair[j][0],i]+B[self.tfpair[j][0],i])/(2.*A[self.tfpair[j][1],i])
                #amp[i] = (A[self.tfpair[j][0],i]+B[self.tfpair[j][0],i])/(A[self.tfpair[j][1],i]+B[self.tfpair[j][1],i])
            self.tf.append(amp)
            
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
        CORE = zeros((nlayer*2,nlayer*2),dtype='complex64')
        
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
        tf = []
        for tfp in range(ntf):
            tft = np.zeros((fnum,1),dtype='complex64')
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
                tft[nf] = (As[tfpair[tfp][0]*2+1]-As[tfpair[tfp][0]*2])/ \
                    (2.*As[tfpair[tfp][1]*2+1])
                    
            tf.append(tft)
        return tf
        
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
        iCORE = np.zeros((nlayer,nlayer),dtype='complex64')
        
        iD[0] = np.sin(iang)
        iCORE[0,slayer] = 1.
        
        for nl in range(nlayer-1):
            iCORE[nl+1,nl] = 1./vs[nl]
            iCORE[nl+1,nl+1] = -1./vs[nl+1]
            
        iA = solve(iCORE,iD)
        
        iS = np.arcsin(iA)
        
        # Lame Parameter(s)
        mu = np.zeros((nlayer,1),dtype='complex64')
        
        for nl in range(nlayer):
            mu[nl]=dn[nl]*(vs[nl]**2)
        
        # horizontal and vertical slowness
        ns = np.zeros((nlayer,1),dtype='complex64')
        
        for nl in range(nlayer):
            ns[nl]=np.cos(iS[nl])/vs[nl]
            
        # building data vector
        A = np.zeros((nlayer*2,1))
        D = np.zeros((nlayer*2,1))
        D[-1] = 1.
        
        # Dispacement and transfer function initialization
        
        fnum = len(freq)
        
        # loop over frequencies and tfpair
        tf = []
        for tfp in range(ntf):
            tft = np.zeros((fnum,1),dtype='complex64')
            for nf in range(fnum):
                # building core matrix
                CORE = np.zeros((nlayer*2,nlayer*2),dtype='complex64')
                
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
                except ValueError:
                    A[:] = np.nan
                    
                # transfer function
                tft[nf] = (A[tfpair[tfp][0]*2+1]+A[tfpair[tfp][0]*2])/ \
                    (2.*A[tfpair[tfp][1]*2+1])
                    
            tf.append(tft)
        return tf
        
        