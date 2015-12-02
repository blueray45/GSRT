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

        if self.mode==cd.mode[4] or self.mode==cd.mode[5] or \
            self.mode==cd.mode[6] or self.mode==cd.mode[7]:
                
            self.vp = np.array(data['vp'])
            self.qp = np.array(data['qp'])
            self.comp = np.array(data['comp'])
        
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
        print 'angle of incidence '+str(np.degrees(self.iang))
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
        
        iD = np.zeros((nlayer*2),dtype='complex64')
        iCORE = np.zeros((nlayer*2,nlayer*2),dtype='complex64')
        
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
            
        iA = solve(iCORE,iD)

        iP = np.zeros((nlayer),dtype='complex64')
        iS = np.zeros((nlayer),dtype='complex64')
        
        for i in range(nlayer):
            row = nl*2
            iP[nl] = np.arcsin(iA[row])            
            iS[nl] = np.arcsin(iA[row+1])
        
        # Lame Parameter(s)
        mu = np.zeros((nlayer),dtype='complex64')
        la = np.zeros((nlayer),dtype='complex64')
        
        for nl in range(nlayer):
            mu[nl]=dn[nl]*(vs[nl]**2)
            la[nl]=dn[nl]*((vp[nl]**2)-(2.*(vs[nl]**2)))            
        
        # horizontal and vertical slowness
        rp = np.zeros((nlayer),dtype='complex64')
        rs = np.zeros((nlayer),dtype='complex64')
        np1 = np.zeros((nlayer),dtype='complex64')
        ns = np.zeros((nlayer),dtype='complex64')
        
        for nl in range(nlayer):
            rp[nl] = np.sin(iP[nl])/vp[nl]
            rs[nl] = np.sin(iS[nl])/vs[nl]
            np1[nl] = np.cos(iP[nl])/vp[nl]
            ns[nl] = np.cos(iS[nl])/vs[nl]
        
        # building data vector
        A = np.zeros((nlayer*4),dtype='complex64')
        D = np.zeros((nlayer*4),dtype='complex64')
        
        if comp=='p':
            D[-2] = vp[-1]
        elif comp=='s':
            D[-1] = vs[-1]
        
        # Dispacement and transfer function initialization
        
        fnum = len(freq)
        
        # loop over frequencies and tfpair
        htf = []; vtf = []
        for tfp in range(ntf):
            htft = np.zeros((fnum),dtype='complex64')
            vtft = np.zeros((fnum),dtype='complex64')
            for nf in range(fnum):
                # building core matrix
                CORE = np.zeros((nlayer*4,nlayer*4),dtype='complex64')
                
                # free surface constraints
                CORE[0,0] = la[0]*rp[0]**2+(la[0]+2.*mu[0])*np1[0]**2
                CORE[0,1] = 2.*mu[0]*rs[0]*ns[0]
                CORE[0,2] = la[0]*rp[0]**2+(la[0]+2.*mu[0])*np1[0]**2
                CORE[0,3] = -2.*mu[0]*rs[0]*ns[0]
                
                CORE[1,0] = -2.*mu[0]*rp[0]*np1[0]
                CORE[1,1] = -mu[0]*(rs[0]**2-ns[0]**2)
                CORE[1,2] = 2.*mu[0]*rp[0]*np1[0]
                CORE[1,3] = -mu[0]*(rs[0]**2-ns[0]**2)

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
                    A=solve(CORE,D)
                except ValueError:
                    A[:] = np.nan
                    
                # transfer function
                htft[nf] = (rp[tfpair[tfp][0]]*(A[tfpair[tfp][0]+2]+A[tfpair[tfp][0]+0])+ \
                            ns[tfpair[tfp][0]]*(A[tfpair[tfp][0]+3]-A[tfpair[tfp][0]+1]))/ \
                            (2.*(rp[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+2])+ \
                            ns[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+3])))
                vtft[nf] = (np1[tfpair[tfp][0]]*(-A[tfpair[tfp][0]+2]+A[tfpair[tfp][0]+0])+ \
                            rs[tfpair[tfp][0]]*(A[tfpair[tfp][0]+3]+A[tfpair[tfp][0]+1]))/ \
                            (2.*(-np1[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+2])+ \
                            rs[tfpair[tfp][1]]*(A[4*tfpair[tfp][1]+3])))
                    
            htf.append(htft)
            vtf.append(vtft)
        return htf,vtf
    
    def tf_kennett(self):
        """
        porting from geopsy based on kennet formalism
        Pierre-Yves Bard (LGIT, Grenoble, France)
        Cecile Cornou (LGIT, Grenoble, France), minor modifications
        Marc Wathelet (LGIT, Grenoble, France), port to subroutine
        Theodosius Marwan Irnaka (GEM, Pavia, Italy), port to python
        """
        
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
        
        # initialization for response spectrum calculation
        ai = np.complex(0.,1.)   
        period = 1./freq
        
        freq0 = 10.         # reference frequency for Q <-- WHY?
        eq = 0.             # ????
        fmin = min(freq)    # use min(freq) instead of 0 on the original code
        fmax = max(freq)    # 
        q = 1.e+20          # ????
        iang = self.iang    # incidence angle (radians)
        mode = 3            # for SH case <-- try for mode 1 and 2 for PSV case later
        jcas = 0            # jcas
        nr = len(tfpair)    # number of receiver given by the pair of input and output motion
        zr = np.array([tfpair[i][0]*hl for i in range(len(tfpair))])
                            # depth of receiver given as list (BEWARE OF INCOMPATIBILITY!!)
        
        # defining depth of the top layer (can be modified later)
        th = np.zeros(nlayer+1)
        for i in range(1,nlayer+1):
            th[i] += hl[i-1]
            
        # portion of signal definition and calculationg of Fourier Transform
            # omitted
        
        nf = len(freq)
        nf1=nf+1
        dfreq=freq[1]-freq[0]
        aw = -np.pi/q           # I don't understand! it's basically 0!
        om0 = 2.*np.pi*freq0    # omega at reference frequency for Q
        
        # index of receiver
        izr = np.array([tfpair[i][1] for i in range(len(tfpair))])
        
        # iterating over frequencies
        for i,fr in enumerate(freq):
            fr = 0.05*fr if fr==0. else fr  # correction for zero frequency
            rw = fr*np.pi*2.
            omega = np.complex(rw,aw)
            omm = np.sqrt(rw**2+aw**2)
            adr = np.log(omm/om0)
            if rw!=0.:
                phi = np.arctan2(rw,aw)
            else:
                phi = np.pi/2. if aw>0. else -np.pi/2.
            ad = adr + ai*phi
            df = fr - freq0
            df = 0. if df<0. else df
            
            # iterating over layers
            cvp = np.zeros_like(vp,dtype='complex64')              
            cvs = np.zeros_like(cvp)
            wb = np.zeros_like(cvp)
            for li in range(nlayer):
                qsi = qs[li]*(1.+eq*df)
                piqs = 1./(np.pi*qsi)
                cqs = 1. - 0.5*ai/qsi
                cds = 1. - piqs*ad
                cvs[li] = vs[li]/(cds*cqs)      # complex Vs
                qpi=qp[li]*(1.+eq*df)
                piqp = 1./(np.pi*qpi)
                cqp = 1.-0.5*ai/qpi
                cdp = 1.-piqp*ad
                cvp[li] = vp[li]/(cdp*cqp)      # complex Vp
                wb[li] = (omega/vs[li])**2      # omega for vs
                wa[li] = (omega/vp[li])**2      # oemga for vp
            
            if modeID==9:
                wx0 = wa[nlayer]*np.sin(iang)
                mx0 = n*fr*np.sin(iang)/vp[nlayer]
                c1 = -ai/wa[nlayer]
            else:
                wx0 = wb[nlayer]*np.sin(iang)
                mx0 = n*fr*np.sin(iang)/vs[nlayer]
                c1 = ai/wb[nlayer]
            wx02 = wx0**2
            
            print cvp,cvs
            
            # paused to built another function
        
        def kenpsv(jcas,ckx,comega,nlayer,th,dn,cvp,cvs,cfwave):
            """
            kennet formalism
            """
            from copy import deepcopy
            
            # if th[0]=0, given the depth of the interfaces, else given the thickness
            hc = np.zeros(nlayer)            
            if th[0]==0.:                
                for ic in range(nlayer):
                    hc[ic]=th[ic+1]-th[ic]
            else:
                hc = th
            
            cwx = deepcopy(ckx)
            cwx2= cwx**2
            omega = deepcopy(comega)
            omega2= omega**2
            
            

# debugging
import IOfile
fname2 = 'sampleinput_psv_p_linear_elastic_1layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
theclass2 = TFCalculator(data2)
theclass2.tf_kennett()