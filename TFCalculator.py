# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 20:04:33 2015
Program to calculate transfer function using several different methods
@author: irnakat
"""
import commondata as cd
import numpy as np

class TFCalculator:
    def __init__(self,data):
        self.mode = data[0]         # calculation mode
        self.ntf = data[1]          # number of transfer function pairs
        self.tfpair = data[2]       # list of transfer function pairs
        self.nlayer = data[3]       # number of soil layer
        self.hl = data[4][0]        # thickness of each layer
        self.vs = data[4][1]        # shear wave velocity
        self.dn = data[4][2]        # density
        self.qs = data[4][3]        # quality factor
        self.freq = np.linspace(0.01,100.,200)
        
    def kramer286(self):
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
        Gc = np.zeros((len(self.hl),len(self.f)),dtype='complex64')
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
                A[i+1,j] = 0.5*A[i,j]*(1.0+alphaz) * np.exp(Yi*alpha*h[i]) + \
                    0.5*B[i,j]*(1.-alphaz)*np.exp(-Yi*alpha*h[i])
                B[i+1,j] = 0.5*A[i,j]*(1.0-alphaz) * np.exp(Yi*alpha*h[i]) + \
                    0.5*B[i,j]*(1.+alphaz)*np.exp(-Yi*alpha*h[i])
        # calculating transfer function
        tf = []
        for j in range(self.ntf):
            for i in range(len(self.freq)):
                amp[i] = (A[lt,i]+B[lt,i])/(A[lb,i]+B[lb,i])
            tf.append(amp)
            
        return tf
        
    def sh_tf_knopoff(self):
        