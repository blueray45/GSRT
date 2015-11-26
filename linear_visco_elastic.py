# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:46:29 2015
@author: irnakat
"""
# kramer268
import numpy as np
import pylab as plt
def linearTF(f,Vs,rho,h,Qs,lt,lb):
    """
    f   : frequency array
    Vs  : Shear wave velocity
    rho : density
    h   : thickness of layers
    Qs  : Quality factor of Shear wave velocity
    lt  : index of top layer
    lb  : index of bottom layer
    """
    Gc = np.zeros((len(h),len(f)),dtype='complex64')
    A = np.zeros_like(Gc)
    B = np.zeros_like(Gc)
    eta = np.zeros_like(Qs)
    amp = np.zeros((len(f)),dtype='complex64')
    
    Yi = np.complex(0.,1.)
    A[0,:] = np.complex(1.,0.)
    B[0,:] = np.complex(1.,0.)
    ai = 0.0
    # calculating complex modulus
    for i in range(len(h)):
        eta[i] = 1./(2.*Qs[i])
        for j in range(len(f)):
            Gc[i,j] = rho[i]*Vs[i]**2 * np.complex(1.,2.*eta[i])
    # calculating propagator terms (A and B)
    for j in range(len(f)):
        for i in range(len(h)-1):
            fc = np.complex(f[j],-ai)
            alphaz = np.sqrt((rho[i]*Gc[i,j])/(rho[i+1]*Gc[i+1,j]))
            alpha = 2.*np.pi*np.sqrt(rho[i]/Gc[i,j])*fc
            A[i+1,j] = 0.5*A[i,j]*(1.0+alphaz) * np.exp(Yi*alpha*h[i]) + \
                0.5*B[i,j]*(1.-alphaz)*np.exp(-Yi*alpha*h[i])
            B[i+1,j] = 0.5*A[i,j]*(1.0-alphaz) * np.exp(Yi*alpha*h[i]) + \
                0.5*B[i,j]*(1.+alphaz)*np.exp(-Yi*alpha*h[i])
    # calculating transfer function
    for i in range(len(f)):
        amp[i] = (A[lt,i]+B[lt,i])/(A[lb,i]+B[lb,i])
        
    return amp
    
# parameter
h = np.array([100.0,1000.0])
rho = np.array([1500.,2500.])
Vs = np.array([50.,1000.])
Qs = np.array([500000.,500000.])
f = np.linspace(0.1,100.,100)
lt = 0
lb = 1


amp = linearTF(f,Vs,rho,h,Qs,lt,lb)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(f,np.abs(amp))
ax.axis('tight')

# parameter
h = np.array([100.0,1000.0])
rho = np.array([1500.,2500.])
Vs = np.array([50.,1000.])
Qs = np.array([5000.,5000.])
f = np.linspace(0.1,100.,100)
lt = 0
lb = 1


amp = linearTF(f,Vs,rho,h,Qs,lt,lb)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(f,np.abs(amp))
ax.axis('tight')