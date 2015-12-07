# -*- coding: utf-8 -*-
"""
Created on Mon Dec 07 13:26:28 2015

Display package
to display TF and friends nicely and well done

@author: irnakat
"""
import numpy as np
import pylab as plt
from matplotlib.ticker import ScalarFormatter

def TFPlot(*arg,**kwargs):
    # check given axis name
    try:
        axname = kwargs['axname']
    except KeyError:
        axname = None
    try:
        label = kwargs['label']
    except KeyError:
        label = ['']*(len(arg)/2)
    # create new figure is axis name is not given
    if axname==None:
        f = plt.figure(figsize=(8.,4.))
        a = f.add_subplot(111)
    else:
        a = axname
    
    # set label and properties of axis
    a.set_xlabel('Frequency (Hz)')
    a.set_ylabel('Amplification')
    a.set_yscale('log')
    a.set_xscale('log')
    a.grid(True,which='major',color='k')
    a.grid(True,which='minor',color='grey')
    a.minorticks_on()
    a.tick_params(axis='both', which='major', labelsize=10, labelcolor='k')
    a.tick_params(axis='both', which='minor', labelsize=7.5, labelcolor='grey')
    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(ScalarFormatter())
    
    # check number of input
    if len(arg)%2!=0:
        raise InputError('Number or input pairs if not correct! Detected input pairs is %.1f.'%len(arg)/2.)
    # check length of each pairs and plot data
    minx = []; maxx = []; miny = []; maxy = []
    for i in range(len(arg)/2):
        if len(arg[2*i])!=len(arg[2*i+1]):
            raise LengthError('Length of data on pair number : %d is not the same'%i)
    
        a.plot(arg[2*i],arg[2*i+1],label=label[i])
        minx.append(np.min(arg[2*i]))
        maxx.append(np.max(arg[2*i]))
        miny.append(np.min(arg[2*i+1]))
        maxy.append(np.max(arg[2*i+1]))
        
    minmax = [[np.min(minx),np.max(maxx)],
              [np.min(miny),np.max(maxy)]]
    percspan = 0.1
    minmaxspan = [np.exp(np.log(minmax[0][1]-minmax[0][0])*percspan),
                  np.exp(np.log(minmax[1][1]-minmax[1][0])*percspan)]
    a.set_xlim(minmax[0][0],minmax[0][1])
    a.set_ylim(minmax[1][0]-minmaxspan[1],minmax[1][1]+minmaxspan[1])
    a.legend(loc='Best',fancybox=True,framealpha=0.5)
    
def PhasePlot(*arg,**kwargs):
    # check given axis name
    try:
        axname = kwargs['axname']
    except KeyError:
        axname = None
    try:
        label = kwargs['label']
    except KeyError:
        label = ['']*(len(arg)/2)
    # create new figure is axis name is not given
    if axname==None:
        f = plt.figure(figsize=(8.,4.))
        a = f.add_subplot(111)
    else:
        a = axname
    
    # set label and properties of axis
    a.set_xlabel('Frequency (Hz)')
    a.set_ylabel('Phase (Degrees)')
    #a.set_yscale('log')
    a.set_xscale('log')
    a.grid(True,which='major',color='k')
    a.grid(True,which='minor',color='grey')
    a.minorticks_on()
    a.tick_params(axis='both', which='major', labelsize=10, labelcolor='k')
    a.tick_params(axis='both', which='minor', labelsize=7.5, labelcolor='grey')
    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(ScalarFormatter())
    
    # check number of input
    if len(arg)%2!=0:
        raise InputError('Number or input pairs if not correct! Detected input pairs is %.1f.'%len(arg)/2.)
    # check length of each pairs and plot data
    minx = []; maxx = []; miny = []; maxy = []
    for i in range(len(arg)/2):
        if len(arg[2*i])!=len(arg[2*i+1]):
            raise LengthError('Length of data on pair number : %d is not the same'%i)
            
        y = np.rad2deg(arg[2*i+1])
        a.plot(arg[2*i],y,label=label[i])
        minx.append(np.min(arg[2*i]))
        maxx.append(np.max(arg[2*i]))
        miny.append(np.min(y))
        maxy.append(np.max(y))
        
    minmax = [[np.min(minx),np.max(maxx)],
              [np.min(miny),np.max(maxy)]]
    percspan = 0.1
    minmaxspan = [np.exp(np.log(minmax[0][1]-minmax[0][0])*percspan),
                  (minmax[1][1]-minmax[1][0])*percspan]
    a.set_xlim(minmax[0][0],minmax[0][1])
    a.set_ylim(minmax[1][0]-minmaxspan[1],minmax[1][1]+minmaxspan[1])
    a.legend(loc='Best',fancybox=True,framealpha=0.5)
    
def SpectroPlot(x,y,z,nx=100,ny=100,ylabel='incidence angle',cmap='rainbow'):
    import scipy.interpolate
    
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    
    xi,yi = np.logspace(np.log10(x.min()),np.log10(x.max()),nx), np.linspace(y.min(),y.max(),ny)
    xi,yi = np.meshgrid(xi,yi)
    
    # interpolation
    zi = scipy.interpolate.griddata((x,y),z,(xi,yi),method='linear')
    
    # plot the data
    f = plt.figure(figsize=(8.,4.))
    a = f.add_subplot(111)
    
    am = a.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
             aspect = 'auto')

    a.set_xlabel('Frequency (Hz)')
    a.set_ylabel(ylabel)
    a.set_xscale('log')
    a.minorticks_on()
    a.tick_params(axis='x', which='major', labelsize=10, labelcolor='k')
    a.tick_params(axis='x', which='minor', labelsize=7.5, labelcolor='grey')
    for axis in [a.xaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(ScalarFormatter())
    plt.colorbar(am)


# test display
import IOfile
from TFCalculator import TFCalculator as TFC

fname = 'sampleinput_linear_elastic_6layer_halfspace.dat'
data = IOfile.parsing_input_file(fname)
theclass = TFC(data)
tf1 = theclass.tf_kramer286_sh()

fname2 = 'sampleinput_linear_elastic_1layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
theclass2 = TFC(data2)
tf2 = theclass2.tf_kramer286_sh()

TFPlot(theclass.freq,np.abs(tf1[0]),theclass2.freq,np.abs(tf2[0]),
       label=['line 1','line 2'])

PhasePlot(theclass.freq,np.angle(tf1[0]),theclass2.freq,np.angle(tf2[0]))

fname3 = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
data3 = IOfile.parsing_input_file(fname2)
x = np.array([])
y = np.array([])
z = np.array([])
ianglist = np.linspace(0.0,90.0,91)
for i in range(len(ianglist)):
    data3['iang'] = np.deg2rad(ianglist[i])
    theclass3 = TFC(data3)
    tf3 = theclass3.tf_knopoff_sh_adv()
    #print data3['iang']
    #print np.shape(np.abs(tf3[0][:,0])),np.shape(z)
    x = np.concatenate((x,theclass3.freq))
    z = np.concatenate((z,np.abs(tf3[0][:,0])))
    y = np.concatenate((y,np.zeros_like(theclass3.freq)+ianglist[i]))
    
SpectroPlot(x,y,z,nx=100,ny=100,ylabel='incidence angle',cmap='rainbow')
