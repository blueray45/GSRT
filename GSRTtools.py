# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 04:42:31 2016

@author: irnakat
"""
import numpy as np

def disp2vel(data,dt):
    datafft = np.fft.fft(data)
    datafreq= np.fft.fftfreq(len(data),dt)
    #datafreq = [(i+1)*datafreq[1] for i in range(len(datafft))]
    return np.real(np.fft.ifft(datafft*(2.*np.pi*np.abs(datafreq))))

def vel2acc(data,dt):
    datafft = np.fft.fft(data)
    datafreq= np.fft.fftfreq(len(data),dt)
    #datafreq = [(i+1)*datafreq[1] for i in range(len(datafft))]
    return np.real(np.fft.ifft(datafft*(2.*np.pi*np.abs(datafreq))))
    
def acc2vel(data,dt):
    datafft = np.fft.fft(data)
    datafreq= np.fft.fftfreq(len(data),dt)
    datafreq = [(i+1)*datafreq[1] for i in range(len(datafft))]
    return np.real(np.fft.ifft(datafft/(2.*np.pi*np.abs(datafreq))))*2.

def vel2disp(data,dt):
    datafft = np.fft.fft(data)
    datafreq= np.fft.fftfreq(len(data),dt)
    datafreq = [(i+1)*datafreq[1] for i in range(len(datafft))]
    return np.real(np.fft.ifft(datafft/(2.*np.pi*np.abs(datafreq))))*2
    
def acc2disp(data,dt):
    return vel2disp(acc2vel(data,dt),dt)
    
def disp2acc(data,dt):
    return vel2acc(disp2vel(data,dt),dt)
    
#
#import IOfile
#import pylab as plt
#from copy import deepcopy
#time,data = IOfile.read_ascii_seismogram('ferndale2.asc')
#
#disp1 = deepcopy(data)
#vel1 = disp2vel(data,time[1]-time[0])
#acc1 = vel2acc(vel1,time[1]-time[0])
#
#fig = plt.figure(figsize=(12,10))
#a1 = fig.add_subplot(311)
#a1.plot(time,data)
##a1.set_xlim(9.5,10.5)
#a1.set_ylabel('$m$')
#a1.set_title('displacement')
#
#a2 = fig.add_subplot(312)
#a2.plot(time,vel1)
##a2.set_xlim(9.5,10.5)
#a2.set_ylabel('$m/s$')
#a2.set_title('velocity')
#
#a3 = fig.add_subplot(313)
#a3.plot(time,acc1)
##a3.set_xlim(9.5,10.5)
#a3.set_ylabel('$m/s^2$')
#a3.set_xlabel('s')
#a3.set_title('acceleration')
#
#time,data = IOfile.read_ascii_seismogram('ferndale2.asc')
#data = np.array(data)/100.
##data= deepcopy(acc1)
#vel = acc2vel(data,time[1]-time[0])
#disp = vel2disp(vel,time[1]-time[0])
#disp2= acc2disp(data,time[1]-time[0])
#
#fig = plt.figure(figsize=(12,10))
#a1 = fig.add_subplot(311)
#a1.plot(time,disp)
##a1.set_xlim(9.5,10.5)
#a1.set_ylabel('$m$')
#a1.set_title('displacement')
#
#a2 = fig.add_subplot(312)
#a2.plot(time,vel)
##a2.set_xlim(9.5,10.5)
#a2.set_ylabel('$m/s$')
#a2.set_title('velocity')
#
#a3 = fig.add_subplot(313)
#a3.plot(time,data)
##a3.set_xlim(9.5,10.5)
#a3.set_ylabel('$m/s^2$')
#a3.set_xlabel('s')
#a3.set_title('acceleration')