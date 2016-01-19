# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:06:03 2016

@author: irnakat
"""
from TFCalculator import TFCalculator as TFC
import numpy as np
import IOfile
import scipy.signal as signal

# temporary import for debugging
import pylab as plt

class TSCalculator:
    """
    class which manages calculation of TimeSeries (TS) given the Transfer Function (TF)
    """
    def __init__(self,tfclass):
        self.tfclass = tfclass
        if self.tfclass.inputmotion[1]=='ascii':
            self.inp_time,self.inp_signal = IOfile.read_ascii_seismogram(self.tfclass.inputmotion[0])
        self.fs = 1/(self.inp_time[1]-self.inp_time[0])
        self.parameters = IOfile.parsing_input_file(self.tfclass.parfile)
        
    def cosine_tapering(self,inp_signal,alpha=0.10):
        window = signal.tukey(len(inp_signal),alpha=alpha,sym=True)
        return inp_signal*window        
        
    def zeropadding(self,inp,newlength):
        # create zero padding as the size of the original signal
        if newlength>len(inp):
            out = np.zeros((newlength),dtype='complex128')
            out[:len(inp)] = inp
        else:
            out = inp
        return out
    
    def butter_highpass_filter(self,data, lowcut, fs, order=4):
        # perform highpass filter of the signal
        def butter_highpass(lowcut, fs, order=4):
            nyq = 0.5 * fs
            low = lowcut / nyq
            b, a = signal.butter(order, [low], btype='highpass')
            return b, a
        b, a = butter_highpass(lowcut, fs, order=order)
        return signal.lfilter(b, a, data)
        
    def butter_lowpass_filter(self,data, highcut, fs, order=4):
        # perform highpass filter of the signal
        def butter_lowpass(highcut, fs, order=4):
            nyq = 0.5 * fs
            high = highcut / nyq
            b, a = signal.butter(order, [high], btype='lowpass')
            return b, a
        b, a = butter_lowpass(highcut, fs, order=order)
        return signal.lfilter(b, a, data)
        
    def ft_mirroring(self,tf,mirror=True):
        if mirror:
            tmp = tf
            tf = np.zeros((len(tmp)*2),dtype='complex128')
            tf[:len(tmp)] = tmp
            tf[len(tmp):2*len(tmp)] = np.flipud(tmp)
            return tf
        else:
            return tf[:len(tf)/2]
        
    def circ_convolution(self,data,IR,domain='time',disp=False):
        # perform circular convolution on time domain
        # if domain is time, then :
        # data --> input data in time domain
        # IR --> impulse response
        # if domain is frequency, then : 
        # data --> DFT of data
        # IR --> DFT of IR
        
        if domain!='time' or domain!='frequency':
            domain= 'time'
        
        # determine final length
        if domain=='frequency':
            data = np.fft.ifft(data)
            IR = np.fft.ifft(IR)
        
        newlength=len(data)+len(IR)
        dat1 = self.zeropadding(data,newlength)
        dat2 = self.zeropadding(IR,newlength)
        
        # highpass filter
        dt = 1./self.fs
        df = 1./((len(dat1)-1)*dt)
        #dat1 = self.butter_highpass_filter(dat1,2.*df,self.fs)
        #dat2 = self.butter_highpass_filter(dat2,2.*df,self.fs)
        
        out = np.fft.ifft(np.fft.fft(dat1)*np.fft.fft(dat2))

        # correction of freq vector
        fmax = 1./(2.*dt)
        self.freq = np.linspace(0.,fmax,len(dat1)/2.)
        
        if disp:
            # plot the result(s)
            fig = plt.figure(figsize=(12.,4.))
            a1 = fig.add_subplot(121)
            a1.set_title('Time Series')
            a1.plot(self.inp_time,self.inp_signal,label='input')
            a1.plot(self.inp_time,out[:len(self.inp_time)],lw=2.,label='output',alpha=0.5)
            a1.legend(loc='best',fancybox=True,framealpha=0.5)
            a1.grid(True)
            #a1.set_xlim(0,2)
            
            
            # plot the tf
            a2 = fig.add_subplot(122)
            a2.set_title('Frequency domain')
            tmp = np.abs(np.fft.fft(self.inp_signal))
            a2.plot(self.tfclass.freq,tmp[:len(self.tfclass.freq)],label='input signal')
            a2.plot(self.tfclass.freq,np.abs(self.tfclass.tf[0]),label='transfer function')
            dftts = np.abs(np.fft.fft(out))
            a2.plot(self.freq,dftts[:len(self.freq)],lw=2.,label='output',alpha=0.5)
            xx = self.freq[::2]
            yy = dftts[::2]
            a2.plot(xx,yy[:len(xx)],lw=1.5,label='output-corrected')
            a2.legend(loc='best',fancybox=True,framealpha=0.5)
            a2.set_xscale('log')
            a2.grid(True)        
        
        return out
    
    def linear_TF2TS(self):
        # frequency list from original inp_signal
        self.inp_freq = np.fft.fftfreq(len(self.inp_signal),d=self.inp_time[1]-self.inp_time[0])
        
        # data treatment
        self.inp_signal = self.cosine_tapering(self.inp_signal)
        
        # number of tf
        ntf = self.tfclass.ntf
        self.time_series = []
        self.convolved_horz = []
        self.convolved_vert = []
        self.convolved_horz_resampled = []
        self.convolved_vert_resampled = []
        self.freq_resampled = self.freq[::2]
        for ii in range(ntf):
            # circular convolution for horizontal tf
            self.time_series.append(self.circ_convolution(self.inp_signal,np.fft.ifft(self.ft_mirroring(self.tfclass.tf[ii*ntf])),disp=True))
            self.convolved_horz.append(np.fft.fft(self.time_series[-1]))
            self.convolved_horz_resampled(self.convolved_horz[-1][::2])
            # circular convolution for vertical tf
            self.time_series.append(self.circ_convolution(self.inp_signal,np.fft.ifft(self.ft_mirroring(self.tfclass.tf[ii*ntf+1])),disp=True))   
            self.convolved_vert.append(np.fft.fft(self.time_series[-1]))
            self.convolved_vert_resampled(self.convolved_vert[-1][::2])

"""
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
# input file reading
datash = IOfile.parsing_input_file(fname)
# kramer sh
print 'TF calculation using kramer approach'
theclass1 = TFC(datash)
theclass1.tf_kramer286_sh()
print 'calculation has been finished!'

TSclass01 = TSCalculator(theclass1)
TSclass01.linear_TF2TS()

# filename
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
# input file reading
datash = IOfile.parsing_input_file(fname)
# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass1 = TFC(datash)
theclass1.tf_knopoff_sh()
print 'calculation has been finished!'

TSclass01 = TSCalculator(theclass1)
TSclass01.linear_TF2TS()

# filename
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
# input file reading
datash = IOfile.parsing_input_file(fname)
# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass1 = TFC(datash)
theclass1.tf_kennet_sh()
print 'calculation has been finished!'

TSclass01 = TSCalculator(theclass1)
TSclass01.linear_TF2TS()

fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
# input file reading
datapsvs = IOfile.parsing_input_file(fname2)
# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass2 = TFC(datapsvs)
theclass2.tf_knopoff_psv_adv()
print 'calculation has been finished!'

TSclass02 = TSCalculator(theclass2)
TSclass02.linear_TF2TS()

fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
# input file reading
datapsvs = IOfile.parsing_input_file(fname2)
datapsvs['iang']=30.
# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass2 = TFC(datapsvs)
theclass2.tf_knopoff_psv_adv()
print 'calculation has been finished!'

TSclass02 = TSCalculator(theclass2)
TSclass02.linear_TF2TS()
"""