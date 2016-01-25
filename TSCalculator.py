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
    def __init__(self,parfile,method='knopoff_sh'):
        self.parfile = parfile
        self.parameters = IOfile.parsing_input_file(parfile)
        if self.parameters['inputmotion'][1]=='ascii':
            self.inp_time,self.inp_signal = IOfile.read_ascii_seismogram(self.parameters['inputmotion'][0])
        self.fs = 1/(self.inp_time[1]-self.inp_time[0])
        self.method = method
        if self.parameters['modeID']==11 or self.parameters['modeID']==12:
            self.linear_equivalent_TF2TS()
        else:
            self.linear_TF2TS()
        
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
        
    def circ_convolution(self,data,tf,disp=False):
        # perform circular convolution on time domain
        # if domain is time, then :
        # data --> input data in time domain
        # IR --> impulse response
        # if domain is frequency, then : 
        # data --> DFT of data
        # IR --> DFT of IR
        """
        tmp = np.real(np.fft.ifft(tf))
        tmp2= np.zeros_like(tmp)
        tmp2[:2000]=tmp[2000:]
        tmp2[2000:]=tmp[:2000]
        fig1 = plt.figure(figsize=(12,8))
        ax = fig1.add_subplot(311)
        ts = np.linspace(-10,10,4000)
        ax.plot(ts,tmp2,label='original impulse response')
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,1.7)
        ax.grid()
        ax.legend()
        
        ax = fig1.add_subplot(312)
        tmp2[:2000]=0.
        ax.plot(ts,tmp2,label='put 0 for negative time')
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,1.7)
        ax.grid()
        ax.legend() 
        
        ax = fig1.add_subplot(313)
        tmp2 = tmp2*(2**2)
        ax.plot(ts,tmp2,label='multiply the signal by $2^2$')
        ax.set_xlim(-0.5,0.5)
        ax.set_ylim(-0.5,1.7)
        ax.grid()
        ax.legend()        
        """
        
        
        newlength=len(data)*2
        dat1 = self.zeropadding(data,newlength)
        dat2 = np.real(np.fft.ifft(tf))*4.
        
        dat2[len(dat2)/2:]=0.

        # highpass filter
        dt = 1./self.fs
        df = 1./((len(dat1)-1)*dt)
        #dat1 = self.butter_highpass_filter(dat1,2.*df,self.fs)
        #dat2 = self.butter_highpass_filter(dat2,2.*df,self.fs)
        
        #out = np.fft.ifft(np.fft.fft(dat1)*tf)
        out = np.real(np.fft.ifft(np.fft.fft(dat1)*np.fft.fft(dat2)))

        # correction of freq vector
        fmax = 1./(2.*dt)
        self.freq = np.linspace(0.,fmax,len(dat1)/2.)
        
        if disp:
            # plot the result(s)
            fig = plt.figure(figsize=(12.,4.))
            a1 = fig.add_subplot(121)
            a1.set_title('Time Series')
            a1.plot(self.inp_time,self.inp_signal,label='input')
            a1.plot(self.inp_time,np.real(out[:len(self.inp_time)]),lw=2.,label='output',alpha=0.5)
            a1.legend(loc='best',fancybox=True,framealpha=0.5)
            a1.grid(True)
            #a1.set_xlim(0,2)
            
            
            # plot the tf abs
            a2 = fig.add_subplot(122)
            a2.set_title('Frequency domain')
            tmp = np.abs(np.fft.fft(dat1))
            a2.plot(self.freq,tmp[:len(self.tfclass.freq)],label='input signal')
            a2.plot(self.tfclass.freq,np.abs(self.tfclass.tf[0]),label='transfer function')
            dftts = np.abs(np.fft.fft(out))
            a2.plot(self.freq,dftts[:len(self.freq)],lw=2.,label='output',alpha=0.5)
            a2.legend(loc='best',fancybox=True,framealpha=0.5)
            a2.set_xscale('log')
            a2.grid(True)  
            
            
            
            # plot the result(s)
            fig = plt.figure(figsize=(12.,4.))
            a1 = fig.add_subplot(121)
            a1.set_title('Time Series')
            a1.plot(self.inp_time,self.inp_signal,label='input')
            a1.plot(self.inp_time,np.real(out[:len(self.inp_time)]),lw=2.,label='output',alpha=0.5)
            a1.legend(loc='best',fancybox=True,framealpha=0.5)
            a1.grid(True)
            #a1.set_xlim(0,2)
            
            
            # plot the tf angle
            a2 = fig.add_subplot(122)
            a2.set_title('Frequency domain')
            tmp = np.angle(np.fft.fft(dat1))
            a2.plot(self.freq,tmp[:len(self.tfclass.freq)],label='input signal')
            a2.plot(self.tfclass.freq,np.angle(self.tfclass.tf[0]),label='transfer function')
            dftts = np.angle(np.fft.fft(out))
            a2.plot(self.freq,dftts[:len(self.freq)],lw=2.,label='output',alpha=0.5)
            a2.legend(loc='best',fancybox=True,framealpha=0.5)
            a2.set_xscale('log')
            a2.grid(True)  
        
        return out
        
    def find_nearest(self,array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx],idx
    
    def linear_TF2TS(self,parameters = [],disp=False):
        if parameters == []: #pure linear calculation
            tfclass = TFC(self.parameters)
            eval('tfclass.tf_'+self.method+'()')
        else:
            tfclass = TFC(parameters)
            eval('tfclass.tf_'+self.method+'()')
        
        # frequency list from original inp_signal
        self.inp_freq = np.fft.fftfreq(len(self.inp_signal),d=self.inp_time[1]-self.inp_time[0])
        
        # data treatment
        self.inp_signal = self.cosine_tapering(self.inp_signal)
        
        # number of tf
        ntf = tfclass.ntf
        self.time_series = []
        self.convolved_horz = []
        self.convolved_vert = []
        self.convolved_horz_resampled = []
        self.convolved_vert_resampled = []
        for ii in range(ntf):
            # circular convolution for horizontal tf
            self.time_series.append(self.circ_convolution(self.inp_signal,self.ft_mirroring(tfclass.tf[ii*2]),disp=disp))
            self.convolved_horz.append(np.fft.fft(self.time_series[-1]))
            self.convolved_horz_resampled.append(self.convolved_horz[-1][::2])
            # circular convolution for vertical tf
            self.time_series.append(self.circ_convolution(self.inp_signal,self.ft_mirroring(tfclass.tf[ii*2+1]),disp=disp))   
            self.convolved_vert.append(np.fft.fft(self.time_series[-1]))
            self.convolved_vert_resampled.append(self.convolved_vert[-1][::2])
        self.freq_resampled = self.freq[::2]
        self.tfclass = tfclass
            
    def linear_equivalent_TF2TS(self):
        """
        Calculation of linear equivalent method for transfer function
        """
        # read G/Gmax and Damping Ratio
        try:
            self.nonlinpar = IOfile.parsing_nonlinear_parameter(self.parameters['GoverGmaxfile'][0])
        except:
            raise KeyError('GoveGmaxfile is not detected! Unable to run linear equivalent calculation!')
        
        # perform sublayer addition
        sublayercriteria = 5. # maximum thickness of layer
        newhl = []; newvs = []; newqs = []; newvp = []; newqp = []; newdn = []; newst = []; nli = []
        for i in range(len(self.parameters['hl'])-1):
            if self.parameters['hl'][i]>sublayercriteria:
                nlayer = np.ceil(self.parameters['hl'][i]/sublayercriteria)
                newhl = np.concatenate((newhl,[self.parameters['hl'][i]/nlayer for j in range(int(nlayer))]))
                newvs = np.concatenate((newvs,[self.parameters['vs'][i] for j in range(int(nlayer))]))
                newqs = np.concatenate((newqs,[self.parameters['qs'][i] for j in range(int(nlayer))]))
                newdn = np.concatenate((newdn,[self.parameters['dn'][i] for j in range(int(nlayer))]))
                newst = np.concatenate((newst,[self.parameters['soiltype'][i] for j in range(int(nlayer))]))
                nli.append(nlayer)
                if self.parameters['modeID']==12:
                    newvp = np.concatenate((newvp,[self.parameters['vp'][i] for j in range(int(nlayer))]))
                    newqp = np.concatenate((newqp,[self.parameters['qp'][i] for j in range(int(nlayer))]))
            else:
                newhl = np.concatenate((newhl,[self.parameters['hl'][i]]))
                newvs = np.concatenate((newvs,[self.parameters['vs'][i]]))
                newqs = np.concatenate((newqs,[self.parameters['qs'][i]]))
                newdn = np.concatenate((newdn,[self.parameters['dn'][i]]))
                nli.append(1.)
                newst = np.concatenate((newst,[self.parameters['soiltype'][i]]))
                if self.parameters['modeID']==12:
                    newvp = np.concatenate((newvp,[self.parameters['vp'][i]]))
                    newqp = np.concatenate((newqp,[self.parameters['qp'][i]]))
        # for the last layer
        newhl = np.concatenate((newhl,[0]))
        newvs = np.concatenate((newvs,[self.parameters['vs'][-1]]))
        newqs = np.concatenate((newqs,[self.parameters['qs'][-1]]))
        newdn = np.concatenate((newdn,[self.parameters['dn'][-1]]))
        newst = np.concatenate((newst,[self.parameters['soiltype'][-1]]))
        nli.append(1.)
        if self.parameters['modeID']==12:
            newvp = np.concatenate((newvp,[self.parameters['vp'][-1]]))
            newqp = np.concatenate((newqp,[self.parameters['qp'][-1]]))
        # assign sublayer to parameter
        self.parameters['nlayer']=len(newhl)
        self.parameters['hl']=newhl.tolist()
        self.parameters['vs']=newvs.tolist()
        self.parameters['qs']=newqs.tolist()
        self.parameters['dn']=newdn.tolist()
        self.parameters['soiltype']=newst.tolist()
        if self.parameters['modeID']==12:
            self.parameters['vp']=newvp.tolist()
            self.parameters['qp']=newqp.tolist()
        # correction of tfpair
        oldpair = self.parameters['tfPair']  
        for i in range(len(oldpair)):
            oldpair[i][0] = int(np.sum(nli[:oldpair[i][0]]))
            oldpair[i][1] = int(np.sum(nli[:oldpair[i][1]]))
        # correction of source location
        self.parameters['sourceloc']=int(np.sum(nli[:self.parameters['sourceloc']]))
        # modification of tfpair
        newpair = [[i,i+1] for i in range(self.parameters['nlayer']-1)]
        self.parameters['tfPair'] = newpair
        self.parameters['ntf'] = len(newpair)
        
        # initial calculation G/Gmax = 1, damping ratio = 1st value
        Gmax = [self.parameters['dn'][i]*self.parameters['vs'][i]**2 for i in range(len(newpair))]
        self.linear_TF2TS(parameters=self.parameters)
        
        # check non linearity and update parameter
        for i in range(len(newpair)):
            val,idx =  self.find_nearest(self.nonlinpar['nonlin strain'][int(self.parameters['soiltype'][i])],np.max(self.time_series[i*2]))
            self.parameters['vs'][i]=np.sqrt((Gmax[i]*self.nonlinpar['nonlin G/Gmax'][int(self.parameters['soiltype'][i])][idx])/self.parameters['dn'][i])
            self.parameters['qs'][i]=1./(2.*self.nonlinpar['nonlin damping'][int(self.parameters['soiltype'][i])][idx]/100.)
            
        numiter = 3
        for j in range(numiter):     
            self.linear_TF2TS(parameters=self.parameters)
            # check non linearity and update parameter
            for i in range(len(newpair)):
                val,idx =  self.find_nearest(self.nonlinpar['nonlin strain'][int(self.parameters['soiltype'][i])],np.max(self.time_series[i*2]))
                self.parameters['vs'][i]=np.sqrt((Gmax[i]*self.nonlinpar['nonlin G/Gmax'][int(self.parameters['soiltype'][i])][idx])/self.parameters['dn'][i])
                self.parameters['qs'][i]=1./(2.*self.nonlinpar['nonlin damping'][int(self.parameters['soiltype'][i])][idx]/100.)
        # correction of tfpair
        self.parameters['tfPair']=oldpair
        self.parameters['ntf']=len(oldpair)
        self.linear_TF2TS(parameters=self.parameters)

        
# test calculation
method = 'knopoff_sh'
TSclass00 = TSCalculator('sampleinput_linear_elastic_1layer_halfspace.dat',method)
TSclass01 = TSCalculator('sampleinput_linear_eq_elastic_1layer_halfspace.dat',method)

#plot the outputs
dt = TSclass01.inp_time[1]-TSclass01.inp_time[0]
fig = plt.figure(figsize=(12,8))
a = fig.add_subplot(211)
a.plot(TSclass01.inp_time,TSclass01.inp_signal,label='in')
a.plot([i*dt for i in range(len(TSclass00.time_series[0]))],TSclass00.time_series[0],label='out linear')
a.plot([i*dt for i in range(len(TSclass01.time_series[0]))],TSclass01.time_series[0],label='out lin-eq')
a.legend(loc='best')
a.set_xlim(9.5,12.5)
a.grid(True)

a = fig.add_subplot(212)
a.plot(np.fft.fftfreq(len(TSclass01.inp_signal),dt)[:len(TSclass01.inp_signal)/2],np.abs(np.fft.fft(TSclass01.inp_signal))[:len(TSclass01.inp_signal)/2],label='in')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass00.time_series[0])),dt)[:len(TSclass00.time_series[0])/2],np.abs(np.fft.fft(TSclass00.time_series[0]))[:len(TSclass00.time_series[0])/2],label='out lin')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass01.time_series[0])),dt)[:len(TSclass01.time_series[0])/2],np.abs(np.fft.fft(TSclass01.time_series[0]))[:len(TSclass01.time_series[0])/2],label='out lin-eq')
a.plot(TSclass00.tfclass.freq,np.abs(TSclass00.tfclass.tf[0])*0.1,label='(transfer function)/10. lin')
a.plot(TSclass01.tfclass.freq,np.abs(TSclass01.tfclass.tf[0])*0.1,label='(transfer function)/10. lin-eq')
a.legend(loc='best')
a.set_xscale('log')
a.set_xlim(0.1,50)
a.grid(True)