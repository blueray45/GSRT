# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:06:03 2016

@author: irnakat
"""
from TFCalculator import TFCalculator as TFC
import numpy as np
import IOfile
import scipy.signal as signal
import GSRTtools as GT
from copy import deepcopy

# temporary import for debugging
import pylab as plt

class TSCalculator:
    """
    class which manages calculation of TimeSeries (TS) given the parameter file and input file
    """
    def __init__(self,parfile,method='auto',sublayercriteria = 5.,numiter = 10,conv_level = 0.01,verbose=False):
        # parameter file initialization
        self.parfile = parfile
        # read file parameter
        self.parameters = IOfile.parsing_input_file(parfile)
        # method is automatically defined
        if method=='auto':
            if self.parameters['type']=='PSV':
                self.method = 'knopoff_psv_adv'
            else:
                if self.parameters['nlayer']<=5:
                    if self.parameters['iang']==0.:
                        self.method='knopoff_sh'
                    else:
                        self.method='knopoff_sh_adv'
                else:
                    self.method = 'kennet_sh'
        # checking input file
        if self.parameters['inputmotion'][1]=='ascii':
            self.inp_time,self.inp_signal = IOfile.read_ascii_seismogram(self.parameters['inputmotion'][0])
        else:
            raise KeyError('Input motion other than ascii format is not yet supported! Please convert it to displacement on another software first!')
        self.inp_signal = [i/100. for i in self.inp_signal]
        self.fs = 1/(self.inp_time[1]-self.inp_time[0])
        self.dt = self.inp_time[1]-self.inp_time[0]
        self.df = 1./((len(self.inp_signal)-1)*self.dt)
        
        # baseline correction for input signal
        self.inp_signal = self.inp_signal-np.mean(self.inp_signal)        
        
#        if self.parameters['inputmotion'][2]=='vel':
#            self.inp_signal = GT.vel2disp(self.inp_signal,self.dt)
#            self.inp_signal = self.cosine_tapering(self.inp_signal)
#            self.inp_signal = self.butter_highpass_filter(self.inp_signal,2.*self.df,self.fs)
#        elif self.parameters['inputmotion'][2]=='acc':
#            self.inp_signal = GT.acc2disp(self.inp_signal,self.dt)
#            self.inp_signal = self.cosine_tapering(self.inp_signal)
#            self.inp_signal = self.butter_highpass_filter(self.inp_signal,2.*self.df,self.fs)
            
        self.method = method
        if self.parameters['modeID']==11 or self.parameters['modeID']==12:
            self.linear_equivalent_TF2TS(sublayercriteria,numiter,conv_level,verbose)
        else:
            self.linear_TF2TS()
            self.lastiter = 1
        
        if self.parameters['inputmotion'][2]=='acc':
            for i in range(len(self.time_series)):
                self.time_series[i] = GT.disp2acc(self.time_series[i],self.dt)
        elif self.parameters['inputmotion'][2]=='vel':
            for i in range(len(self.time_series)):
                self.time_series[i] = GT.disp2vel(self.time_series[i],self.dt)
        
    def cosine_tapering(self,inp_signal,alpha=0.10):
        #window = signal.tukey(len(inp_signal),alpha=alpha,sym=True)
        window = np.zeros_like(inp_signal)
        N = len(inp_signal)
        for i in range(N):
            if i < alpha*(N-1)/2:
                window[i]=0.5*(1.+np.cos(np.pi*(((2.*i)/(alpha*(N-1)))-1.)))
            elif i > (N-1)*(1.-alpha/2):
                window[i]=0.5*(1.+np.cos(np.pi*(((2.*i)/(alpha*(N-1)))-((2./alpha)+1.))))
            else:
                window[i]=1.0
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
        df = 1./((1./self.fs)*len(tf))
        datafreq = [(i+1)*df for i in range(len(tf))]
        if self.parameters['inputmotion'][2]=='acc':
            tf = GT.acc2dispfreq(tf,datafreq)
        elif self.parameters['inputmotion'][2]=='vel':
            tf = GT.vel2dispfreq(tf,datafreq)
        
        newlength=len(data)*2
        dat1 = self.zeropadding(data,newlength)
        dat2 = np.real(np.fft.ifft(tf))*4.
        
        dat2[len(dat2)/2:]=0.
        
        # highpass filter
        dt = 1./self.fs
        df = 1./((len(dat1)-1)*dt)
        dat1 = self.butter_highpass_filter(dat1,2.*df,self.fs)
        dat2 = self.butter_highpass_filter(dat2,2.*df,self.fs)
        
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
        magnitude = 7.5
        ratio = (magnitude-1.)/10.
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
            self.time_series.append(self.circ_convolution(self.inp_signal,self.ft_mirroring(tfclass.tf[ii*2]),disp=disp)*ratio)
            self.convolved_horz.append(np.fft.fft(self.time_series[-1]))
            self.convolved_horz_resampled.append(self.convolved_horz[-1][::2])
            # circular convolution for vertical tf
            self.time_series.append(self.circ_convolution(self.inp_signal,self.ft_mirroring(tfclass.tf[ii*2+1]),disp=disp)*ratio)   
            self.convolved_vert.append(np.fft.fft(self.time_series[-1]))
            self.convolved_vert_resampled.append(self.convolved_vert[-1][::2])
        self.freq_resampled = self.freq[::2]
        self.tfclass = tfclass
            
    def linear_equivalent_TF2TS(self,sublayercriteria = 5.,numiter = 10,conv_level = 0.01,verbose=False):
        """
        Calculation of linear equivalent method for transfer function
        """
        # set up linear equivalent parameters
        #sublayercriteria = 5.   # maximum thickness of layer
        #numiter = 10            # maximum number of iteration
        #conv_level = 0.01       # convergence limit
        
        # read G/Gmax and Damping Ratio
        try:
            self.nonlinpar = IOfile.parsing_nonlinear_parameter(self.parameters['GoverGmaxfile'][0])
        except:
            raise KeyError('GoveGmaxfile is not detected! Unable to run linear equivalent calculation!')
        
        # perform sublayer addition
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
        Gerr = np.zeros((len(newpair)))
        Derr = np.zeros((len(newpair)))
        
        # print original model information
        if verbose:
            print 'iteration %3d'%0
            print 'Thickness\tVs\t Qs\t convG\t convD'
            for i in range(len(newpair)):
                print '%.2f\t\t%.2f\t%.2f\tn/a\tn/a'%(self.parameters['hl'][i],self.parameters['vs'][i],self.parameters['qs'][i])        
        
        for j in range(numiter):     
            if verbose:
                print 'iteration %3d'%(j+1)
                print 'Thickness\tVs\t Qs\t convG\t convD'
            self.linear_TF2TS(parameters=self.parameters)
            # check non linearity and update parameter
            for i in range(len(newpair)):
                # retrieve old G and D
                oldG = self.parameters['dn'][i]*self.parameters['vs'][i]**2                
                oldD = 1./(2.*self.parameters['qs'][i])
                # update G and D parameter
                if int(self.parameters['soiltype'][i])>=0:
                    val,idx =  self.find_nearest(self.nonlinpar['nonlin strain'][int(self.parameters['soiltype'][i])-1],np.max(self.time_series[i*2]))
                    self.parameters['vs'][i]=np.sqrt((Gmax[i]*self.nonlinpar['nonlin G/Gmax'][int(self.parameters['soiltype'][i])-1][idx])/self.parameters['dn'][i])
                    self.parameters['qs'][i]=1./(2.*self.nonlinpar['nonlin damping'][int(self.parameters['soiltype'][i]-1)][idx]/100.)
                # retrieve new G and D
                newG = self.parameters['dn'][i]*self.parameters['vs'][i]**2
                newD = 1./(2.*self.parameters['qs'][i])
                # calculate rate of change as the convergence level
                Gerr[i] = np.abs(newG-oldG)/oldG
                Derr[i] = np.abs(newD-oldD)/oldD
                if verbose:
                    print '%.2f\t\t%.2f\t%.2f\t%.4f\t%.4f'%(self.parameters['hl'][i],self.parameters['vs'][i],self.parameters['qs'][i],Gerr[i],Derr[i])
            if np.max(Gerr)<=conv_level and np.max(Derr)<=conv_level:
                if verbose:
                    print('convergence has been reached! Calculation is stopped!')
                break
            
        self.lastiter = deepcopy(j)+1

        # correction of tfpair
        self.parameters['tfPair']=oldpair
        self.parameters['ntf']=len(oldpair)
        self.linear_TF2TS(parameters=self.parameters)

"""
# calculation test (circular convolution method)
print('knopoff calculation')
method = 'knopoff_sh'
TSclass00 = TSCalculator('sampleinput_linear_elastic_1layer_halfspace.dat',method,verbose=True)
"""
"""
# calculation test (Comparison between linear and linear equivalent)
fname1 = 'sampleinput_linear_elastic_1layer_halfspace.dat'
fname2 = 'sampleinput_linear_eq_elastic_1layer_halfspace.dat'
fname1 = 'ferndale2_parameter.dat'
fname2 = 'ferndale2_parameter_lineq.dat'
#print('knopoff calculation')
#method = 'knopoff_sh'
#TSclass00 = TSCalculator(fname1,method,verbose=True)
#TSclass01 = TSCalculator(fname2,method,verbose=True)
#print('kennet calculation')
#method2 = 'kennet_sh'
#TSclass02 = TSCalculator(fname1,method2,verbose=True)
#TSclass03 = TSCalculator(fname2,method2,verbose=True)
print('kramer calculation')
method2 = 'kramer286_sh'
TSclass04 = TSCalculator(fname1,method2,verbose=True)
TSclass05 = TSCalculator(fname2,method2,sublayercriteria = 2.,verbose=True)

#plot the outputs
dt = TSclass04.inp_time[1]-TSclass04.inp_time[0]
fig = plt.figure(figsize=(12,8))
a = fig.add_subplot(211)
a.plot(TSclass04.inp_time,GT.disp2acc(TSclass04.inp_signal,TSclass04.dt),label='in')
#a.plot([i*dt for i in range(len(TSclass00.time_series[0]))],TSclass00.time_series[0],label='out linear knopoff')
#a.plot([i*dt for i in range(len(TSclass01.time_series[0]))],TSclass01.time_series[0],label='out lin-eq knopoff')
#a.plot([i*dt for i in range(len(TSclass02.time_series[0]))],TSclass02.time_series[0],label='out linear kennet')
#a.plot([i*dt for i in range(len(TSclass03.time_series[0]))],TSclass03.time_series[0],label='out lin-eq kennet')
a.plot([i*dt+TSclass04.inp_time[0] for i in range(len(TSclass04.time_series[0]))],GT.disp2acc(TSclass04.time_series[0],TSclass04.dt),label='out linear kramer')
a.plot([i*dt+TSclass04.inp_time[0] for i in range(len(TSclass05.time_series[0]))],GT.disp2acc(TSclass05.time_series[0],TSclass05.dt),label='out lin-eq kramer')
a.legend(loc='upper right')
#a.set_xlim(9.5,12.5)
#a.set_xlim(5,7)
a.grid(True)

a = fig.add_subplot(212)
#a.plot(np.fft.fftfreq(len(TSclass01.inp_signal),dt)[:len(TSclass01.inp_signal)/2],np.abs(np.fft.fft(TSclass01.inp_signal))[:len(TSclass01.inp_signal)/2],label='in')
#a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass00.time_series[0])),dt)[:len(TSclass00.time_series[0])/2],np.abs(np.fft.fft(TSclass00.time_series[0]))[:len(TSclass00.time_series[0])/2],label='out lin knopoff')
#a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass01.time_series[0])),dt)[:len(TSclass01.time_series[0])/2],np.abs(np.fft.fft(TSclass01.time_series[0]))[:len(TSclass01.time_series[0])/2],label='out lin-eq knopoff')
#a.plot(TSclass00.tfclass.freq,np.abs(TSclass00.tfclass.tf[0])*0.01,label='(transfer function)/100. lin knopoff')
#a.plot(TSclass01.tfclass.freq,np.abs(TSclass01.tfclass.tf[0])*0.01,label='(transfer function)/100. lin-eq knopoff')
#a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass02.time_series[0])),dt)[:len(TSclass02.time_series[0])/2],np.abs(np.fft.fft(TSclass02.time_series[0]))[:len(TSclass02.time_series[0])/2],label='out lin kennet')
#a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass03.time_series[0])),dt)[:len(TSclass03.time_series[0])/2],np.abs(np.fft.fft(TSclass03.time_series[0]))[:len(TSclass03.time_series[0])/2],label='out lin-eq kennet')
#a.plot(TSclass02.tfclass.freq,np.abs(TSclass02.tfclass.tf[0])*0.1,label='(transfer function)/10. lin kennet')
#a.plot(TSclass03.tfclass.freq,np.abs(TSclass03.tfclass.tf[0])*0.1,label='(transfer function)/10. lin-eq kennet')
a.plot(np.fft.fftfreq(len(TSclass04.inp_signal),dt)[:len(TSclass04.inp_signal)/2],np.abs(np.fft.fft(TSclass04.inp_signal))[:len(TSclass04.inp_signal)/2],label='in')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass04.time_series[0])),dt)[:len(TSclass04.time_series[0])/2],np.abs(np.fft.fft(TSclass04.time_series[0]))[:len(TSclass04.time_series[0])/2],label='out lin knopoff')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass05.time_series[0])),dt)[:len(TSclass05.time_series[0])/2],np.abs(np.fft.fft(TSclass05.time_series[0]))[:len(TSclass05.time_series[0])/2],label='out lin-eq knopoff')
a.plot(TSclass04.tfclass.freq,np.abs(TSclass04.tfclass.tf[0])*10.*np.average(np.abs(np.fft.fft(TSclass05.time_series[0]))),label='(transfer function)/%.2e lin knopoff'%(10.*np.average(np.abs(np.fft.fft(TSclass05.time_series[0])))))
a.plot(TSclass05.tfclass.freq,np.abs(TSclass05.tfclass.tf[0])*10.*np.average(np.abs(np.fft.fft(TSclass05.time_series[0]))),label='(transfer function)/%.2e lin-eq knopoff'%(10.*np.average(np.abs(np.fft.fft(TSclass05.time_series[0])))))

a.legend(loc='best')
a.set_xscale('log')
a.set_xlim(0.1,50)
a.grid()
"""
"""
# Calculation test (sublayering sensitivity)
method2 = 'kramer286_sh'
TSclass04 = TSCalculator('sampleinput_linear_elastic_1layer_halfspace.dat',method2,verbose=False)
TSclass05 = TSCalculator('sampleinput_linear_eq_elastic_1layer_halfspace.dat',method2,sublayercriteria = 5.,numiter = 10,conv_level = 0.01,verbose=False)
TSclass06 = TSCalculator('sampleinput_linear_eq_elastic_1layer_halfspace.dat',method2,sublayercriteria = 2.,numiter = 10,conv_level = 0.01,verbose=False)

#plot the outputs
dt = TSclass04.inp_time[1]-TSclass04.inp_time[0]
fig = plt.figure(figsize=(12,8))
a = fig.add_subplot(211)
a.plot(TSclass04.inp_time,TSclass04.inp_signal,label='in')
a.plot([i*dt for i in range(len(TSclass04.time_series[0]))],TSclass04.time_series[0],label='out linear')
a.plot([i*dt for i in range(len(TSclass05.time_series[0]))],TSclass05.time_series[0],label='out lin-eq 5m')
a.plot([i*dt for i in range(len(TSclass06.time_series[0]))],TSclass06.time_series[0],label='out lin-eq 2m')

a.legend(loc='upper right')
a.set_xlim(9.5,12.5)
a.grid(True)

a = fig.add_subplot(212)
a.plot(np.fft.fftfreq(len(TSclass04.inp_signal),dt)[:len(TSclass04.inp_signal)/2],np.abs(np.fft.fft(TSclass04.inp_signal))[:len(TSclass04.inp_signal)/2],label='in')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass04.time_series[0])),dt)[:len(TSclass04.time_series[0])/2],np.abs(np.fft.fft(TSclass04.time_series[0]))[:len(TSclass04.time_series[0])/2],label='out lin kramer')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass05.time_series[0])),dt)[:len(TSclass05.time_series[0])/2],np.abs(np.fft.fft(TSclass05.time_series[0]))[:len(TSclass05.time_series[0])/2],label='out lin-eq kramer 5')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass06.time_series[0])),dt)[:len(TSclass06.time_series[0])/2],np.abs(np.fft.fft(TSclass06.time_series[0]))[:len(TSclass06.time_series[0])/2],label='out lin-eq kramer 2')
a.plot(TSclass04.tfclass.freq,np.abs(TSclass04.tfclass.tf[0])*0.1,label='(transfer function)/10. lin')
a.plot(TSclass05.tfclass.freq,np.abs(TSclass05.tfclass.tf[0])*0.1,label='(transfer function)/10. lin-eq 5m')
a.plot(TSclass06.tfclass.freq,np.abs(TSclass06.tfclass.tf[0])*0.1,label='(transfer function)/10. lin-eq 2m')

a.legend(loc='best')
a.set_xscale('log')
a.set_xlim(0.1,50)
a.grid()

fig = plt.figure(figsize=(4,12))
a = fig.add_subplot(111)

a.plot(TSclass04.parameters['vs'],np.cumsum(TSclass04.parameters['hl']),label='linear')
a.plot(TSclass05.parameters['vs'],np.cumsum(TSclass05.parameters['hl']),label='5m sublayer')
a.plot(TSclass06.parameters['vs'],np.cumsum(TSclass06.parameters['hl']),label='2m sublayer')
a.set_ylim(0,30)
a.set_xlim(350,410)
a.set_xlabel('Vs (m/s)')
a.set_ylabel('depth (m)')
a.invert_yaxis()
a.legend(loc='best')
a.grid()
"""
"""
# Calculation test (borehole outcrop correction test for kennet)
print('knopoff calculation')
method = 'knopoff_sh'
TSclass00 = TSCalculator('sampleinput_linear_elastic_1layer_halfspace.dat',method,verbose=True)
print('kennet calculation')
method2 = 'kennet_sh'
TSclass02 = TSCalculator('sampleinput_linear_elastic_1layer_halfspace.dat',method2,verbose=True)

#plot the outputs
dt = TSclass00.inp_time[1]-TSclass00.inp_time[0]
fig = plt.figure(figsize=(12,8))
a = fig.add_subplot(211)
a.plot(TSclass00.inp_time,TSclass00.inp_signal,label='in')
a.plot([i*dt for i in range(len(TSclass00.time_series[0]))],TSclass00.time_series[0],label='out linear knopoff 1')
a.plot([i*dt for i in range(len(TSclass02.time_series[0]))],TSclass02.time_series[0],label='out linear kennet 1')
a.legend(loc='upper right')
#a.set_xlim(9.5,12.5)
a.set_xlim(4,15)
a.grid(True)

a = fig.add_subplot(212)
a.plot(np.fft.fftfreq(len(TSclass00.inp_signal),dt)[:len(TSclass00.inp_signal)/2],np.abs(np.fft.fft(TSclass00.inp_signal))[:len(TSclass00.inp_signal)/2],label='in')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass00.time_series[0])),dt)[:len(TSclass00.time_series[0])/2],np.abs(np.fft.fft(TSclass00.time_series[0]))[:len(TSclass00.time_series[0])/2],label='out lin knopoff')
a.plot(TSclass00.tfclass.freq,np.abs(TSclass00.tfclass.tf[0])*0.1,label='(transfer function)/10. lin knopoff')
a.plot(np.fft.fftfreq(len(np.fft.fft(TSclass02.time_series[0])),dt)[:len(TSclass02.time_series[0])/2],np.abs(np.fft.fft(TSclass02.time_series[0]))[:len(TSclass02.time_series[0])/2],label='out lin kennet')
a.plot(TSclass02.tfclass.freq,np.abs(TSclass02.tfclass.tf[0])*0.1,label='(transfer function)/10. lin kennet')
a.legend(loc='best')
a.set_xscale('log')
a.set_xlim(0.1,50)
a.grid()
"""