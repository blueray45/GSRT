# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:06:03 2016

@author: irnakat
"""

class TSCalculator:
    """
    class which manages calculation of TimeSeries (TS) given the Transfer Function (TF)
    """
    import IOfile
    from TFCalculator import TFCalculator as TFC
    import numpy as np
    def __init__(self,tf,inputsignal,parfile):
        self.htf = tf[0]
        self.vtf = tf[1]
        self.inputsignal = inputsignal
        self.inp_time,self.inp_signal = IOfile.read_ascii_seismogram(self.inputsignal)
        self.parameters = IOfile.parsing_input_file(parfile)
    
    def TF2TS(self):
        self.inp_ft = np.fft.fft(self.inp_signal)
        self.inp_freq = np.fft.fftfreq(len(self.inp_signal),d=self.inp_time[1]-self.inp_time[0])
        # using knopoff SH for the time being <-- need to be changed whenever possible to dynamic choice
        TFclass = TFC(datash,self.inp_freq)
        TFclass.tf_knopoff_sh() #<-- need to be changed in the near future
        case = 'sh'        
        # convolution in frequency domain
        if case=='sh':
            self.out_hft = TFclass.ft[0]*self.inp_ft
        else:
            self.out_hft = TFclass.ft[0]*self.inp_ft
            self.out_vft = TFclass.ft[1]*self.inp_ft
        self.out_hamp = np.fft.ifft(self.out_hft)
        return self.inp_time,self.out_hamp