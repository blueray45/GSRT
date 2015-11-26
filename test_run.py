# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:55:33 2015

@author: irnakat
"""

# exporting external helper module
import numpy as np
import pylab as plt

# test IOfile
import IOfile
from TFCalculator import TFCalculator as TFC

fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
mode1 = 'linear-elastic'
mode2 = 'linear-viscoelastic'

data = IOfile.parsing_input_file(fname,mode2)

print data
print('data reading was OK! proceed!')

print 'creating class'
theclass = TFC(data)
print 'class creation was succeed! proceed!'

print 'TF calculatoin using kramer approach'
tf1 = theclass.tf_kramer286_sh() # check/verufy kramer calculation
print 'calculation was succeed!'

plt.plot(theclass.freq,np.abs(tf1[0]),'b',label='kramer286')

print 'TF calculation using simple knopoff approach'
tf2 = theclass.tf_knopoff_sh()
print 'calculation has been finished!'

plt.plot(theclass.freq,np.abs(tf2[0]),'r--',label='knopoff sh')







plt.xlabel('frequency (Hz)')
plt.ylabel('Amplification')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1.,50.)
plt.ylim(0.8,10.)
plt.legend()