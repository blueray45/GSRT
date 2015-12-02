# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:55:33 2015

@author: irnakat
"""

# exporting external helper module
import numpy as np
import pylab as plt
import commondata as cd

# test IOfile
import IOfile
from TFCalculator import TFCalculator as TFC

fname = 'sampleinput_linear_elastic_6layer_halfspace.dat'
fname = 'sampleinput_psv_p_linear_elastic_1layer_halfspace.dat'
mode1 = cd.mode[0]
mode2 = cd.mode[1]

data = IOfile.parsing_input_file(fname)
f = plt.figure(figsize=(10.,10.))
ax = f.add_subplot(111)
print data
print('data reading was OK! proceed!')

print 'creating class'
theclass = TFC(data)
print 'class creation was succeed! proceed!'

print 'TF calculatoin using kramer approach'
#tf1 = theclass.tf_kramer286_sh() # check/verufy kramer calculation
print 'calculation was succeed!'

#ax.plot(theclass.freq,np.abs(tf1[0]),'b',label='kramer286')

print 'TF calculation using simple knopoff approach'
tf2 = theclass.tf_knopoff_sh()
print 'calculation has been finished!'

ax.plot(theclass.freq,np.abs(tf2[0]),'r',label='knopoff sh')

print 'TF calculation using complete knopoff sh approach'
fname3 = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
data3 = IOfile.parsing_input_file(fname3)
theclass3 = TFC(data3)
tf3 = theclass3.tf_knopoff_sh_adv()
print 'Done! Amazing!'

ax.plot(theclass.freq,np.abs(tf3[0]),'g--',lw=4.,label='knopoff sh complete')

fname2 = 'sampleinput_psv_p_linear_elastic_1layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
print 'creating class'
theclass2 = TFC(data2)
print 'class creation was succeed! proceed!'
print 'TF calculation using complete knopoff psv approach'
htf4,vtf4 = theclass.tf_knopoff_psv_adv()
print 'Done! Amazing!'
plt.plot(theclass2.freq,np.abs(htf4[0]),'r-.',lw=2.,label='knopoff psv-s complete')


plt.xlabel('frequency (Hz)')
plt.ylabel('Amplification')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1.,50.)
plt.ylim(0.8,10.)
plt.grid(True,which='both')
plt.legend(loc='best',fancybox=True,framealpha=0.5)