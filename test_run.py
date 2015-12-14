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
"""
print 'creating class'
theclass = TFC(data)
print 'class creation was succeed! proceed!'

print 'TF calculatoin using kramer approach'
tf1 = theclass.tf_kramer286_sh() # check/verufy kramer calculation
print 'calculation was succeed!'

ax.plot(theclass.freq,np.abs(tf1[0]),'b',label='kramer286')

print 'TF calculation using simple knopoff approach'
tf2 = theclass.tf_knopoff_sh()
print 'calculation has been finished!'

ax.plot(theclass.freq,np.abs(tf2[0]),'r',label='knopoff sh')
"""
"""
print 'TF calculation using complete knopoff sh approach'
fname3 = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
data3 = IOfile.parsing_input_file(fname3)
theclass3 = TFC(data3)
tf3 = theclass3.tf_knopoff_sh_adv()
print 'Done! Amazing!'

ax.plot(theclass3.freq,np.abs(tf3[0]),'g--',lw=4.,label='knopoff sh complete')
"""
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
data2['iang'] = np.deg2rad(85.)
print 'creating class'
theclass2 = TFC(data2)
print 'class creation was succeed! proceed!'
print 'TF calculation using complete knopoff psv approach'
tf4 = theclass2.tf_knopoff_psv_adv()
print 'Done! Amazing!'
ax.plot(theclass2.freq,np.abs(tf4[0]),'r-.',lw=2.,label='htf knopoff psv-s complete')
ax.plot(theclass2.freq,np.abs(tf4[1]),'g-.',lw=2.,label='vtf knopoff psv-s complete')
#ax.plot(theclass2.freq,np.sqrt(np.abs(tf4[0])**2+np.abs(tf4[0])**2),'b-.',lw=2.,label='amp knopoff psv-s complete')
"""
fname5 = 'sampleinput_psv_p_linear_elastic_1layer_halfspace.dat'
data5 = IOfile.parsing_input_file(fname5)
theclass5 = TFC(data5)
htf5,vtf5 = theclass5.tf_kennett()
plt.plot(theclass5.freq,np.abs(htf5),'r-.',lw=2.,label='kennet sh')
"""

plt.xlabel('frequency (Hz)')
plt.ylabel('Amplification')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1.,50.)
plt.ylim(0.5,5.)
plt.grid(True,which='both')
plt.legend(loc='best',fancybox=True,framealpha=0.5)