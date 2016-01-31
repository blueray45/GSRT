# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 07:33:50 2016

@author: irnakat
"""
import numpy as np
import pylab as plt

fname = 'log_2016-01-31_07:25:23.105485.txt'

nlayer = np.logspace(0,2,11)
nlayer = [int(i) for i in nlayer]

tm_kramer = [2.50e-01,3.49e-02,5.64e-02,3.89e-01,4.57e-01,5.37e-01,6.31e-01,8.26e-01,1.10e+00,1.56e+00,4.06e+00]
tm_knopoff_sh = [1.07e-01,2.92e-02,3.38e-02,5.78e-02,4.90e-01,7.50e-01,1.13e+00,1.88e+00,3.23e+00,6.78e+00,2.06e+01]
tm_knopoff_sh2 = [4.09e-02,3.93e-02,4.89e-02,6.43e-02,9.55e-01,1.33e+00,1.81e+00,3.13e+00,5.35e+00,9.98e+00,2.93e+01]
tm_kennet_sh = [8.27e-02,8.20e-02,9.73e-02,1.08e-01,4.88e-01,5.54e-01,6.19e-01,8.11e-01,1.02e+00,1.41e+00,2.19e+00]
tm_knopoff_psv = [7.99e-02,7.84e-02,9.55e-01,1.24e+00,2.27e+00,3.58e+00,5.17e+00,8.99e+00,1.76e+01,4.65e+01,1.37e+02]

fig = plt.figure(figsize=(14,6))
a = fig.add_subplot(111)
a.plot(nlayer,tm_kramer,'o-',label='kramer sh')
a.plot(nlayer,tm_knopoff_sh,'o-',label='knopoff sh')
a.plot(nlayer,tm_knopoff_sh2,'o-',label='knopoff sh2')
a.plot(nlayer,tm_kennet_sh,'o-',label='kennet sh')
#a.plot(nlayer,tm_knopoff_psv,'o-',label='knopoff psv')
a.grid(which='both')
a.legend(loc='best',fancybox=True,framealpha=0.5)
a.set_yscale('log')
a.set_xscale('log')
a.set_ylabel('time')
a.set_xlabel('number of layers')

N = 1
bar_kramer = [3.79e-02]
bar_knopoff_sh = [3.02e-02]
bar_knopoff_sh2= [4.05e-02]
bar_kennet_sh = [8.33e-02]
bar_knopoff_psv=[2.55e+00]

width =0.05
ind = np.arange(N)
fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(111)
ax.grid()
rects1 = ax.bar(ind, bar_kramer, width,color='b')
rects2 = ax.bar(ind+width, bar_knopoff_sh, width,color='g')
rects3 = ax.bar(ind+2.*width, bar_knopoff_sh2, width,color='r')
rects4 = ax.bar(ind+3.*width, bar_kennet_sh, width,color='cyan')
#rects5 = ax.bar(ind+4.*width, bar_knopoff_psv, width,color='k')

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%.3f' %height,
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
#autolabel(rects5)

ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Kramer', 'Knopoff SH','Knopoff SH2','Kennet SH'),loc='best')

N = 1
bar_kramer = [4.73e-01]
bar_knopoff_sh = [1.01e-01]
bar_knopoff_sh2= [1.42e-01]
bar_kennet_sh = [2.82e-01]
bar_knopoff_psv=[7.89e+00]

width =0.05
ind = np.arange(N)
fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(111)
ax.grid()
rects1 = ax.bar(ind, bar_kramer, width,color='b')
rects2 = ax.bar(ind+width, bar_knopoff_sh, width,color='g')
rects3 = ax.bar(ind+2.*width, bar_knopoff_sh2, width,color='r')
rects4 = ax.bar(ind+3.*width, bar_kennet_sh, width,color='cyan')
#rects5 = ax.bar(ind+4.*width, bar_knopoff_psv, width,color='k')

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%.3f' %height,
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
#autolabel(rects5)

ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Kramer', 'Knopoff SH','Knopoff SH2','Kennet SH'),loc='best')