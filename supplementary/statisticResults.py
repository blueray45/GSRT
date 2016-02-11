# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 07:33:50 2016

@author: irnakat
"""
import os
import numpy as np
import pylab as plt

dirname = '/media/data/gdriveubuntu/Kuliah/Semester III/04-Individual_Project/Site_Response_Toolkit/02-Code/BatchTest/logcollections'
fname = 'log_2016-02-08_14:32:20.571967.txt'

TfBlockModel = []
TfBlockElapsedRead = []
TfBlockStdRead = []
TfBlockElapsed = []
TfBlockStd = []

TsBlockModel = []
TsBlockIter = []
TsBlockElapsed = []
TsBlockStd = []

with open(os.path.join(dirname,fname)) as f:
    counter = -1
    for ln in f.readlines():
        if not ln.rstrip('\n'):
            # ignore empty rows
            continue
        
        if ' - ' in ln:
            counter +=1
            internalCounter = 0
#            print('reading %03d'%counter)
            continue
        
        if counter > 2:
            # skipping header
        
            if internalCounter==0 :
                # skipping header row
                internalCounter+=1
                continue
            
            if counter<=7 and internalCounter>=0:
                # read first data block
                tmp = ln.split()
                TfBlockModel.append(int(tmp[0]))
                TfBlockElapsedRead.append(float(tmp[1]))
                TfBlockStdRead.append(float(tmp[2]))
                TfBlockElapsed.append(float(tmp[3]))
                TfBlockStd.append(float(tmp[4]))
                internalCounter+=1
                continue
                
            if counter>7 and counter<=12 and internalCounter>=0:
                # read second data block
                tmp = ln.split()
                TsBlockModel.append(int(tmp[0]))
                TsBlockIter.append(float(tmp[1]))
#                print TsBlockIter[-1],float(tmp[2]),float(tmp[2])/TsBlockIter[-1]
                TsBlockElapsed.append(float(tmp[2])/TsBlockIter[-1])
                TsBlockStd.append(float(tmp[3]))
                internalCounter+=1
                continue
            
            

#print TsBlockElapsed
#print TsBlockStd

# TF

N = 3
bar_kramer = [0.0962, 0.39, 0.57]
bar_knopoff_sh = [0.0506, 0.112, 0.153]
bar_knopoff_sh2= [ 0.0674, 0.189, 0.259]
bar_kennet_sh = [0.143, 0.345, 0.466]
#bar_knopoff_psv=[0.176, 0.572]

bar_kramer_std = [0.000511, 0.00105, 0.00112]
bar_knopoff_sh_std = [0.000604,0.000377, 0.00209]
bar_knopoff_sh2_std= [0.000213, 0.00187, 0.000424]
bar_kennet_sh_std = [0.000122, 0.000376, 0.00139]
#bar_knopoff_psv=[0.000616, 0.00226]

width =0.2
ind = np.arange(N)
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111)
ax.grid()
rects1 = ax.bar(ind, bar_kramer, width,color='b',yerr=bar_kramer_std)
rects2 = ax.bar(ind+width, bar_knopoff_sh, width,color='g',yerr=bar_knopoff_sh_std)
rects3 = ax.bar(ind+2.*width, bar_knopoff_sh2, width,color='r',yerr=bar_knopoff_sh2_std)
rects4 = ax.bar(ind+3.*width, bar_kennet_sh, width,color='cyan',yerr=bar_kennet_sh_std)
#rects5 = ax.bar(ind+4.*width, bar_knopoff_psv, width,color='k')

ax.set_xticks(ind + 2*width)
ax.set_xticklabels(('1 layer', '6 layers', 'ferndale'))

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

ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Kramer', \
          'Knopoff SH','Knopoff SH2','Kennet SH'),loc='upper left')
          
# TS

N = 6
bar_kramer = [0.419, 1.57, 4.3, 3.05, 3.46, 26.75]
bar_knopoff_sh = [0.239, 0.482, 4.51,3.075, 1.02, 107.5]
bar_knopoff_sh2= [0.306, 0.773, 7.140000000000001, 4.95,1.65, 190.0]
bar_kennet_sh = [0.607, 1.39, 12.7, 8.9, 2.85, 309.5]
#bar_knopoff_psv=[0.739,2.26]

bar_kramer_std = [0.000825, 0.00152, 0.679, 0.0144, 0.0159, 0.0546]
bar_knopoff_sh_std = [0.00116, 0.00142,0.0472, 0.0123, 0.00228, 1.1]
bar_knopoff_sh2_std= [ 0.00144, 0.000921, 0.0519, 0.0298, 0.00205,0.438]
bar_kennet_sh_std = [0.000272, 0.00191, 0.0676, 0.0138, 0.00361, 0.641]
#bar_knopoff_psv=[0.0135, 0.00371]

# linear
N = 3
bar_kramer = [0.419, 1.57, 3.46]
bar_knopoff_sh = [0.239, 0.482, 1.02]
bar_knopoff_sh2= [0.306, 0.773,1.65]
bar_kennet_sh = [0.607, 1.39, 2.85]
#bar_knopoff_psv=[0.739,2.26]

bar_kramer_std = [0.000825, 0.00152, 0.0159]
bar_knopoff_sh_std = [0.00116, 0.00142, 0.00228]
bar_knopoff_sh2_std= [ 0.00144, 0.000921, 0.00205]
bar_kennet_sh_std = [0.000272, 0.00191, 0.00361]
#bar_knopoff_psv=[0.0135, 0.00371]

width =0.2
ind = np.arange(N)
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111)
ax.grid()
rects1 = ax.bar(ind, bar_kramer, width,color='b',yerr=bar_kramer_std)
rects2 = ax.bar(ind+width, bar_knopoff_sh, width,color='g',yerr=bar_knopoff_sh_std)
rects3 = ax.bar(ind+2.*width, bar_knopoff_sh2, width,color='r',yerr=bar_knopoff_sh2_std)
rects4 = ax.bar(ind+3.*width, bar_kennet_sh, width,color='cyan',yerr=bar_kennet_sh_std)
#rects5 = ax.bar(ind+4.*width, bar_knopoff_psv, width,color='k')

ax.set_xticks(ind + 2*width)
#ax.set_xticklabels(('1 layer', '6 layers', '1 lineq', '6 lineq', 'ferndale', 'fern_lineq'))
ax.set_xticklabels(('1 layer', '6 layers', 'ferndale'))

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

ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Kramer', \
          'Knopoff SH','Knopoff SH2','Kennet SH'),loc='upper left')
#ax.set_ylim(0,4)

# TS

N = 6
bar_kramer = [0.419, 1.57, 4.3, 3.05, 3.46, 26.75]
bar_knopoff_sh = [0.239, 0.482, 4.51,3.075, 1.02, 107.5]
bar_knopoff_sh2= [0.306, 0.773, 7.140000000000001, 4.95,1.65, 190.0]
bar_kennet_sh = [0.607, 1.39, 12.7, 8.9, 2.85, 309.5]
#bar_knopoff_psv=[0.739,2.26]

bar_kramer_std = [0.000825, 0.00152, 0.679, 0.0144, 0.0159, 0.0546]
bar_knopoff_sh_std = [0.00116, 0.00142,0.0472, 0.0123, 0.00228, 1.1]
bar_knopoff_sh2_std= [ 0.00144, 0.000921, 0.0519, 0.0298, 0.00205,0.438]
bar_kennet_sh_std = [0.000272, 0.00191, 0.0676, 0.0138, 0.00361, 0.641]
#bar_knopoff_psv=[0.0135, 0.00371]

# linear
N = 3
bar_kramer = [0.419, 1.57, 3.46]
bar_knopoff_sh = [0.239, 0.482, 1.02]
bar_knopoff_sh2= [0.306, 0.773,1.65]
bar_kennet_sh = [0.607, 1.39, 2.85]
#bar_knopoff_psv=[0.739,2.26]

bar_kramer_std = [0.000825, 0.00152, 0.0159]
bar_knopoff_sh_std = [0.00116, 0.00142, 0.00228]
bar_knopoff_sh2_std= [ 0.00144, 0.000921, 0.00205]
bar_kennet_sh_std = [0.000272, 0.00191, 0.00361]
#bar_knopoff_psv=[0.0135, 0.00371]

# linear equivalent
N = 2
bar_kramer = [ 4.3, 3.05]
bar_knopoff_sh = [ 4.51,3.075]
bar_knopoff_sh2= [7.140000000000001, 4.95,]
bar_kennet_sh = [ 12.7, 8.9]
#bar_knopoff_psv=[0.739,2.26]

bar_kramer_std = [0.679, 0.0144]
bar_knopoff_sh_std = [0.0472, 0.0123]
bar_knopoff_sh2_std= [0.0519, 0.0298]
bar_kennet_sh_std = [0.0676, 0.0138]
#bar_knopoff_psv=[0.0135, 0.00371]

width =0.2
ind = np.arange(N)
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111)
ax.grid()
rects1 = ax.bar(ind, bar_kramer, width,color='b',yerr=bar_kramer_std)
rects2 = ax.bar(ind+width, bar_knopoff_sh, width,color='g',yerr=bar_knopoff_sh_std)
rects3 = ax.bar(ind+2.*width, bar_knopoff_sh2, width,color='r',yerr=bar_knopoff_sh2_std)
rects4 = ax.bar(ind+3.*width, bar_kennet_sh, width,color='cyan',yerr=bar_kennet_sh_std)
#rects5 = ax.bar(ind+4.*width, bar_knopoff_psv, width,color='k')

ax.set_xticks(ind + 2*width)
#ax.set_xticklabels(('1 layer', '6 layers', '1 lineq', '6 lineq', 'ferndale', 'fern_lineq'))
#ax.set_xticklabels(('1 layer', '6 layers', 'ferndale'))
ax.set_xticklabels(('1 lineq', '6 lineq'))

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

ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]), ('Kramer', \
          'Knopoff SH','Knopoff SH2','Kennet SH'),loc='upper left')
#ax.set_ylim(0,4)
# many layers

nlayer = np.logspace(0,2,11)
nlayer = [int(i) for i in nlayer]

tm_kramer = [9.55e-02,9.54e-02,1.55e-01,2.13e-01,3.89e-01,6.24e-01,9.19e-01,1.50e+00,2.31e+00,3.70e+00,5.83e+00]
tm_knopoff_sh = [4.97e-02,4.99e-02,6.31e-02,7.52e-02,1.11e-01,1.59e-01,2.24e-01,3.54e-01,5.53e-01,9.65e-01,1.61e+00]
tm_knopoff_sh2 = [6.73e-02,6.72e-02,9.12e-02,1.16e-01,1.86e-01,2.79e-01,4.00e-01,6.43e-01,9.96e-01,1.67e+00,2.75e+00]
tm_kennet_sh = [1.43e-01,1.43e-01,1.83e-01,2.23e-01,3.42e-01,5.02e-01,6.96e-01,1.08e+00,1.63e+00,2.57e+00,3.99e+00]
tm_knopoff_psv = [1.76e-01,1.75e-01,2.53e-01,3.30e-01,5.73e-01,8.87e-01,1.30e+00,2.12e+00,3.39e+00,6.19e+00,9.76e+00]

tm_kramer_std = [1.44e-04,2.14e-04,4.10e-04,1.66e-04,2.67e-04,1.92e-03,3.88e-03,3.47e-03,2.79e-03,5.46e-03,1.31e-02]
tm_knopoff_sh_std = [6.80e-05,4.60e-05,8.54e-04,1.07e-04,1.40e-04,1.27e-04,9.05e-04,2.69e-04,4.70e-04,1.48e-02,4.18e-03]
tm_knopoff_sh2_std = [1.42e-04,1.18e-04,8.19e-04,1.70e-04,2.44e-04,4.18e-04,6.25e-04,3.39e-03,9.72e-04,1.37e-03,4.13e-03]
tm_kennet_sh_std = [1.65e-04,9.97e-05,2.04e-04,1.03e-04,2.40e-04,3.17e-03,6.72e-04,1.17e-03,8.67e-03,8.41e-03,5.09e-03]
tm_knopoff_psv_std = [4.31e-04,3.89e-04,1.10e-03,8.79e-04,9.91e-03,2.31e-03,4.71e-03,5.91e-03,1.45e-02,1.61e-02,2.39e-02]

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


# memory leak plot
nlayer = [1,2,5,10,20,30,50,70,80,90,100]

memory_kramer = [0.203,0.988,1.535,2.406,0.383,0.371,0.449,0.379,0.379, 0.449,0.402,]
memory_knopoff_sh = [2.496,2.434,2.156,2.469,3.012,3.574,3.770,5.941,6.461,4.480,4.965,]
memory_knopoff_sh_adv = [0.184,0.137,0.004,0.000,0.242,0.309,0.668,0.926,1.211,1.527,0.629,]
memory_knopoff_psv = [0.316,0.324,0.707,0.270,0.961,1.258,2.418,3.770,6.227,5.410,8.074,]
memory_kennet = [0.121,0.125,0.121,0.,0.012,0.,0.0,0.,0.,0.,0.004]

fig3 = plt.figure(dpi=300)
ax1 = fig3.add_subplot(111)
ax1.plot(nlayer,memory_kramer,label='kramer')
ax1.plot(nlayer,memory_knopoff_sh,label='knopoff_sh')
ax1.plot(nlayer,memory_knopoff_sh_adv,label='knopoff_sh2')
ax1.plot(nlayer,memory_knopoff_psv,label='knopoff_psv')
ax1.plot(nlayer,memory_kennet,label='kennet')
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlabel('number of layer')
ax1.set_ylabel('memory leak (Mb)')


# memory usage plot
fname = ['mprofile_20160211080207_001layer.dat',
         'mprofile_20160211074931_020layer.dat',
         'mprofile_20160211075542_050layer.dat',
         'mprofile_20160211075938_070layer.dat',
         'mprofile_20160211085440_090layer.dat',
         'mprofile_20160211085222_100layer.dat',
         'mprofile_20160211090325_200layer.dat',]
fname = [os.path.pardir+os.sep+i for i in fname]

memlist = []
timelist = []
relativeTimeList = []
fig3 = plt.figure(figsize=(10,8),dpi=300)
ax1 = fig3.add_subplot(111)

nlayer = []

memory_kramer = []
memory_knopoff_sh = []
memory_knopoff_sh_adv = []
memory_knopoff_psv = []
memory_kennet = []

for i in fname:
    label = i[i.index('layer')-3:i.index('layer')]+'layers'
    with open(i,'rb') as f:
        nlayer.append(int(i[i.index('layer')-3:i.index('layer')]))
        memory1 = []
        time1 = []
        memory = []
        time = []
        for ln in f.readlines():
            if 'MEM' in ln:
                tmp = ln.split(' ')
                memory1.append(float(tmp[1]))
                time1.append(float(tmp[2]))
            if 'FUNC __main__.calc_kramer' in ln:
                tmp = ln.split(' ')
                starttime = float(tmp[3])
                memory_kramer.append(float(tmp[4])-float(tmp[2]))
            if 'FUNC __main__.calc_kennet' in ln:
                tmp = ln.split(' ')
                endtime = float(tmp[5])
                memory_kennet.append(float(tmp[4])-float(tmp[2]))
            if 'FUNC __main__.calc_knopoff_sh_adv' in ln:
                tmp = ln.split(' ')
                memory_knopoff_sh_adv.append(float(tmp[4])-float(tmp[2]))
            if 'FUNC __main__.calc_knopoff_sh ' in ln:
                tmp = ln.split(' ')
                memory_knopoff_sh.append(float(tmp[4])-float(tmp[2]))                
            if 'FUNC __main__.calc_knopoff_psv' in ln:
                tmp = ln.split(' ')
                memory_knopoff_psv.append(float(tmp[4])-float(tmp[2]))                
        for i in range(len(time1)):
            if time1[i]>=starttime and time1[i]<=endtime:
                time.append(time1[i])
                memory.append(memory1[i])
        timelist.append(np.array(time)-time[0])
        memlist.append(np.array(memory))
        relativeTimeList.append(np.linspace(0.,1.,len(time)))
    ax1.plot(relativeTimeList[-1],memlist[-1],label=label)


ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlabel('normalized time')
ax1.set_ylabel('memory consumption (Mb)')
ax1.set_ylim(70,110)


fig3 = plt.figure(dpi=300)
ax1 = fig3.add_subplot(111)
ax1.plot(nlayer,memory_kramer,label='kramer')
ax1.plot(nlayer,memory_knopoff_sh,label='knopoff_sh')
ax1.plot(nlayer,memory_knopoff_sh_adv,label='knopoff_sh2')
ax1.plot(nlayer,memory_knopoff_psv,label='knopoff_psv')
ax1.plot(nlayer,memory_kennet,label='kennet')
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlabel('number of layer')
ax1.set_ylabel('memory usage (Mb)')
ax1.set_xlim(1,1000)
ax1.set_ylim(1,250)
ax1.set_xscale('log')
ax1.set_yscale('log')

# memory usage theory

nlayer = np.logspace(0,3,16)
kramer = [(24*i*0.)/(1024.*1024.) for i in nlayer]
knopoff_sh = [24*(2*i)**2/(1024.*1024.) for i in nlayer]
knopoff_sh2 = [24*(2*i)**2/(1024.*1024.) for i in nlayer]
knopoff_psv = [24*(4*i)**2/(1024.*1024.) for i in nlayer]
kennet = [24*60*i/(1024.*1024.) for i in nlayer]

fig3 = plt.figure(dpi=300)
ax1 = fig3.add_subplot(111)
ax1.plot(nlayer,kramer,label='kramer')
ax1.plot(nlayer,knopoff_sh,label='knopoff_sh')
ax1.plot(nlayer,knopoff_sh2,'-.',label='knopoff_sh2')
ax1.plot(nlayer,knopoff_psv,label='knopoff_psv')
ax1.plot(nlayer,kennet,label='kennet')
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlabel('number of layer')
ax1.set_ylabel('theoretical memory usage (Mb)')
ax1.set_xlim(1,1000)
ax1.set_ylim(1,250)
ax1.set_xscale('log')
ax1.set_yscale('log')