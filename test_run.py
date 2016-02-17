# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:55:33 2015

@author: irnakat
"""

# exporting external helper module
import numpy as np
import pylab as plt
import commondata as cd
import os

# test IOfile
import IOfile
from TFCalculator import TFCalculator as TFC
import TFDisplayTools
from TSCalculator import TSCalculator as TSC

# single layer test case

# filename
basedir = 'Example/Input'
fname = os.path.join(basedir,'sampleinput_linear_elastic_1layer_halfspace.dat')
fname2 =  os.path.join(basedir,'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat')
fname3 =  os.path.join(basedir,'sampleinput_psv_p_linear_elastic_1layer_halfspace.dat')
fname4 =  os.path.join(basedir,'GoverGmax.dat')

# input file reading
datash = IOfile.parsing_input_file(fname)
datapsvs = IOfile.parsing_input_file(fname2)
# datapsvp = IOfile.parsing_input_file(fname3)
datanonlin = IOfile.parsing_nonlinear_parameter(fname4,True)

def modnlayer(data,nlayer):
    nl = nlayer
    data['sourceloc']=nl
    data['nlayer']=nl+1
    data['tfPair'][0][1]=nl
    newhl = []; newvs = []; newdn = []; newqs = []; newvp = []; newqp = []
    newsoiltype = []
    for j in range(nl):
        newhl.append(data['hl'][0]/nl)
        newvs.append(data['vs'][0])
        newqs.append(data['qs'][0])
        newdn.append(data['dn'][0])
    newhl.append(data['hl'][-1])
    newvs.append(data['vs'][-1])
    newqs.append(data['qs'][-1])
    newdn.append(data['dn'][-1])
    data['hl']=newhl
    data['vs']=newvs
    data['qs']=newqs
    data['dn']=newdn
    try:
        for j in range(nl):
            newvp.append(data['vp'][0])
            newqp.append(data['qp'][0])
        newvp.append(data['vp'][-1])
        newqp.append(data['qp'][-1])
        data['vp']=newvp
        data['qp']=newqp
    except:
        pass
    try:
        for j in range(nl):
            newsoiltype.append(data['soiltype'][0])
        newsoiltype.append(data['soiltype'][-1])
        data['soiltype']=newsoiltype
    except:
        pass
    return data

nlayer = 1
datash = modnlayer(datash,nlayer)
datapsvs = modnlayer(datapsvs,nlayer)

#@profile
def calc_kramer(data):
    theclass = TFC(data)
    theclass.tf_kramer286_sh()
    return theclass

#@profile
def calc_knopoff_sh(data):
    theclass = TFC(data)
    theclass.tf_knopoff_sh()
    return theclass

#@profile
def calc_knopoff_sh_adv(data):
    theclass = TFC(data)
    theclass.tf_knopoff_sh_adv()
    return theclass

#@profile
def calc_knopoff_psv(data):
    theclass = TFC(data)
    theclass.tf_knopoff_psv_adv()
    return theclass

#@profile
def calc_kennet(data):
    theclass = TFC(data)
    theclass.tf_kennet_sh()
    return theclass

# kramer
print 'TF calculatoin using kramer approach'
theclass1 = calc_kramer(datash)
print 'calculation has been finished!'

# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass2 = calc_knopoff_sh(datash)
print 'calculation has been finished!'

# knopoff sh complete
print 'TF calculation using complete knopoff sh approach'
theclass3 = calc_knopoff_sh_adv(datash)
print 'calculation has been finished!'

# knopoff psv-s
print 'TF calculation using complete knopoff psv-s approach'
theclass4 = calc_knopoff_psv(datapsvs)
print 'calculation has been finished!'

# kennet sh
print 'TF calculation using kennet sh method'
theclass5 = calc_kennet(datash)
print 'calculation has been finished!'


TFDisplayTools.TFPlot(theclass1,theclass2,theclass3,theclass5,theclass4, \
    label=['Kramer SH','Knopoff SH','Knopoff SH2','Kennet SH','Knopoff PSV'])
    
TFDisplayTools.PhasePlot(theclass1,theclass2,theclass3,theclass5,theclass4, \
    label=['Kramer SH','Knopoff SH','Knopoff SH2','Kennet SH','Knopoff PSV'])


"""
from TSCalculator import TSCalculator as TSC
TSclass = TSC('ricker.dat',fname)
time, amp =TSclass.TF2TS()
"""
"""
# creating ricker signal
from scipy import signal
points = 4000
a = 10.0
inputmotion = signal.ricker(points,a)/1000.
inputtime = np.linspace(0.,20.,points)
with open('ricker.dat','w+') as f:
    for i in range(len(inputtime)):
        f.write('%2.5f %2.5f\n'%(inputtime[i],inputmotion[i]))



# creating dirac signal
points = 4000
inputmotion = np.zeros(points)
inputmotion[512] = 0.001
inputtime = np.linspace(0.,20.,points)
with open('dirac.dat','w+') as f:
    for i in range(len(inputtime)):
        f.write('%2.5f %2.5f\n'%(inputtime[i],inputmotion[i]))


"""
"""

# creating dirac signal for paolucci
points = 2000
inputmotion = np.zeros(points)
inputmotion[512] = 1.
inputtime = np.linspace(0.,20.,points)
with open('diracpaolucci.dat','w+') as f:
    for i in range(len(inputtime)/8):
        j=i*8
        f.write('%2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f \n'%(inputmotion[j],inputmotion[j+1],inputmotion[j+2],inputmotion[j+3],
                                                                      inputmotion[j+4],inputmotion[j+5],inputmotion[j+6],inputmotion[j+7]))

"""

"""
# creating ricker signal for paolucci
from scipy import signal
points = 2000
a = 100.0
inputmotion = signal.ricker(points,a)
inputtime = np.linspace(0.,20.,points)
with open('rickerpaolucci.dat','w+') as f:
    for i in range(len(inputtime)/8):
        j=i*8
        f.write('%2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f \n'%(inputmotion[j],inputmotion[j+1],inputmotion[j+2],inputmotion[j+3],
                                                                      inputmotion[j+4],inputmotion[j+5],inputmotion[j+6],inputmotion[j+7]))
"""

"""
# read output from paolucci and plot it
from TSCalculator import TSCalculator
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
datash = IOfile.parsing_input_file(fname)
theclass2 = TFC(datash)
theclass2.tf_knopoff_sh()
TSclass01 = TSCalculator(theclass2)
TSclass01.linear_TF2TS()

points = 4000
inputmotion = np.zeros(points/2)
inputmotion[512] = 1.
data = []
with open('psvq.out','r') as f:
    next(f)
    for line in f:
        tmp2= line.split()
        for i in range(8):
            data.append(tmp2[i])
fig = plt.figure(figsize=(12,12))
a = fig.add_subplot(311)
a.plot(np.linspace(0.,20.,points),data,label='paolucci code')
a.plot(np.linspace(0.,20.,points),TSclass01.time_series[0],lw=2,label='knopoff SH')
a.plot(np.linspace(0.,10.,points/2),inputmotion,label='input motion')
a.grid(True)
a.legend(loc='best',fancybox=True,framealpha=0.5)
a.set_xlim(1,4)
a = fig.add_subplot(312)
tmp = np.fft.fftfreq(points,d=0.01)
a.plot(tmp[:2000],np.abs(np.fft.fft(data))[:2000])
a.plot(theclass2.freq,np.abs(theclass2.tf[0]))
a.set_xlim(0,50)
a.grid(True)
a = fig.add_subplot(313)
tmp = np.fft.fftfreq(points,d=0.01)
a.plot(tmp[:2000],np.angle(np.fft.fft(data))[:2000])
a.plot(theclass2.freq,np.angle(theclass2.tf[0]))
a.set_xlim(0,50)
a.grid(True)
"""

"""
# additional calculations
fname3 = 'sampleinput_linear_elastic_2layer_halfspace.dat'
datash2 = IOfile.parsing_input_file(fname3)
theclass6 = TFC(datash2)
theclass6.tf_knopoff_sh_adv()


theclass7 = TFC(datash2)
theclass7.tf_kennet_sh()

fname4 = 'sampleinput_linear_elastic_6layer_halfspace_400vs30.dat'
datash3 = IOfile.parsing_input_file(fname4)
theclass8 = TFC(datash3)
theclass8.tf_knopoff_sh_adv()
theclass9 = TFC(datash3)
theclass9.tf_kennet_sh()


TFDisplayTools.TFPlot(theclass3,theclass5, \
    label=['1 layer - Knopoff SH','1 layer - Kennet SH'])
    
TFDisplayTools.PhasePlot(theclass3,theclass5, \
    label=['1 layer - Knopoff SH','1 layer - Kennet SH'])
    
TFDisplayTools.TFPlot(theclass6,theclass7, \
    label=['2 layers- Knopoff SH','2 layers- Kennet SH'])
    
TFDisplayTools.PhasePlot(theclass6,theclass7, \
    label=['2 layers- Knopoff SH','2 layers- Kennet SH'])
    
TFDisplayTools.TFPlot(theclass8,theclass9, \
    label=['6 layers- Knopoff SH','6 layers- Kennet SH'])
    
TFDisplayTools.PhasePlot(theclass8,theclass9, \
    label=['6 layers- Knopoff SH','6 layers- Kennet SH'])

TFDisplayTools.TFPlot(theclass3,theclass5,theclass6,theclass7,theclass8,theclass9, \
    label=['1 layer - Knopoff SH','1 layer - Kennet SH','2 layers- Knopoff SH','2 layers- Kennet SH','6 layers- Knopoff SH','6 layers- Kennet SH'])
    
TFDisplayTools.PhasePlot(theclass3,theclass5,theclass6,theclass7,theclass8,theclass9, \
    label=['1 layer - Knopoff SH','1 layer - Kennet SH','2 layers- Knopoff SH','2 layers- Kennet SH','6 layers- Knopoff SH','6 layers- Kennet SH'])

TFDisplayTools.TFPlot(theclass1,theclass2,theclass3,theclass4,theclass5,theclass10, \
    label=['Kramer SH','Knopoff Simple SH','Knopoff Complete SH','Knopoff PSV','Kennet SH','Kennet PSV'])
    
TFDisplayTools.PhasePlot(theclass1,theclass2,theclass3,theclass4,theclass5,theclass10, \
    label=['Kramer SH','Knopoff Simple SH','Knopoff Complete SH','Knopoff PSV','Kennet SH','Kennet PSV'])

TFDisplayTools.TFPlot(theclass4,theclass10, \
    label=['Knopoff PSV','Kennet PSV'])
    
TFDisplayTools.PhasePlot(theclass4,theclass10, \
    label=['Knopoff PSV','Kennet PSV'])
    
TFDisplayTools.TFPlot(theclass3,theclass6,theclass8, \
    label=['1 layer','2 layers','6 layers'])
"""

"""
# kennet sh
print 'TF calculation using kennet sh method'
datash['inputtype'][0]='outcrop'
theclass11 = TFC(datash)
theclass11.tf_kennet_sh()
print 'calculation has been finished!'


TFDisplayTools.TFPlot(theclass1,theclass2,theclass3,theclass4,theclass5,theclass11, \
    label=['Kramer SH - Borehole','Knopoff Simple SH - Borehole','Knopoff Complete SH - Borehole', \
    'Knopoff PSV - Borehole','Kennet SH - Borehole','Kennet SH - Outcrop'])
    
TFDisplayTools.PhasePlot(theclass1,theclass2,theclass3,theclass4,theclass5,theclass11, \
    label=['Kramer SH - Borehole','Knopoff Simple SH - Borehole','Knopoff Complete SH - Borehole', \
    'Knopoff PSV - Borehole','Kennet SH - Borehole','Kennet SH - Outcrop'])
"""
"""
TFDisplayTools.TFPlot(theclass4, \
    label=['Knopoff PSV-S'],tfid=0)
    
TFDisplayTools.PhasePlot(theclass4, \
    label=['Knopoff PSV-S'],tfid=0)
    
TFDisplayTools.TFPlot(theclass4, \
    label=['Knopoff PSV-S'],tfid=1)
    
TFDisplayTools.PhasePlot(theclass4, \
    label=['Knopoff PSV-S'],tfid=1)
"""
"""
# test knopoff
# filename
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'

# input file reading
datash = IOfile.parsing_input_file(fname)

# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass2 = TFC(datash)
#theclass2.tf_knopoff_sh()
print 'calculation has been finished!'
"""
"""
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
    
# kramer
print 'TS calculatoin using kramer approach'
TS1 = TSC(fname,method='kramer286_sh')
print 'calculation has been finished!'

# knopoff sh
print 'TS calculation using simple knopoff approach'
TS2 = TSC(fname,method='knopoff_sh')
print 'calculation has been finished!'

# knopoff sh complete
print 'TS calculation using complete knopoff sh approach'
TS3 = TSC(fname,method='knopoff_sh_adv')
print 'calculation has been finished!'

# knopoff psv-s
print 'TS calculation using complete knopoff psv-s approach'
TS4 = TSC(fname2,method='knopoff_psv_adv')
print 'calculation has been finished!'

# kennet sh
print 'TS calculation using kennet sh method'
TS5 = TSC(fname,method='kennet_sh')
print 'calculation has been finished!'
"""