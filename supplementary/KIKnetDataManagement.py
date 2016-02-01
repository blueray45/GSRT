# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 01:45:57 2016

@author: irnakat
"""
import tarfile
import os
import pickle
import numpy as np
import sys
sys.path.insert(0,'..')
from TSCalculator import TSCalculator as TSC

foldername = '/media/DATA/VirtualMachines/SharedFolder/JAPANESEDATA/Japan_Data/alldata'
stationcomponent = ['HYGH10','FKSH12']
case = 'parlinear'
masterdirpar = 'masterparameter'
mincomponent = 6

targetextract = [[],[]]
idfind = 0

outname = ''.join(stationcomponent)
if not os.path.isfile(outname):
    listyear = os.listdir(foldername)
    for year in listyear:
        listmonth = os.listdir(foldername+os.sep+year)
        for month in listmonth:
            print 'processing --> %s%s'%(year,month)
            listfile = os.listdir(foldername+os.sep+year+os.sep+month)
            for thefile in listfile:
                try:
                    tar = tarfile.open(foldername+os.sep+year+os.sep+month+os.sep+thefile)
                    filecount = 0
                    for member in tar.getnames():
                        for target in stationcomponent:
                            if member.find(target)>0:
                                filecount+=1
                    if filecount>0:
                        targetextract[0].append(foldername+os.sep+year+os.sep+month+os.sep+thefile)
                        targetextract[1].append(filecount)
                        idfind+=1
                except:
                    print('Error reading tar file. Move on.')
    print('Total %d events with intended station has been found!'%idfind)
    with open(outname,'wb') as f:
        pickle.dump(targetextract, f)
        print('data has been writen in file for further usage!')
    
else:
    print('previous scanning is detected! Loading previous data instead!')
    with open(outname,'rb') as f:      
        targetextract = pickle.load(f)
        
def kiknetdataread(fileobject):
    seismo = []
    lnid=0
    for ln in f:
        if lnid==10:
            samplingRate = float(ln[18:].replace('Hz',''))
        if lnid==13:
            tmp = ln[18:]
            scalingfactor = eval(tmp[:17].replace('(gal)','.')+'.')
        if lnid>16:
            for dat in ln.split():
                seismo.append(float(dat)*scalingfactor*0.01) #0.01 is gal to m/s^2 conversion
        lnid+=1
    return seismo,samplingRate
    
def kiknet2ascii(fin,fout,plot=False):
    import pylab as plt
    sig,sR = kiknetdataread(fin)
    s = sig
    lensignal = len(s)
    t = [i*(1./sR) for i in range(lensignal)]
    if plot:
        plt.figure()
        plt.plot(t,s)
        plt.title(fout)
        plt.grid()
        plt.ylabel('$m/s^2$')
        plt.xlabel('Time (s)')
    with open(fout,'wb') as f:
        for i in range(lensignal):
            f.write('%.6f\t%.6f\n'%(t[i],s[i]))

#for i in range(len(targetextract[0])):
#for i in range(500,501):
for i in range(2,3):
    if targetextract[1][i]>=mincomponent:
        tar = tarfile.open(targetextract[0][i])
        for member in tar.getmembers():
            for target in stationcomponent:
                tmp1 = member.name.find(target)
                tmp2 = member.name.find('UD2')
                if tmp1>=0 and tmp2>=0:
                    # managing seismogram
                    f=tar.extractfile(member)
                    dirprocess = member.name[tmp1:tmp2-1]
                    if not os.path.isdir(member.name[tmp1:tmp2-1]):
                        os.mkdir(dirprocess)
                    kiknet2ascii(f,dirprocess+os.sep+member.name[tmp1:tmp2+3],True)
                    # managing parameter
                    tmp = os.listdir(masterdirpar+os.sep)
                    for filepar in tmp:
                        if filepar.find(target)>=0 and filepar.find(case)>=0:
                            with open(masterdirpar+os.sep+filepar,'rb') as fp:
                                with open(dirprocess+os.sep+filepar,'wb') as fp2:
                                    for ln in fp:
                                        fp2.write(ln.replace('%variable%',dirprocess+os.sep+member.name[tmp1:tmp2+3]))
                            break
                    # running tscalculator
                    print dirprocess+os.sep+filepar
                    tsclass = TSC(dirprocess+os.sep+filepar,method='knopoff_sh',verbose=True)
                    plt.figure()                    
                    plt.plot(tsclass.inp_time,tsclass.time_series[0][:len(tsclass.inp_time)])
                    plt.ylim(-30,30)
                    plt.grid()
                    plt.ylabel('$m/s^2$')
                    plt.xlabel('Time (s)')
                    plt.title(dirprocess)
                    
                    # read amplification model
                    mod_f = []
                    mod_amp = []
                    with open(masterdirpar+os.sep+target+'_model.amp') as famp:
                        for line in famp:
                            tmp = line.split()
                            mod_f.append(tmp[0])
                            mod_amp.append(tmp[1])
                    
                    plt.figure()
                    plt.plot(tsclass.tfclass.freq,np.abs(tsclass.tfclass.tf[0]),'b',label='this study')
                    plt.plot(mod_f,mod_amp,'r',label='from kik-net folder')
                    plt.xscale('log')
                    plt.grid(which='both')
                    plt.ylabel('Amplification')
                    plt.title(target)
                    plt.xlabel('Frequency (Hz)')
                    plt.legend(loc='best')