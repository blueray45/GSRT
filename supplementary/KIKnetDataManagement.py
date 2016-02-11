# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 01:45:57 2016

@author: irnakat
"""
# disable plot in ipython
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import tarfile
import os
import pickle
import numpy as np
import sys
sys.path.insert(0,'..')
from TSCalculator import TSCalculator as TSC
from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing
import scipy.signal as signal

foldername = '/media/data/VirtualMachines/SharedFolder/JAPANESEDATA/Japan_Data/alldata'
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
        if 'Lat.' in ln[:6]:
            eventLatitude = float(ln[18:])
        if 'Long.' in ln[:6]:
            eventLongitude = float(ln[18:])
        if 'Depth. (km)' in ln:
            eventDepth = float(ln[18:])
        if 'Mag.' in ln:
            eventMagnitude = float(ln[18:])
        if 'Station Lat.' in ln:
            stationLatitude = float(ln[18:])
        if 'Station Long.' in ln:
            stationLongitude = float(ln[18:])
        if 'Station Height(m)' in ln:
            stationHeight = float(ln[18:])
        if 'Sampling Freq(Hz)' in ln:
            samplingRate = float(ln[18:].replace('Hz',''))
        if 'Scale Factor' in ln:
            tmp = ln[18:]
            scalingfactor = eval(tmp[:17].replace('(gal)','.')+'.')
        if lnid>16:
            for dat in ln.split():
                seismo.append(float(dat)*scalingfactor*0.01) #0.01 is gal to m/s^2 conversion
        lnid+=1
    return seismo,[samplingRate,eventLatitude,eventLongitude,eventDepth,eventMagnitude,
        stationLatitude,stationLongitude,stationHeight]
    
def kiknet2ascii(fin,fout,plot=False,targetlength=24000):
    sig,signalProperties = kiknetdataread(fin)
    sR = signalProperties[0]
    # tapering & baseline correction
    window = signal.tukey(len(sig),alpha=0.15,sym=True)
    s = (np.array(sig)-np.mean(sig))*window
    # time vector
    lensignal = len(s)
    t = [i*(1./sR) for i in range(targetlength)]
    if plot:
        plt.figure()
        plt.plot(t,s-np.mean(s))
        plt.title(fout)
        plt.grid()
        plt.ylabel('$m/s^2$')
        plt.xlabel('Time (s)')
        #plt.ylim(-0.25,0.25)
    with open(fout,'wb') as f:
        for i in range(targetlength):
            if i<lensignal:
                f.write('%.6f\t%.6f\n'%(t[i],s[i]))
            else:
                f.write('%.6f\t%.6f\n'%(t[i],0.))
    return signalProperties
                
def butter_lowpass_filter(data, highcut, fs, order=4):
    # perform highpass filter of the signal
    def butter_lowpass(highcut, fs, order=4):
        nyq = 0.5 * fs
        high = highcut / nyq
        b, a = signal.butter(order, [high], btype='lowpass')
        return b, a
    b, a = butter_lowpass(highcut, fs, order=order)
    return signal.lfilter(b, a, data)
    
def rs_NigamJennings(data,dt,damping,period):
    from numpy import pi,sqrt,exp,sin,cos,zeros,zeros_like
    from numpy import max as maxnp
    from numpy import abs as absnp
    
    #data = data/100. # convert from cm to m
    data = np.array(data)
    #print('Starting NigamJennings Calculation...')
    Sa = zeros((len(period),len(damping)))
    for ii in range(len(damping)):
        for jj in range(len(period)):
            # simplyfying variable
            ksi = damping[ii]
            #print('Simplifying variable...')
            
            # a bit compensation for T=0 s which is impossible to calculate
            if period[jj]==0:
                period[jj]=period[jj+1]
                #print('Period 0 compensation...')
                
            # convert period to angular frequency
            w = 2.*pi/period[jj]
            #print('Converting to angular frequency...')
            
            # calculate damped period in angular frequency
            wd = w*sqrt(1-ksi**2)
            #print('Calculate damped period in angular frequency...')
            
            # parameter
            #print('Calculating parameter...')
            G = exp(-ksi*w*dt)
            F = sin(wd*dt)*G
            E = cos(wd*dt)*G
            
            # calculating alfa and beta value
            A11 = E+(ksi*w/wd)*F
            A12 = F/wd
            A21 = -w**2*A12
            A22 = E-(ksi*w/wd)*F
            B11 = (A11-1)/w**2
            B12 = (A12-2*ksi*w*B11-dt)/(w**2*dt)
            B21 = -A12
            B22 = B11/dt
            
            # Xi+1 = AXi + Bai
            # set up initial condition
            u = zeros((len(data),1))
            v = zeros((len(data),1))
            a = zeros((len(data),1))
            
            # loop for each timestep
            #print('Time stepping...')
            for kk in range(1,len(data)-1):
                u[kk] = A11*u[kk-1]+A12*v[kk-1]+B11*data[kk-1]+B12*data[kk]
                v[kk] = A21*u[kk-1]+A22*v[kk-1]+B21*data[kk-1]+B22*data[kk]
                a[kk] = -(2*ksi*w*v[kk-1]+w**2*u[kk-1])-data[kk-1]
                
            # find the max spectral values
            Sa[jj,ii] = maxnp(absnp(a)) # convert to g
            #print('Finishing Sa calculation...')
    return Sa

empiricalTF = [[],[]]
theoreticalTF = [[],[]]
pgareal = [[],[]]
pgatheory = [[],[]]
fftfreqTF = [[],[]]

# for i in range(len(targetextract[0])):
#for i in range(len(targetextract[0])-2,len(targetextract[0])):
#for i in range(500,501):
for i in range(2,3):
    print('processing data %d from total %d\n'%(i+1,len(targetextract[0])))
    if targetextract[1][i]>=mincomponent:
        tar = tarfile.open(targetextract[0][i])
        for member in tar.getmembers():
            for target in stationcomponent:
                tmp1 = member.name.find(target)
                tmp2 = member.name.find('UD1')
                tmp3 = member.name.find('UD2')
                
                if tmp1>=0 and tmp2>=0:
                    # managing seismogram
                    f=tar.extractfile(member)
                    dirprocess = member.name[tmp1:tmp2-1]
                    print dirprocess
                    if not os.path.isdir(member.name[tmp1:tmp2-1]):
                        os.mkdir(dirprocess)
                    signalProperties = kiknet2ascii(f,dirprocess+os.sep+member.name[tmp1:tmp2+3])
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
                    tsclass = TSC(dirprocess+os.sep+filepar,method='knopoff_sh',verbose=True)
                    
                    # draft plot
                    ts = []
                    ds = []
                    with open(dirprocess+os.sep+member.name[tmp1:tmp2+3],'rb') as fs:
                        for line in fs:
                            tmp = line.split()
                            ts.append(float(tmp[0]))
                            ds.append(float(tmp[1]))
                            
                    #calcsignal = butter_lowpass_filter(tsclass.time_series[0][:len(tsclass.inp_time)], 80., 200.)
                    calcsignal = tsclass.time_series[0][:len(tsclass.inp_time)]
                    fig3 = plt.figure()
                    a3 = fig3.add_subplot(111)
                    a3.plot(ts,ds-np.mean(ds),'g',label='real borehole recording')
                    a3.set_title(member.name[tmp1:tmp2-1])
                    a3.grid()
                    a3.set_ylabel('$m/s^2$')
                    a3.set_xlabel('Time (s)')                        
                    
                    fig1 = plt.figure()
                    a1 = fig1.add_subplot(111)
                    a1.plot(tsclass.inp_time,calcsignal,label='modeled surface recording')
                    #a1.set_ylim(-0.4,0.4)
                    a1.grid()
                    a1.set_ylabel('$m/s^2$')
                    a1.set_xlabel('Time (s)')
                    a1.set_title(dirprocess)
                    
                    # read amplification model
#                    mod_f = []
#                    mod_amp = []
#                    with open(masterdirpar+os.sep+target+'_model.amp') as famp:
#                        for line in famp:
#                            tmp = line.split()
#                            mod_f.append(tmp[0])
#                            mod_amp.append(tmp[1])
                    
#                    pgatheory[stationcomponent.index(target)].append(np.max(np.abs(calcsignal)))
                    tmp =rs_NigamJennings(calcsignal,ts[1]-ts[0],[0.05],[0.0001])
                    pgatheory[stationcomponent.index(target)].append(tmp[0,0])


                if tmp1>=0 and tmp3>=0:
                    # managing seismogram
                    f=tar.extractfile(member)
                    dirprocess = member.name[tmp1:tmp3-1]
                    if not os.path.isdir(member.name[tmp1:tmp3-1]):
                        os.mkdir(dirprocess)
                    signalProperties = kiknet2ascii(f,dirprocess+os.sep+member.name[tmp1:tmp3+3])
                    ts2 = []
                    ds2 = []
                    with open(dirprocess+os.sep+member.name[tmp1:tmp3+3],'rb') as fs:
                        for line in fs:
                            tmp = line.split()
                            ts2.append(float(tmp[0]))
                            ds2.append(float(tmp[1]))
                    
                    a1.plot(ts2,ds2-np.mean(ds2),'r',label='real surface recording',alpha=0.6)
                    a1.legend(loc='best')
                    fig1.tight_layout()
#                    fig1.savefig(dirprocess+os.sep+dirprocess+'_02.png')
                    
                    a3.plot(ts2,ds2-np.mean(ds2),'r',label='real surface recording',alpha=0.6)
                    a3.legend(loc='best')
                    fig3.tight_layout()
#                    fig3.savefig(dirprocess+os.sep+dirprocess+'_01.png')
                    
                    # calculation of observed transfer function
                    tf2 = np.abs(np.fft.fft(ds2-np.mean(ds2))/np.fft.fft(ds-np.mean(ds)))
                    tf2_freq = np.linspace(1./(ts2[-1]-ts2[0]),100.,len(ds2)/2)
                    tf2_smooth = konnoOhmachiSmoothing(tf2[:len(tf2_freq)],tf2_freq)
                    
                    # testing to treat the data
                    def datatreatment(inputdata,inputtime):
                        import scipy.signal as signal
                        # baseline correction
                        inputdata = inputdata-np.mean(inputdata)
                        # tapering
                        window = signal.tukey(len(inputdata),alpha=0.15,sym=True)
                        inputdata = inputdata*window
                        # fft
                        dt = inputtime[1]-inputtime[0]
                        fftdata = np.fft.fft(inputdata)
#                        fftfreq = np.fft.fftfreq(len(fftdata),dt)
                        # cut the length
                        fftdata = fftdata[:int(len(fftdata)/2)]
                        df = 1./(inputtime[-1]-inputtime[0])
                        fftfreq = np.array([(i)*df for i in range(len(fftdata))])
                        fftamp = np.abs(fftdata)
                        
                        #konno ohmachi smoothing
                        fftsmooth = konnoOhmachiSmoothing(fftamp,fftfreq)
                        
                        return fftamp,fftfreq,fftsmooth
                    
                    fftamp,fftfreq,fftsmooth = datatreatment(ds2,ts2)
                    fftamp2,fftfreq2,fftsmooth2 = datatreatment(ds,ts)
                    fftres = fftsmooth/fftsmooth2
                    empiricalTF[stationcomponent.index(target)].append(fftres)
                    fftfreqTF[stationcomponent.index(target)].append(fftfreq)
                    if not theoreticalTF[stationcomponent.index(target)]:
                        theoreticalTF[stationcomponent.index(target)].append(tsclass.tfclass.freq)
                        theoreticalTF[stationcomponent.index(target)].append(np.abs(tsclass.tfclass.tf[0]))
                        
                        
#                    pgareal[stationcomponent.index(target)].append(np.max(np.abs(ds2-np.mean(ds2))))
                    tmp = rs_NigamJennings(ds2,ts2[1]-ts2[0],[0.05],[0.0001])
                    pgareal[stationcomponent.index(target)].append(tmp[0,0])
                    
                    fig2 = plt.figure()
                    a2 = fig2.add_subplot(111)
                    a2.plot(tsclass.tfclass.freq,np.abs(tsclass.tfclass.tf[0]),'b',label='this study')
                    #a2.plot(mod_f,mod_amp,'r',label='from kik-net folder')
                    #a2.plot(tsclass.tfclass.freq,tf2[:len(tsclass.tfclass.freq)],'r',label='observation transfer function')
                    a2.plot(fftfreq,fftres,'r',label='observation transfer function')
                    a2.set_xscale('log')
                    a2.grid(which='both')
                    a2.set_ylabel('Amplification')
                    a2.set_title('Transfer Function - '+target)
                    a2.set_xlabel('Frequency (Hz)')
                    a2.legend(loc='best')
                    #a2.set_xlim(0.1,10.)
                    a2.set_yscale('log')
                    #a2.set_ylim(0.1,200)
                    a2.set_xlim(0.1,50)
                    a2.set_ylim(0.1,50)
                    fig2.tight_layout()
                    # fig2.savefig(dirprocess+os.sep+dirprocess+'_'+str(i)+'_03.png')
                    print('processing station %s with PGA calculated is %.4f and PGA observed is %.4f\n' \
                        %(target,pgatheory[stationcomponent.index(target)][-1],pgareal[stationcomponent.index(target)][-1]))
                plt.close("all")


with open('empiricalTF.bin','wb') as f:
    pickle.dump(empiricalTF, f)
    print('empiricalTF has been writen in file for further usage!')
    
with open('fftfreqTF.bin','wb') as f:
    pickle.dump(fftfreqTF, f)
    print('fftfreqTF has been writen in file for further usage!')
    
with open('theoreticalTF.bin','wb') as f:
    pickle.dump(theoreticalTF, f)
    print('theoreticalTF has been writen in file for further usage!')

#import seaborn as sns
for k in range(len(stationcomponent)):
    print 'Number of data for station %s is %d'%(stationcomponent[k],len(empiricalTF[k]))
    if empiricalTF[k]:
        empiricalTFFinal = []
        empiricalTFFinalstd = []
        nn=len(empiricalTF[k])
        for j in range(len(fftfreq)):
            tmp = []
            for i in range(nn):
                tmp.append(empiricalTF[k][i][j])
            tmp = np.log(np.array(tmp))
            empiricalTFFinal.append(np.mean(tmp))
            empiricalTFFinalstd.append(np.std(tmp))
            
        empiricalTFFinal = np.array(empiricalTFFinal)
        empiricalTFFinalstd = np.array(empiricalTFFinalstd)
        fig4 = plt.figure()
        a4 = fig4.add_subplot(111)
        a4.fill_between(fftfreq,10**(empiricalTFFinal-empiricalTFFinalstd),10**(empiricalTFFinal+empiricalTFFinalstd),facecolor='red',alpha=0.3)
        a4.plot(fftfreq,10**empiricalTFFinal,'r',label='observation transfer function')
        a4.plot(theoreticalTF[k][0],theoreticalTF[k][1],'b',label='this study')
        a4.set_xscale('log')
        a4.grid(which='both')
        a4.set_ylabel('Amplification')
        a4.set_title('Transfer Function - '+stationcomponent[k])
        a4.set_xlabel('Frequency (Hz)')
        a4.legend(loc='best')
        a4.set_yscale('log')
        fig4.tight_layout()
        # fig4.savefig('TF'+stationcomponent[k]+'.png')
        
        fig5 = plt.figure(figsize=(5,5))
        a5 = fig5.add_subplot(111)
        a5.plot([0,1], [0,1],color='grey',linestyle='--')
        a5.plot(pgareal[k],pgatheory[k],'bo')
        a5.grid(which='both')
        a5.set_ylabel('PGA calculated ($m/s^2$)')
        a5.set_title('PGA - '+stationcomponent[k])
        a5.set_xlabel('PGA observation ($m/s^2$)')
        a5.set_aspect('equal')
        a5.set_xlim(0,np.max([np.max(pgareal),np.max(pgatheory)]))
        a5.set_ylim(0,np.max([np.max(pgareal),np.max(pgatheory)]))
        fig5.tight_layout()
        # fig5.savefig('PGA'+stationcomponent[k]+'.png')