# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 02:34:44 2016

@author: irnakat
"""

import pickle
import pylab as plt
import sys
import numpy as np
sys.path.insert(0,'..')
from TSCalculator import TSCalculator as TSC

print('previous scanning is detected! Loading previous data instead!')

with open('empiricalTF.bin','rb') as f:
    empiricalTF = pickle.load(f)
    print('empiricalTF has been read in file for further usage!')
    
theoreticalTF = [[[],[]],[[],[]]]
    
# calculate theoreticalTF for FKSH12
fname = '/media/data/gdriveubuntu/Kuliah/Semester III/04-Individual_Project/Site_Response_Toolkit/02-Code/Github/GSRT/FKSH120007120129/FKSH12_parlinear.dat'
FKSH12 = TSC(fname,method='knopoff_sh',verbose=True)
theoreticalTF[1][0] = FKSH12.tfclass.freq
theoreticalTF[1][1] = np.abs(FKSH12.tfclass.tf[0])

# calculate theoreticalTF for HYGH10
fname = '/media/data/gdriveubuntu/Kuliah/Semester III/04-Individual_Project/Site_Response_Toolkit/02-Code/Github/GSRT/HYGH100101200724/HYGH10_parlinear.dat'
HYGH10 = TSC(fname,method='knopoff_sh',verbose=True)
theoreticalTF[0][0] = HYGH10.tfclass.freq
theoreticalTF[0][1] = np.abs(HYGH10.tfclass.tf[0])


fftfreq = HYGH10.tfclass.freq[:len(HYGH10.tfclass.freq)/2]
stationcomponent = ['HYGH10','FKSH12']
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
            tmp = np.log10(np.array(tmp))
            empiricalTFFinal.append(np.mean(tmp))
            empiricalTFFinalstd.append(np.std(tmp))
            
        empiricalTFFinal = np.array(empiricalTFFinal)
        empiricalTFFinalstd = np.array(empiricalTFFinalstd)
        #empiricalTFFinal = np.array(np.log(empiricalTF[k][19]))
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
        a4.set_xlim(0.1,50)
        a4.set_ylim(0.1,50)
        fig4.tight_layout()
        fig4.savefig('TF'+stationcomponent[k]+'.png')
        
#        fig5 = plt.figure(figsize=(5,5))
#        a5 = fig5.add_subplot(111)
#        a5.plot([0,1], [0,1],color='grey',linestyle='--')
#        a5.plot(pgareal[k],pgatheory[k],'bo')
#        a5.grid(which='both')
#        a5.set_ylabel('PGA calculated ($m/s^2$)')
#        a5.set_title('PGA - '+stationcomponent[k])
#        a5.set_xlabel('PGA observation ($m/s^2$)')
#        a5.set_aspect('equal')
#        a5.set_xlim(0,np.max([np.max(pgareal),np.max(pgatheory)]))
#        a5.set_ylim(0,np.max([np.max(pgareal),np.max(pgatheory)]))
#        fig5.tight_layout()
#        fig5.savefig('PGA'+stationcomponent[k]+'.png')