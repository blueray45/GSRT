# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 03:24:06 2016

@author: irnakat
"""
import pylab as plt
import numpy as np
# test read file

def kiknetdataread(fname):
    components = ['EW1','NS1','UD1','EW2','NS2','UD2']
    seismo = []
    for comp in components:
        thefile = fname+'.'+comp
        lnid = 0
        seismo.append([])
        with open(thefile,'rb') as f:
            for ln in f:
                if lnid==10:
                    samplingRate = float(ln[18:].replace('Hz',''))
                if lnid==13:
                    tmp = ln[18:]
                    scalingfactor = eval(tmp[:17].replace('(gal)','.')+'.')
                if lnid>16:
                    for dat in ln.split():
                        seismo[-1].append(float(dat)*scalingfactor*0.01) #0.01 is gal to m/s^2 conversion
                lnid+=1
    return seismo,samplingRate
    
def kiknet2ascii(fin,fout,outid=5):
    sig,sR = kiknetdataread(fin)
    s = sig[outid]
    lensignal = len(s)
    t = [i*(1./sR) for i in range(lensignal)]
    plt.figure()
    plt.plot(t,s)
    with open(fout,'wb') as f:
        for i in range(lensignal):
            f.write('%.6f\t%.6f\n'%(t[i],s[i]))
    
    

fname = 'AICH049710081906'
#seismo = kiknetdataread(fname)
kiknet2ascii(fname,fname)
print np.shape(seismo)