# -*- coding: utf-8 -*-
"""
Created on Mon Dec 07 13:26:28 2015

Display package
to display TF and friends nicely and well done

@author: irnakat
"""
import numpy as np
import pylab as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.cm as cm

def colorcycle(ncolorinput):
    ncolor = 256
    clist = cm.rainbow(np.arange(ncolor))
    
    cclist = []
    if ncolorinput>=2:
        for i in range(ncolorinput-1):
            n = ncolor/(ncolorinput-1)
            cclist.append(clist[n*i])
    cclist.append(clist[-1])
    return cclist
    
def velocityprofileplot(inp,a2,cclist):
        xmax = []
        xmin = []
        for i in range(len(inp)):
            hl = inp[i].hl
            depthtemp = np.concatenate(([0.],np.cumsum(hl)))
            try:
                vptemp = inp[i].vp
                novp = False
            except:
                print('vp is not found for argument : %d!'%i)
                novp = True
            vstemp = inp[i].vs
            depth = [0.]
            vs = [vstemp[0]/1000.]
            if not novp:
                vp = [vptemp[0]/1000.]
            for j in range(1,len(hl)):
                depth.append(depthtemp[j])
                depth.append(depthtemp[j])
                vs.append(vstemp[j-1]/1000.)
                vs.append(vstemp[j]/1000.)  
                if not novp:
                    vp.append(vptemp[j-1]/1000.)
                    vp.append(vptemp[j]/1000.) 
                  
            depth.append(depth[-1]+0.1*depth[-1])
            vs.append(vs[-1])
            if not novp:                   
                vp.append(vp[-1])
                if i==0:
                    a2.plot(vp,depth,color=cclist[i],linestyle=':',lw=3.,label='vp')
                    a2.plot(vs,depth,color=cclist[i],lw=3.,label='vs')
                else:
                    a2.plot(vp,depth,color=cclist[i],linestyle=':',lw=3.)
                    a2.plot(vs,depth,color=cclist[i],lw=3.)
                xmin.append(np.min(vs))
                xmax.append(np.max(vp))
            else:
                if i==0:
                    a2.plot(vs,depth,color=cclist[i],lw=3.)
                else:
                    a2.plot(vs,depth,color=cclist[i],lw=3.)
                xmin.append(np.min(vs))
                xmax.append(np.max(vs))
                
        a2.set_ylim(np.min(depth),np.max(depth)) 
        a2.set_xlim(min(xmin)-0.05,max(xmax)+0.05)
        a2.legend(loc='best',fancybox=True,framealpha=0.5)
        a2.invert_yaxis()
        a2.grid(True)
        a2.set_xlabel('Velocity (km/s)')
        a2.set_ylabel('Depth (m)')
        a2.set_title('Velocity profile')

def TFPlot(*arg,**kwargs):
    # check given axis name
    try:
        axname = kwargs['axname']
    except KeyError:
        axname = None
    try:
        label = kwargs['label']
    except KeyError:
        label = ['']*(len(arg))
    # create new figure is axis name is not given
    if axname==None:
        f = plt.figure(figsize=(14.,7.),dpi=300)
        
        # create default color cycle
        cclist = colorcycle(len(arg))        
        
        a2= f.add_subplot(1,5,5)
        velocityprofileplot(arg,a2,cclist)
        
        a = f.add_subplot(1,5,(1,4))
    else:
        a = axname
        
    # set label and properties of axis
    a.set_xlabel('Frequency (Hz)')
    a.set_ylabel('Amplification')
    a.set_yscale('log')
    a.set_xscale('log')
    a.grid(True,which='major',color='k')
    a.grid(True,which='minor',color='grey')
    a.minorticks_on()
    a.tick_params(axis='both', which='major', labelsize=11, labelcolor='k')
    a.tick_params(axis='both', which='minor', labelsize=10, labelcolor='grey')
    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(ScalarFormatter())
    
    # check number of input
    #if len(arg)%2!=0:
    #    raise InputError('Number or input pairs if not correct! Detected input pairs is %.1f.'%len(arg)/2.)
    # check length of each pairs and plot data
    minx = []; maxx = []; miny = []; maxy = []
    tfid=0
    for i in range(len(arg)):    
        freq = arg[i].freq
        tf = arg[i].tf       
        a.plot(freq,np.abs(tf[tfid]),color=cclist[i],label=label[i])
        minx.append(np.min(arg[i].freq))
        maxx.append(np.max(arg[i].freq))
        miny.append(np.min(np.abs(arg[i].tf[tfid])))
        maxy.append(np.max(np.abs(arg[i].tf[tfid])))
        
    minmax = [[np.min(minx),np.max(maxx)],
              [np.min(miny),np.max(maxy)]]
    percspan = 0.1
    minmaxspan = [np.exp(np.log(minmax[0][1]-minmax[0][0])*percspan),
                  np.exp(np.log(minmax[1][1]-minmax[1][0])*percspan)]
    a.set_xlim(minmax[0][0],minmax[0][1])
    a.set_ylim(minmax[1][0]-minmaxspan[1],minmax[1][1]+minmaxspan[1])
    a.legend(loc='best',fancybox=True,framealpha=0.5)
    f.tight_layout()
    
def PhasePlot(*arg,**kwargs):
    # check given axis name
    try:
        axname = kwargs['axname']
    except KeyError:
        axname = None
    try:
        label = kwargs['label']
    except KeyError:
        label = ['']*(len(arg))
    # create new figure is axis name is not given
    if axname==None:
        f = plt.figure(figsize=(14.,7.),dpi=300)
        
        # create default color cycle
        cclist = colorcycle(len(arg))    
        
        a2= f.add_subplot(1,5,5)
        velocityprofileplot(arg,a2,cclist)
        
        a = f.add_subplot(1,5,(1,4))
    else:
        a = axname
    
    # set label and properties of axis
    a.set_xlabel('Frequency (Hz)')
    a.set_ylabel('Phase (Rad)')
    #a.set_yscale('log')
    a.set_xscale('log')
    a.grid(True,which='major',color='k')
    a.grid(True,which='minor',color='grey')
    #a.minorticks_on()
    a.tick_params(axis='both', which='major', labelsize=11, labelcolor='k')
    #a.tick_params(axis='both', which='minor', labelsize=10, labelcolor='grey')
    for axis in [a.xaxis, a.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
        #axis.set_minor_formatter(ScalarFormatter())
    
    # check length of each pairs and plot data
    tfid = 0
    minx = []; maxx = []; miny = []; maxy = []
    for i in range(len(arg)):
        freq = arg[i].freq
        tf = arg[i].tf      
        y = np.angle(tf[tfid])
        a.plot(freq,y,label=label[i],color=cclist[i])
        minx.append(np.min(arg[i].freq))
        maxx.append(np.max(arg[i].freq))
        miny.append(np.min(y))
        maxy.append(np.max(y))
        
    minmax = [[np.min(minx),np.max(maxx)],
              [np.min(miny),np.max(maxy)]]
    percspan = 0.1
    minmaxspan = [np.exp(np.log(minmax[0][1]-minmax[0][0])*percspan),
                  (minmax[1][1]-minmax[1][0])*percspan]
    a.set_xlim(minmax[0][0],minmax[0][1])
    a.set_ylim(minmax[1][0]-minmaxspan[1],minmax[1][1]+minmaxspan[1])
    a.legend(loc='best',fancybox=True,framealpha=0.5)
    f.tight_layout()
    
def SpectroPlot(data,nx=100,ny=100,ylabel='incidence angle',zlabel='Amplification',cmap='rainbow'):
    import scipy.interpolate
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import MaxNLocator
    
    if data['sensitivity']==False:
        raise IOError("Sensitivity calculation hasn't been performed! Please do so.")
    
    x = np.asarray(data['x'])
    y = np.asarray(data['y'])
    z = np.asarray(data['z'])
    
    xi,yi = np.logspace(np.log10(x.min()),np.log10(x.max()),nx), np.linspace(y.min(),y.max(),ny)
    xi,yi = np.meshgrid(xi,yi)
    
    # interpolation
    zi = scipy.interpolate.griddata((x,y),z,(xi,yi),method='linear')
    
    # plot the data
    f = plt.figure(figsize=(16.,8.),dpi=300)
    
    # plot velocity profile
    #print data
    #print data.hl
    #print data['hl']
    a2= f.add_subplot(1,5,5)
    if type(data['hl'][0])==float:
        xmin = []
        xmax = []
        # plot vp
        if data['modeID']>=5:
            if type(data['vp'][0])==list:
                lvp = len(data['vp'][0])
                clistvp = cm.cool(np.arange(lvp))
                for i in range(len(data['vp'][0])):
                    vptemp = np.concatenate(([data['vp'][0][i]],data['vp'][1:]))
                    depthtemp = np.concatenate(([0.],np.cumsum(data['hl'])))
                    
                    depth = [0.]
                    vp = [vptemp[0]/1000.]
                    for j in range(1,len(depthtemp)-1):
                        depth.append(depthtemp[j])
                        vp.append(vptemp[j-1]/1000.)
                        depth.append(depthtemp[j])
                        vp.append(vptemp[j]/1000.)   
                    depth.append(depth[-1]+0.1*depth[-1])
                    vp.append(vp[-1])                    
                    
                    a2.plot(vp,depth,color=clistvp[i])
                    xmax.append(np.max(vp))
            else:
                vptemp = data['vp']
                depthtemp = np.concatenate(([0.],np.cumsum(data['hl'])))
                
                depth = [0.]
                vp = [vptemp[0]/1000.]
                for j in range(1,len(depthtemp)-1):
                    depth.append(depthtemp[j])
                    vp.append(vptemp[j-1]/1000.)
                    depth.append(depthtemp[j])
                    vp.append(vptemp[j]/1000.)   
                depth.append(depth[-1]+0.1*depth[-1])
                vp.append(vp[-1])                
                
                a2.plot(vp,depth,color='b')
                xmax.append(np.max(vp))
        # plot vs
        if type(data['vs'][0])==list:
            lvs = len(data['vs'][0])
            clistvs = cm.hot(np.arange(lvs))
            for i in range(len(data['vs'][0])):
                vstemp = np.concatenate(([data['vs'][0][i]],data['vs'][1:]))
                depthtemp = np.concatenate(([0.],np.cumsum(data['hl'])))
                
                depth = [0.]
                vs = [vstemp[0]/1000.]
                for j in range(1,len(depthtemp)-1):
                    depth.append(depthtemp[j])
                    vs.append(vstemp[j-1]/1000.)
                    depth.append(depthtemp[j])
                    vs.append(vstemp[j]/1000.)   
                depth.append(depth[-1]+0.1*depth[-1])
                vs.append(vs[-1])                     
                
                a2.plot(vs,depth,color=clistvs[i])
                xmin.append(np.min(vs))
                xmax.append(np.max(vs))
        else:
            vstemp = data['vs']
            depthtemp = np.concatenate(([0.],np.cumsum(data['hl'])))
            
            depth = [0.]
            vs = [vstemp[0]/1000.]
            for j in range(1,len(depthtemp)-1):
                depth.append(depthtemp[j])
                vs.append(vstemp[j-1]/1000.)
                depth.append(depthtemp[j])
                vs.append(vstemp[j]/1000.)   
            depth.append(depth[-1]+0.1*depth[-1])
            vs.append(vs[-1])             

            a2.plot(vs,depth,color='r')
            xmin.append(np.min(vs))
            xmax.append(np.max(vs))
        a2.set_xlim(np.min(xmin)-0.05,np.max(xmax)+0.05)
        a2.set_ylim(np.min(depth),np.max(depth))
    else:
        if data['modeID']>=5:
            ld = len(data['hl'][0])
            clistvp = cm.cool(np.arange(ld))
            clistvs = cm.hot(np.arange(ld))
            for i in range(len(data['hl'][0])):
                hl = np.concatenate(([data['hl'][0][i]],data['hl'][1:]))
                depthtemp = np.concatenate(([0.],np.cumsum(hl)))
                vptemp = data['vp']
                vstemp = data['vs']
                depth = [0.]
                vs = [vstemp[0]/1000.]
                vp = [vptemp[0]/1000.]
                for j in range(1,len(hl)):
                    depth.append(depthtemp[j])
                    vs.append(vstemp[j-1]/1000.)
                    vp.append(vptemp[j-1]/1000.)
                    depth.append(depthtemp[j])
                    vs.append(vstemp[j]/1000.)  
                    vp.append(vptemp[j]/1000.)   
                depth.append(depth[-1]+0.1*depth[-1])
                vs.append(vs[-1])
                vp.append(vp[-1])

                a2.plot(vp,depth,color=clistvp[i])
                a2.plot(vs,depth,color=clistvs[i])
                a2.set_xlim(np.min(vs)-0.05,np.max(vp)+0.05)
                a2.set_ylim(np.min(depth),np.max(depth))
        else:
            ld = len(data['hl'][0])
            clistvs = cm.hot(np.arange(ld))
            for i in range(len(data['hl'][0])):
                hl = np.concatenate(([data['hl'][0][i]],data['hl'][1:]))
                vstemp = data['vs']
                depthtemp = np.concatenate(([0.],np.cumsum(hl)))
                depth = [0.]
                vs =[vstemp[0]/1000.]
                for j in range(1,len(hl)):
                    depth.append(depthtemp[j])
                    vs.append(vstemp[j-1]/1000.)
                    depth.append(depthtemp[j])
                    vs.append(vstemp[j]/1000.)     
                depth.append(depth[-1]+0.1*depth[-1])
                vs.append(vs[-1])
                a2.plot(vs,depth,color=clistvs[i])
                a2.set_xlim(np.min(vs)-0.05,np.max(vs)+0.05)
                a2.set_ylim(np.min(depth),np.max(depth))
    
    a2.invert_yaxis()
    a2.set_xlabel('Velocity (km/s)')
    a2.set_ylabel('Depth (m)')
    a2.set_title('Velocity profile')
    # plot data
    a = f.add_subplot(1,5,(1,4))
    #zi = np.log10(zi)
    #z = np.log10(z)
    am = a.imshow(zi, vmin=0.1, vmax=z.max(), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()],
             aspect = 'auto',norm=LogNorm())

    a.set_xlabel('Frequency (Hz)')
    a.set_ylabel(ylabel)
    a.set_xscale('log')
    a.minorticks_on()
    a.tick_params(axis='x', which='major', labelsize=11, labelcolor='k')
    a.tick_params(axis='x', which='minor', labelsize=10, labelcolor='grey')
    for axis in [a.xaxis]:
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(ScalarFormatter())
    cb = plt.colorbar(am,label=zlabel)
    cb.locator = MaxNLocator(20)
    cb.formatter = ScalarFormatter()
    cb.update_ticks()
    

"""
# test display
import IOfile
from TFCalculator import TFCalculator as TFC
from sensitivityTools import sensitivityTools as sT

fname = 'sampleinput_linear_elastic_6layer_halfspace.dat'
data = IOfile.parsing_input_file(fname)
theclass = TFC(data)
tf1 = theclass.tf_kramer286_sh()

fname2 = 'sampleinput_linear_elastic_1layer_halfspace.dat'
data2 = IOfile.parsing_input_file(fname2)
theclass2 = TFC(data2)
tf2 = theclass2.tf_kramer286_sh()

TFPlot(theclass,theclass2,label=['six layers','one layer'])

PhasePlot(theclass,theclass2)
"""
"""
fname3 = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
data3 = IOfile.parsing_input_file(fname2)
x = np.array([])
y = np.array([])
z = np.array([])
ianglist = np.linspace(0.0,90.0,91)
for i in range(len(ianglist)):
    data3['iang'] = np.deg2rad(ianglist[i])
    theclass3 = TFC(data3)
    tf3 = theclass3.tf_knopoff_sh_adv()
    #print data3['iang']
    #print np.shape(np.abs(tf3[0][:,0])),np.shape(z)
    x = np.concatenate((x,theclass3.freq))
    z = np.concatenate((z,np.abs(tf3[0])))
    y = np.concatenate((y,np.zeros_like(theclass3.freq)+ianglist[i]))
    
SpectroPlot(x,y,z,nx=100,ny=100,ylabel='incidence angle',zlabel='Amplification',cmap='rainbow')
"""