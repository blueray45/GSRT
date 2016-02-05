# -*- coding: utf-8 -*-
"""
Created on Tue Dec 08 01:49:36 2015
@author: irnakat

Sensitivity analysis toolkit
This module allows to perform sensitivity analysis in terms of :
    - incidence angle (given constant velocity and thickness)
    - thickness of layer (given constant incidence angle and velocity)
    - Vs (given constant incidence angle and thickness)
"""

def sensitivityTools(fname,sensitivity='incidence angle',senslinspace=[0.,90.,91],method='tf_kramer286_sh'):
    import IOfile
    import numpy as np
    import pylab as plt
    from TFCalculator import TFCalculator as TFC
    from TFDisplayTools import SpectroPlot
    
    # timing module
    import time
    
    data = IOfile.parsing_input_file(fname)
    data['sensitivity'] = True
    elapsed = []
    x = np.array([])
    y = np.array([])
    z = np.array([])
    ianglist = np.linspace(senslinspace[0],senslinspace[1],senslinspace[2])
    
    if sensitivity=='incidence angle':
        for i in range(len(ianglist)):
            data['iang'] = np.deg2rad(ianglist[i])
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'incidence angle'
            #print np.min(theclass.freq),np.min(np.abs(tf[0]))
        data['iang']=ianglist
            
    elif sensitivity=='incidence angle phase':
        for i in range(len(ianglist)):
            data['iang'] = np.deg2rad(ianglist[i])
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.rad2deg(np.angle(tf[0]))))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'incidence angle'
        data['iang']=ianglist
        
    elif sensitivity=='incidence angle vectorial':
        for i in range(len(ianglist)):
            data['iang'] = np.deg2rad(ianglist[i])
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.sqrt(np.abs(tf[0])**2+np.abs(tf[1])**2)))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'incidence angle'
        data['iang']=ianglist
            
    elif sensitivity=='thickness':   
        for i in range(len(ianglist)):
            data['hl'] = [ianglist[i],0.]
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'thickness (m)'
        data['hl'][0]=ianglist.tolist()
            
    elif sensitivity[:2]=='vp':   
        tmp = sensitivity.split()
        for i in range(len(ianglist)):
            data['vp'][int(tmp[1])-1] = ianglist[i]
            data['comp']='p'
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[1])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'vp (m/s)'
        data['vp'][int(tmp[1])-1]=ianglist.tolist()
            
    elif sensitivity[:2]=='vs':   
        tmp = sensitivity.split()
        for i in range(len(ianglist)):
            data['vs'][int(tmp[1])-1] = ianglist[i]
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'vs (m/s)'
        data['vs'][int(tmp[1])-1]=ianglist.tolist()
        
    elif sensitivity[:2]=='qp':   
        tmp = sensitivity.split()
        for i in range(len(ianglist)):
            data['qp'][int(tmp[1])-1] = ianglist[i]
            data['comp']='p'
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[1])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = '$Q_{P}$'
        data['qp'][int(tmp[1])-1]=ianglist.tolist()
        
    elif sensitivity[:2]=='qs':   
        tmp = sensitivity.split()
        for i in range(len(ianglist)):
            data['qs'][int(tmp[1])-1] = ianglist[i]
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = '$Q_{S}$'
        data['qs'][int(tmp[1])-1]=ianglist.tolist()
        
    elif sensitivity[:2]=='dn':   
        tmp = sensitivity.split()
        for i in range(len(ianglist)):
            data['dn'][int(tmp[1])-1] = ianglist[i]
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq+theclass.freq[1]))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = '$\\rho (kg/m^3)$'
        data['dn'][int(tmp[1])-1]=ianglist.tolist()
            
            
    data['x'] = x
    data['y'] = y
    data['z'] = z

    print('method : %s; average elapsed time : %.6f with standard deviation : %.6f'%(method,np.mean(elapsed),np.std(elapsed)))
    SpectroPlot(data,nx=100,ny=100,ylabel=ylabel,cmap='rainbow')
    #plt.title('Sensitivity analysis on %s using %s method'%(sensitivity,method))


"""
# incidence angle --> plot phase
print('incidence angle sensitivity analysis for phase plot')
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_kramer286_sh')
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_knopoff_sh')
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_knopoff_sh_adv')
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_kennet_sh')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='incidence angle phase',method='tf_knopoff_psv_adv')
"""
"""
# incidence angle
print('incidence angle sensitivity analysis for transfer function plot')
#fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
#sensitivityTools(fname,sensitivity='incidence angle',senslinspace=[87.,89.,88],method='tf_kramer286_sh')
#sensitivityTools(fname,sensitivity='incidence angle',senslinspace=[87.,89.,88],method='tf_knopoff_sh')
#sensitivityTools(fname,sensitivity='incidence angle',senslinspace=[89.5,90.,88],method='tf_knopoff_sh_adv')
#sensitivityTools(fname,sensitivity='incidence angle',senslinspace=[89.5,90.,88],method='tf_kennet_sh')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='incidence angle',senslinspace=[0.,87.,88],method='tf_knopoff_psv_adv')
"""

# Thickness of layer
print('thickness sensitivity analysis')
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
#sensitivityTools(fname,sensitivity='thickness',senslinspace=[1.,100.,100],method='tf_kramer286_sh')
sensitivityTools(fname,sensitivity='thickness',senslinspace=[1000.,10000.,100],method='tf_knopoff_sh')
#sensitivityTools(fname,sensitivity='thickness',senslinspace=[1.,100.,100],method='tf_knopoff_sh_adv')
#sensitivityTools(fname,sensitivity='thickness',senslinspace=[1.,100.,100],method='tf_kennet_sh')
#fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
#sensitivityTools(fname2,sensitivity='thickness',senslinspace=[1.,100.,100],method='tf_knopoff_psv_adv')

"""
# Vs of first layer
print('Vs sensitivity analysis')
#fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
#sensitivityTools(fname,sensitivity='vs 1',senslinspace=[50.,1000.,81],method='tf_kramer286_sh')
#sensitivityTools(fname,sensitivity='vs 1',senslinspace=[50.,1000.,81],method='tf_knopoff_sh')
#sensitivityTools(fname,sensitivity='vs 1',senslinspace=[50.,1000.,81],method='tf_knopoff_sh_adv')
#sensitivityTools(fname,sensitivity='vs 1',senslinspace=[50.,1000.,81],method='tf_kennet_sh')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='vs 1',senslinspace=[50.,1000.,81],method='tf_knopoff_psv_adv')
"""
"""
# Vp of first layer
print('Vp sensitivity analysis')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='vp 1',senslinspace=[1500.,3000.,81],method='tf_knopoff_psv_adv')
"""
"""
# Qp of first layer
print('Qp sensitivity analysis')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='qp 1',senslinspace=[1.,100.,81],method='tf_knopoff_psv_adv')
"""
"""
# Qs of first layer
print('Qs sensitivity analysis')
#fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
#sensitivityTools(fname,sensitivity='qs 1',senslinspace=[1.,100.,81],method='tf_kramer286_sh')
#sensitivityTools(fname,sensitivity='qs 1',senslinspace=[1.,100.,81],method='tf_knopoff_sh')
#sensitivityTools(fname,sensitivity='qs 1',senslinspace=[1.,100.,81],method='tf_knopoff_sh_adv')
#sensitivityTools(fname,sensitivity='qs 1',senslinspace=[1.,100.,81],method='tf_kennet_sh')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='qs 1',senslinspace=[1.,100.,81],method='tf_knopoff_psv_adv')
"""
"""
# Density of first layer
print('Density sensitivity analysis')
#fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
#sensitivityTools(fname,sensitivity='dn 1',senslinspace=[500.,3000.,81],method='tf_kramer286_sh')
#sensitivityTools(fname,sensitivity='dn 1',senslinspace=[500.,3000.,81],method='tf_knopoff_sh')
#sensitivityTools(fname,sensitivity='dn 1',senslinspace=[500.,3000.,81],method='tf_knopoff_sh_adv')
#sensitivityTools(fname,sensitivity='dn 1',senslinspace=[500.,3000.,81],method='tf_kennet_sh')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='dn 1',senslinspace=[500.,3000.,81],method='tf_knopoff_psv_adv')
"""
"""
# incidence angle vectorial
print('incidence angle sensitivity analysis for transfer function plot')
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='incidence angle',senslinspace=[0.,80.,91],method='tf_knopoff_sh_adv')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='incidence angle vectorial',senslinspace=[0.,80.,91],method='tf_knopoff_psv_adv')
"""
"""
# incidence angle
print('incidence angle sensitivity analysis')
fname = 'sampleinput_linear_elastic_6layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='incidence angle',method='tf_kramer286_sh')
sensitivityTools(fname,sensitivity='incidence angle',method='tf_knopoff_sh')
sensitivityTools(fname,sensitivity='incidence angle',method='tf_knopoff_sh_adv')
fname2 = 'sampleinput_psv_s_linear_elastic_6layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='incidence angle',method='tf_knopoff_psv_adv')
"""
"""
# thickness of the layer
print('thickness sensitivity analysis')
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='thickness',method='tf_kramer286_sh')
sensitivityTools(fname,sensitivity='thickness',method='tf_knopoff_sh')
sensitivityTools(fname,sensitivity='thickness',method='tf_knopoff_sh_adv')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='thickness',method='tf_knopoff_psv_adv')
"""
"""
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='thickness',senslinspace=[10.,1000.,100],method='tf_kramer286_sh')
"""