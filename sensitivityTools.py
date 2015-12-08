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

def sensitivityTools(fname,sensitivity='incidence angle',method='tf_kramer286_sh'):
    import IOfile
    import numpy as np
    import pylab as plt
    from TFCalculator import TFCalculator as TFC
    from TFDisplayTools import SpectroPlot
    
    # timing module
    import time
    
    data = IOfile.parsing_input_file(fname)
    elapsed = []
    x = np.array([])
    y = np.array([])
    z = np.array([])
    
    if sensitivity=='incidence angle':
        ianglist = np.linspace(0.0,90.0,91)
        for i in range(len(ianglist)):
            data['iang'] = np.deg2rad(ianglist[i])
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'incidence angle'
    elif sensitivity=='incidence angle phase':
        ianglist = np.linspace(0.0,90.0,91)
        for i in range(len(ianglist)):
            data['iang'] = np.deg2rad(ianglist[i])
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq))
            z = np.concatenate((z,np.rad2deg(np.angle(tf[0]))))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'incidence angle'
    elif sensitivity=='thickness':
        #ianglist = np.linspace(1.0,500.0,500)
        ianglist = np.linspace(1.0,20000.0,500)        
        for i in range(len(ianglist)):
            data['hl'] = [ianglist[i],0.]
            theclass = TFC(data)
            start = time.clock()
            tf = eval('theclass.'+method+'()')
            elapsed.append(time.clock()-start)
            x = np.concatenate((x,theclass.freq))
            z = np.concatenate((z,np.abs(tf[0])))
            y = np.concatenate((y,np.zeros_like(theclass.freq)+ianglist[i]))
            ylabel = 'thickness (m)'
            
    print('average elapsed time : %.6f with standard deviation : %.6f'%(np.mean(elapsed),np.std(elapsed)))
    SpectroPlot(x,y,z,nx=100,ny=100,ylabel=ylabel,cmap='rainbow')
    plt.title('Sensitivity analysis on %s using %s method'%(sensitivity,method))


# incidence angle
print('incidence angle sensitivity analysis')
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='incidence angle',method='tf_kramer286_sh')
sensitivityTools(fname,sensitivity='incidence angle',method='tf_knopoff_sh')
sensitivityTools(fname,sensitivity='incidence angle',method='tf_knopoff_sh_adv')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='incidence angle',method='tf_knopoff_psv_adv')

"""
# incidence angle --> plot phase
print('incidence angle sensitivity analysis')
fname = 'sampleinput_linear_elastic_1layer_halfspace_adv.dat'
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_kramer286_sh')
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_knopoff_sh')
sensitivityTools(fname,sensitivity='incidence angle phase',method='tf_knopoff_sh_adv')
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'
sensitivityTools(fname2,sensitivity='incidence angle phase',method='tf_knopoff_psv_adv')
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