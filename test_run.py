# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:55:33 2015

@author: irnakat
"""

# exporting external helper module
import numpy as np
import pylab as plt
import commondata as cd

# test IOfile
import IOfile
from TFCalculator import TFCalculator as TFC
import TFDisplayTools

# single layer test case

# filename
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'

# input file reading
datash = IOfile.parsing_input_file(fname)
datapsv = IOfile.parsing_input_file(fname2)

# kramer
print 'TF calculatoin using kramer approach'
theclass1 = TFC(datash)
theclass1.tf_kramer286_sh() # check/verufy kramer calculation
print 'calculation has been finished!'

# knopoff sh
print 'TF calculation using simple knopoff approach'
theclass2 = TFC(datash)
theclass2.tf_knopoff_sh()
print 'calculation has been finished!'

# knopoff sh complete
print 'TF calculation using complete knopoff sh approach'
theclass3 = TFC(datash)
theclass3.tf_knopoff_sh_adv()
print 'calculation has been finished!'

# knopoff psv-s
print 'TF calculation using complete knopoff psv-s approach'
theclass4 = TFC(datapsv)
theclass4.tf_knopoff_psv_adv()
print 'calculation has been finished!'

# kennet sh
print 'TF calculation using kennet sh method'
theclass5 = TFC(datash)
theclass5.tf_kennet_sh()
print 'calculation has been finished!'

# kennet psv-s
print 'TF calculation using kennet psv-s approach'
theclass10 = TFC(datapsv)
theclass10.tf_kennett_psv()
print 'calculation has been finished!'
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
"""

"""
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
"""

TFDisplayTools.TFPlot(theclass4,theclass10, \
    label=['Knopoff PSV','Kennet PSV'])
    
TFDisplayTools.PhasePlot(theclass4,theclass10, \
    label=['Knopoff PSV','Kennet PSV'])