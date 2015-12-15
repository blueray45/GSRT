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

TFDisplayTools.TFPlot(theclass1,theclass2,theclass3,theclass4,theclass5, \
    label=['Kramer SH','Knopoff Simple SH','Knopoff Complete SH','Knopoff PSV','Kennet SH'])
    
TFDisplayTools.PhasePlot(theclass1,theclass2,theclass3,theclass4,theclass5, \
    label=['Kramer SH','Knopoff Simple SH','Knopoff Complete SH','Knopoff PSV','Kennet SH'])