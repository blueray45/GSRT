# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 04:00:41 2016

@author: irnakat
"""

# test IOfile
import IOfile
from TFCalculator import TFCalculator as TFC
import TFDisplayTools

# validity test for SH PSV case using S wave as an input

# filename
fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
fname2 = 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat'

# input file reading
datash = IOfile.parsing_input_file(fname)
datapsvs = IOfile.parsing_input_file(fname2)

# kramer
print 'TF calculatoin using kramer approach'
theclass1 = TFC(datash)
theclass1.tf_kramer286_sh() # check/verify kramer calculation
print 'calculation has been finished!'

# knopoff sh complete
print 'TF calculation using complete knopoff sh approach'
theclass3 = TFC(datash)
theclass3.tf_knopoff_sh_adv()
print 'calculation has been finished!'

# knopoff psv-s
print 'TF calculation using complete knopoff psv-s approach'
theclass4 = TFC(datapsvs)
theclass4.tf_knopoff_psv_adv()
print theclass4.tf[1][19]
print 'calculation has been finished!'

# kennet sh
print 'TF calculation using kennet sh method'
theclass5 = TFC(datash)
theclass5.tf_kennet_sh()
print 'calculation has been finished!'

TFDisplayTools.TFPlot(theclass1,theclass3,theclass5,theclass4, \
    label=['Kramer SH','Knopoff SH','Kennet SH','Knopoff PSV'])
    
TFDisplayTools.PhasePlot(theclass1,theclass3,theclass5,theclass4, \
    label=['Kramer SH','Knopoff SH','Kennet SH','Knopoff PSV'])