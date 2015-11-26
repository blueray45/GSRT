# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:55:33 2015

@author: irnakat
"""

# exporting external helper module
import numpy as np

# test IOfile
import IOfile
import sh_tf_fix

fname = 'sampleinput_linear_elastic_1layer_halfspace.dat'
mode = 'linear-elastic'

data = IOfile.parsing_input_file(fname,mode)
freq = np.linspace(0.1,50.,100)
#tf = sh_tf_fix.sh_tf_fix(hl,vs,dn,qs,freq)
print data