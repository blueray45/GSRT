# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 20:07:31 2015
Common lists used on GSRT
@author: irnakat
"""

mode = ['sh-linear-elastic',
        'sh-linear-viscoelastic',
        'sh-linear-elastic-adv',
        'sh-linear-viscoelastic-adv',
        'sh-kennet',
        'psv-p-linear-elastic-adv',
        'psv-p-linear-viscoelastic-adv',
        'psv-s-linear-elastic-adv',
        'psv-s-linear-viscoelastic-adv',
        'psv-p-kennet',
        'psv-s-kennet',
        'sh-linear-equivalent',
        'psv-linear-equivalent',
        'non-linear']
        
def modescan(inp):
    if inp.find('sh')>0:
        case = '1'
    elif inp.find('psv_s')>0:
        case = '2'
    elif inp.find('psv_p')>0:
        case = '3'
    else:
        case = '0'
    if inp.find('viscoelastic')>0:
        elastic = '2'
    else: # elastic
        elastic = '1'
    if inp.find('kramer'):
        method = '2'
    elif inp.find('knopoff'):
        method = '3'
    elif inp.find('kennet'):
        method = '4'
    else:
        method = '1'    # method auto
    return case+elastic+method