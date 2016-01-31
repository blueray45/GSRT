# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 13:01:35 2016

@author: irnakat
"""
# 3n+1 problem

a = 20
b = 10

c = 0

while b != 0:
    if (b&1) != 0:
        c=c+a
    a<<1
    b>>1
    
print c