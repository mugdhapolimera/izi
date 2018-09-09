#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 13:46:56 2018

@author: mugpol
"""
import pickle
f3 = open("/afs/cas.unc.edu/users/m/u/mugpol/Documents/IZI/izi/IZI_results.pkl", "rb") 

while True:
    try:
        d = pickle.load(f3)
        print '\n\n', d
    except EOFError:
        print 'End of File'
        break
    
f3.close()