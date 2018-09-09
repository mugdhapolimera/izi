#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 16:39:26 2018

@author: mugpol
"""
a = np.arange(1,5)
b = np.arange(1,5)
print b[np.where(a <=4)]
print max([1,2,3] or [0])
print max( list(b[np.where(a <=4)]) or [0])