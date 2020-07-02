#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 08:59:19 2017

@author: eder
"""
import numpy as np
a = float(input('a = '))
b = float(input('b = '))
c = float(input('c = '))

delta=b**2-4*a*c
print('\nDelta: ',delta)
if(delta>=0):
    x1=(-b+np.sqrt(delta))/(2*a)
    x2=(-b-np.sqrt(delta))/(2*a)
    print('x1: ',x1)
    print('x2: ',x2)
else:
    print('Não existem raizes reais!!')