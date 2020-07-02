#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 22:51:20 2020

@author: eder
"""

import laminado as lm
import numpy as np
Q=np.array([[20,.7,0],[.7,2,0],[0,0,.7]])
Qs=np.zeros((3,3,2))
thetas=np.array([45, 0])
Qs[:,:,0]=Q
Qs[:,:,1]=Q
hs=np.array([3,5])
ABD=lm.calc_ABD2(Qs,thetas,hs)