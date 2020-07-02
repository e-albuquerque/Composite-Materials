#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:22:53 2020

@author: eder
"""
import laminado as lm
import numpy as np
Q=np.array([[20,.7,0],[.7,2,0],[0,0,.7]])
Qs=np.zeros((3,3,3))
thetas=np.array([45, 0, 45])
Qs[:,:,0]=Q
Qs[:,:,1]=Q
Qs[:,:,2]=Q
hs=np.array([3,6,3])
ABD=lm.calc_ABD2(Qs,thetas,hs)
Nx=1000
Ny=200
Nxy=0
Mx=0
My=0
Mxy=0
NM=np.array([Nx, Ny, Nxy, Mx, My, Mxy])
epgamma=np.linalg.inv(ABD).dot(NM)
Qbar45=lm.calc_Qbar(Q,45)
sig45=Qbar45.dot(epgamma[0:3])
sig0=Q.dot(epgamma[0:3])