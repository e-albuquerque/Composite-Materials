#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 16:09:21 2020

@author: eder
"""
import numpy as np
import CriteriosdeFalha as CF

xt=1260.
xc=2500.
yt=61.
yc=202.
s=67.
Q=np.array([[20., .7, 0],
            [.7, 2, 0],
            [0, 0, .7]])
Qs=np.zeros((3,3,3))
thetas=np.array([45, 0, 45])
Qs[:,:,0]=Q
Qs[:,:,1]=Q
Qs[:,:,2]=Q
hs=np.array([3, 6, 3])
ABD=CF.calc_ABD2(Qs,thetas,hs)
Nx=1000.
Ny=200.
Nxy=0.
Mx=0.
My=0.
Mxy=0.
NM=np.array([Nx, Ny, Nxy, Mx, My, Mxy])
epgamma=np.linalg.inv(ABD).dot(NM)
Qbar45=CF.calc_Qbar(Q,45)
sig45xy=Qbar45.dot(epgamma[0:3])
sig0xy=Q.dot(epgamma[0:3])

# Fator de segurança por Tsai-Wu:

# Lâmina a 0 grau
sig0LT=sig0xy
k=CF.calc_TsaiWu(sig0LT[0],sig0LT[1],sig0LT[2],xt,xc,yt,yc,s)
print("Fator de segurança de Tsai-Wu da lâmina com fibras a 0 grau: ",k)

# Lâmina a 45 graus
T45=CF.calc_T(45)
sig45LT=T45.dot(sig45xy)
k=CF.calc_TsaiWu(sig45LT[0],sig45LT[1],sig45LT[2],xt,xc,yt,yc,s)
print("Fator de segurança de Tsai-Wu da lâmina com fibras a 45 graus: ",k)


# Fator de segurança por Tsai-Hill:

# Lâmina a 0 grau
k=CF.calc_TsaiHill(sig0LT[0],sig0LT[1],sig0LT[2],xt,xc,yt,yc,s)
print("Fator de segurança de Tsai-Hill da lâmina com fibras a 0 grau: ",k)

# Lâmina a 45 graus
k=CF.calc_TsaiHill(sig45LT[0],sig45LT[1],sig45LT[2],xt,xc,yt,yc,s)
print("Fator de segurança de Tsai-Hill da lâmina com fibras a 45 graus: ",k)