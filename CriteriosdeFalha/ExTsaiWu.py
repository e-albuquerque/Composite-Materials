#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:58:40 2020

@author: eder
"""
import numpy as np
import CriteriosdeFalha as CF

xt=1725.
xc=1350.
yt=40.
yc=275.
s=95.


# Propriedades elásticas da lâmina
EL=20000. #MPa
ET=2000.
GLT=700.
nuLT=0.35
espessura=1. # mm
theta=np.array([0, 90])

# ET=EL # Nos materiais isotrópicos EL=ET
# GLT=EL/(2*(1+nuLT)) # GLT para materiais isotrópicos

# Material=[EL ET GLT nuLT espessura theta]
# Número de linhas da matiz Material é igual ao número de camadas do laminado
Material=np.array([[EL, ET, GLT, nuLT, espessura, theta[0]],
           [EL, ET, GLT, nuLT, espessura, theta[1]],
           [EL, ET, GLT, nuLT, espessura, theta[1]],
           [EL, ET, GLT, nuLT, espessura,theta[0]]])
ABD=CF.calc_ABD(Material)


Nx=500. 
Ny=1000.
# Nx=500. *0.8236534127487375
# Ny=1000. *0.8236534127487375
Nxy=0.
Mx=0.
My=0.
Mxy=0.
NM=np.array([Nx, Ny, Nxy, Mx, My, Mxy])
epgamma=np.linalg.inv(ABD).dot(NM)
num_theta=Material.shape[0]
theta=Material[:,5]
Q=CF.calc_Q(EL,ET,GLT,nuLT)

for i in range(num_theta):
    # Lâmina a theta1 grau
    Qbar=CF.calc_Qbar(Q,theta[i])
    sigxy=Qbar.dot(epgamma[0:3])
    T=CF.calc_T(theta[i])
    sigLT=T.dot(sigxy)
    k=CF.calc_TsaiWu(sigLT[0],sigLT[1],sigLT[2],xt,xc,yt,yc,s)
    print("k = ",k)