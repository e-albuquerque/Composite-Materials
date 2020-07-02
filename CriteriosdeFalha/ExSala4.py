#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 16:09:21 2020

@author: eder
"""
import numpy as np
import CriteriosdeFalha as CF


k=2 # Fator de segurança
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
theta=np.array([0,90])

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
Nxy=0.
Mx=0.
My=0.
Mxy=0.
NM=np.array([Nx, Ny, Nxy, Mx, My, Mxy])

epgamma=np.linalg.inv(ABD).dot(NM)
Q=CF.calc_Q(EL,ET,GLT,nuLT)
Qbar1=CF.calc_Qbar(Q,theta[0])
Qbar2=CF.calc_Qbar(Q,theta[1])
sig1xy=Qbar1.dot(epgamma[0:3])
sig2xy=Qbar2.dot(epgamma[0:3])


# Lâmina a theta[0] grau
T1=CF.calc_T(theta[0])
sig1LT=T1.dot(sig1xy)
p1=CF.calc_pressao(sig1LT[0],sig1LT[1],sig1LT[2],xt,xc,yt,yc,s,k)
print("Pressão que provoca a falha da lâmina com fibras a ",theta[0]," grau: ",p1," MPa")


# Lâmina a theta[1] grau
T2=CF.calc_T(theta[1])
sig2LT=T2.dot(sig2xy)
p2=CF.calc_pressao(sig2LT[0],sig2LT[1],sig2LT[2],xt,xc,yt,yc,s,k)
print("Pressão que provoca a falha da lâmina com fibras a ",theta[1]," grau: ",p2," MPa")



if(p2>p1):
    plim=p1
else:
    plim=p2
p=plim*10 # Pressão em bar

print("Pressão em bar que provoca a falha do tubo: ",p," bar")
