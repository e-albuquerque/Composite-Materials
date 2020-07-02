#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 16:34:07 2020

@author: eder
"""


#A thin-walled pressure vessel is manufactured by a filament winding method using glass/epoxy prepregs. Find the optimum angles, θ, if the pressure vessel is made of [±θ] ns sublaminate with 
# 1. Spherical construction for maximum strength
# 2. Cylindrical construction for maximum strength
# 3. Cylindrical construction for no change in the internal diameter
# Apply Tsai–Wu failure theory. Use properties of unidirectional glass/epoxy lamina from Table 2.2.


import numpy as np
import CriteriosdeFalha as CF
import matplotlib.pyplot as plt

xt=1725.
xc=1350.
yt=40.
yc=275.
s=95.

Nx=500.
Ny=1000.
Nxy=0.
Mx=0.
My=0.
Mxy=0.
NM=np.array([Nx, Ny, Nxy, Mx, My, Mxy])

# Propriedades elásticas da lâmina
EL=20000. #MPa
ET=2000.
GLT=700.
nuLT=0.35
espessura=2. # mm
n=50;
theta=np.linspace(0,90,n)
vetk1=np.zeros(n)
vetk2=np.zeros(n)
theta1_otimo=0;
theta2_otimo=0;
# ET=EL # Nos materiais isotrópicos EL=ET
# GLT=EL/(2*(1+nuLT)) # GLT para materiais isotrópicos

# Material=[EL ET GLT nuLT espessura theta]
# Número de linhas da matiz Material é igual ao número de camadas do laminado
kmax=0
for i in range(n):
    theta1=theta[i]
    theta2=-theta[i]
    Material=np.array([[EL, ET, GLT, nuLT, espessura, theta[i]],
           [EL, ET, GLT, nuLT, espessura, -theta[i]],
           [EL, ET, GLT, nuLT, espessura, -theta[i]],
           [EL, ET, GLT, nuLT, espessura,theta[i]]])

    ABD=CF.calc_ABD(Material)

    epgamma=np.linalg.inv(ABD).dot(NM)

    Q=CF.calc_Q(EL,ET,GLT,nuLT)
    Qbar1=CF.calc_Qbar(Q,theta1)
    Qbar2=CF.calc_Qbar(Q,theta2)
    sig1xy=Qbar1.dot(epgamma[0:3])
    sig2xy=Qbar2.dot(epgamma[0:3])

    # Fator de segurança por Tsai-Wu:

    # Lâmina a 0 grau
    T1=CF.calc_T(theta1)
    sig1LT=T1.dot(sig1xy)
    k1=CF.calc_TsaiWu(sig1LT[0],sig1LT[1],sig1LT[2],xt,xc,yt,yc,s)
    
    # Lâmina a 45 graus
    T2=CF.calc_T(theta2)
    sig2LT=T2.dot(sig2xy)
    k2=CF.calc_TsaiWu(sig2LT[0],sig2LT[1],sig2LT[2],xt,xc,yt,yc,s)
    
    if(k1<k2):
        if(k1>kmax):
            kmax=k1
            theta1_otimo=theta1
            theta2_otimo=theta2
    else:
        if(k2>kmax):
            kmax=k2
            theta1_otimo=theta1
            theta2_otimo=theta2
    vetk1[i]=k1
    vetk2[i]=k2

print("kmax = ",kmax)
print("theta1_otimo = ",theta1_otimo)
print("theta2_otimo = ",theta2_otimo)
plt.close("all")
plt.figure()
plt.plot(theta,vetk1,linestyle="-",marker="*",label="k1")
plt.plot(theta,vetk2,linestyle="--",marker="o",label="k2")
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="upper left",fancybox="true") # Create a legend of all the existing plots using their labels as names