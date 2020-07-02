#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 14:28:08 2020

@author: eder
"""

import numpy as np
import CriteriosdeFalha as CF
import matplotlib.pyplot as plt


plt.close("all")
xt=1725.
xc=1350.
yt=40.
yc=275.
s=95.

sigmaxy=np.array([200,  -100, 50])
npontos=100

thetas = np.linspace(0,90,npontos)

f=np.zeros(npontos)
f2=np.zeros(npontos)
for i in range(npontos):
    theta=thetas[i]
    T=CF.calc_T(theta)
    sigmalt=T.dot(sigmaxy)
    sigmal=sigmalt[0]
    sigmat=sigmalt[1]
    tault=sigmalt[2]
    # Tsai-Wu
    f[i]=CF.calc_TsaiWu(sigmal,sigmat,tault,xt,xc,yt,yc,s)
    # Tsai-Hill
    f2[i]=CF.calc_TsaiHill(sigmal,sigmat,tault,xt,xc,yt,yc,s)

plt.figure()
plt.plot(thetas,f,linestyle="-",marker="x",label="Tsai-Wu")
plt.plot(thetas,f2,linestyle="--",marker="o",label="Tsai-Hill")
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="upper center",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel(r"$\theta$")
plt.ylabel("Índice de falha")
plt.title("Comparação dos critérios de falha")
plt.show()