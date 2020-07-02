#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:35:36 2020

@author: eder
"""
import matplotlib.pyplot as plt
import HalpinTsai as HT
import numpy as np

plt.close("all")
npontos=20

Vf = np.linspace(0,1,npontos)

Ef=233
Em=4.62
nuf=0.2
num=0.36
Gm = Em/(2*(1 + num))
Gf= Ef/(2*(1+nuf))
EL=HT.calc_EL(Vf,Ef,Em)
nuLT=HT.calc_nuLT(Vf,nuf,num)
GLT=HT.calc_GLT(Vf,Gf,Gm,1)
GLT2=HT.calc_GLT(Vf,Gf,Gm,2)
ET=HT.calc_ET(Vf,Ef,Em,1)
ET2=HT.calc_ET(Vf,Ef,Em,2)


plt.figure()

p1 = plt.plot(Vf,EL,linestyle='-',marker='x',label='$E_L$')
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="upper center",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel("$V_f$")
plt.ylabel("$E_L$")
plt.title("Módulo de Elasticidade Longitudinal")
plt.show()

plt.figure()
p1 = plt.plot(Vf,nuLT,linestyle="-",marker="x",label='$\nu_{LT}$')
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="upper center",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel("$V_f$")
plt.ylabel("$\nu_{LT}$")
plt.title("Coeficiente de Poisson principal")
plt.show()


plt.figure()
p1 = plt.plot(Vf,GLT,linestyle="-",marker="x",label="$G_{LT}$ simplificado")
varlabel="$G_{LT}$ Halpin-Tsai com $E_f$ =" + str(Ef)
p2 = plt.plot(Vf,GLT2,linestyle="--",marker="o",label=varlabel)
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="upper center",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel("$V_f$")
plt.ylabel("$G_{LT}$")
plt.title("Módulo de cisalhamento")
plt.show()

plt.figure()
p1 = plt.plot(Vf,ET,linestyle="-",marker="x",label="$E_T$ simplificado")
p1 = plt.plot(Vf,ET2,linestyle="--",marker="o",label="$E_T$ Halpin-Tsai")
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="upper center",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel("$V_f$")
plt.ylabel("$E_T$")
plt.title("Módulo de elasticidade transversal")
plt.show()