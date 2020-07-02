#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 22:28:20 2020

@author: eder
"""
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import lamina as lm
import numpy as np


plt.close("all")


# Número de pontos usados para plotar os gráficos
npontos=100

# Cria um vetor theta com npontos de zero a 360
theta=np.linspace(0,360,npontos)

# Inicializa as variáveis que serão plotadas
sigL=np.zeros(npontos)
sigT=np.zeros(npontos)
tauLT=np.zeros(npontos)
epsilonL=np.zeros(npontos)
epsilonT=np.zeros(npontos)
gammaLT=np.zeros(npontos)
epsilonx=np.zeros(npontos)
epsilony=np.zeros(npontos)
gammaxy=np.zeros(npontos)

# Tensões atuantes nas lâminas
sigx=36
sigy=10
tauxy=-2
tensoesxy=np.array([sigx,sigy,tauxy]) # vetor que contém as 3 componentes de tensões

# Propriedades elásticas da lâmina
EL=148.
nuLT=0.3

ET=10.5
GLT=5.61

EL=ET # Nos materiais isotrópicos EL=ET
GLT=EL/(2*(1+nuLT)) # GLT para materiais isotrópicos

# Calcula a matriz de rigidez na direção LT
Q=lm.calc_Q(EL,ET,GLT,nuLT)


# Loop que calcula os vetores que serão plotados
for i in range(npontos):
	T=lm.calc_T(theta[i]) # Matriz de transformação
	tensoesLT=T.dot(tensoesxy) # Tensões no sistema LT
	sigL[i]=tensoesLT[0]
	sigT[i]=tensoesLT[1]
	tauLT[i]=tensoesLT[2]
	defLT=np.linalg.inv(Q).dot(tensoesLT) # deformações no sistema LT
	epsilonL[i]=defLT[0]
	epsilonT[i]=defLT[1]
	gammaLT[i]=defLT[2]	
	Qbar=lm.calc_Qbar(Q,theta[i]) # Matriz de rigidez no sistema xy	
	defxy=np.linalg.inv(Qbar).dot(tensoesxy) # Deformações no sistema xy
	epsilonx[i]=defxy[0]
	epsilony[i]=defxy[1]
	gammaxy[i]=defxy[2]	

# Plota os gráficos
plt.figure()
plt.plot(theta,sigL,linestyle="-",marker="x",label="$\sigma_L$")
plt.plot(theta,sigT,linestyle="-",marker="o",label="$\sigma_T$")
plt.plot(theta,tauLT,linestyle="-",marker="p",label=r"$\tau_{LT}$")
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="lower right",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel(r"$\theta$")
plt.ylabel("Tensões")
plt.title(r"Tensões no sistema LT em função de $\theta$")

plt.figure()
plt.plot(theta,epsilonx,linestyle="-",marker="x",label="$\epsilon_x$")
plt.plot(theta,epsilony,linestyle="-",marker="o",label="$\epsilon_y$")
plt.plot(theta,gammaxy,linestyle="-",marker="p",label="$\gamma_{xy}$")
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="lower right",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel(r"$\theta$")
plt.ylabel("Deformações")
plt.title(r"Deformações no sistema $xy$ em função de $\theta$")

plt.figure()
plt.plot(theta,epsilonL,linestyle="-",marker="x",label="$\epsilon_L$")
plt.plot(theta,epsilonT,linestyle="-",marker="o",label="$\epsilon_T$")
plt.plot(theta,gammaLT,linestyle="-",marker="p",label="$\gamma_{LT}$")
plt.axis("tight") # Fit the axis tightly to the plot
plt.grid("on")
plt.legend(loc="lower right",fancybox="true") # Create a legend of all the existing plots using their labels as names
plt.xlabel(r"$\theta$")
plt.ylabel("Deformações")
plt.title(r"Deformações no sistema $LT$ em função de $\theta$")