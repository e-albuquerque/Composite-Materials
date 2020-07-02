#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 20:29:26 2020

@author: eder
"""
import numpy as np
import laminado as lm
import matplotlib.pyplot as plt

plt.close("all")
# Propriedades elásticas da lâmina
EL=148e9
ET=10.5e9
GLT=5.61e9
nuLT=0.3
espessura=2e-3


# ET=EL # Nos materiais isotrópicos EL=ET
# GLT=EL/(2*(1+nuLT)) # GLT para materiais isotrópicos

# Material=[EL ET GLT nuLT espessura theta]
# Número de linhas da matiz Material é igual ao número de camadas do laminado

Material=np.array([[EL, ET, GLT, nuLT, espessura, 45],
           [EL, ET, GLT, nuLT, espessura, 90],
           [EL, ET, GLT, nuLT, espessura, 90],
           [EL, ET, GLT, nuLT, espessura, 45]])

Nx=10e3
Ny=0
Nxy=0
Mx=0
My=0
Mxy=0
NM=np.array([Nx, Ny, Nxy, Mx, My, Mxy])

ABD=lm.calc_ABD(Material)
epgamma=np.linalg.inv(ABD).dot(NM)
epsilonx=epgamma[0]
epsilony=epgamma[1]
gammaxy=epgamma[2]
kappax=epgamma[3]
kappay=epgamma[4]
kappaxy=epgamma[5]


# Número de nós que serão gerados nas direções x e y
N=20

# NOS=[coordenada x do no, coordenada y do no, coordenada z do no]
# u = deslocamentos dos nos na direção x
# v = deslocamentos dos nos na direção v
# w = deslocamentos dos nos na direção w

NOS,uvw=lm.gera_NOSuvw(N,epgamma)

fator=0.1 # máxima deformação da malha deformada em relação a distmaxima

lm.plota_deslocamentos(NOS,uvw,fator)

lm.plota_tensoes(epgamma,Material)