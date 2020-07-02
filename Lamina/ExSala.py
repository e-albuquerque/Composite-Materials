#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:20:02 2020

@author: eder
"""

import lamina as lm
import numpy as np

# Tensões atuantes nas lâminas

# Tensões atuantes na lâmina
sigx=-3.5e6
sigy=7e6
tauxy=-1.4e6
tensoesxy=np.array([sigx,sigy,tauxy]) # vetor que contém as 3 componentes de tensões

# Propriedades elásticas da lâmina
EL=14.e9
nuLT=0.4
ET=3.5e9
GLT=4.2e9

# ângulo de inclinação das fibras
theta=60

# Calcula as tensões sigL, sigT e tauLT
tensoesxy=np.array([sigx,sigy,tauxy])
T=lm.calc_T(theta)
tensoesLT=T.dot(tensoesxy)
sigL=tensoesLT[0]
sigT=tensoesLT[1]
tauLT=tensoesLT[2]

# Calcula as deformações epsilonL, epsilonT e gammaLT
Q=lm.calc_Q(EL,ET,GLT,nuLT)
defLT=np.linalg.inv(Q).dot(tensoesLT)
epsilonL=defLT[0]
epsilonT=defLT[1]
gammaLT=defLT[2]

# Calcula as deformações epsilonx, epsilony e gammaxy 
# Primeira forma de cálculo
Qbar=lm.calc_Qbar(Q,theta)
defxy=np.linalg.inv(Qbar).dot(tensoesxy)
epsilonx=defxy[0]
epsilony=defxy[1]
gammaxy=defxy[2]	

# # Segunda forma de cálculo (forma usada na sala)
defxy2=np.linalg.inv(T).dot([epsilonL,epsilonT,gammaLT/2])
epsilonx2=defxy2[0]
epsilony2=defxy2[1]
gammaxy2=2*defxy2[2]	