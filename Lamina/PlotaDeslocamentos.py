#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 22:13:04 2020

@author: eder
"""
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import lamina as lm
import numpy as np

plt.close("all")

# Número de nós que serão gerados nas direções x e y
N=6

# inicializa a matriz NOS
# A matriz nós terá 2 colunas. Na primeira coluna
# armazernará as coordenadas x dos nós e na segunda coluna
# as coordenadas y. O número de linhas é N^2
NOS=np.zeros((N**2,2))

# Inicializa os vetores deslocamentos
u=np.zeros(N**2) # deslocamentos na direção x
v=np.zeros(N**2) # deslocamentos na direção y

# coordenadas mínimas e máximas da lâmina
xmin=-1
xmax=1
ymin=-1
ymax=1

# distância entre oois pontos da malha
deltax=(xmax-xmin)/(N-1) # distância na direção x
deltay=(ymax-ymin)/(N-1) # distância na direção y

# Tensões atuantes na lâmina
sigx=36
sigy=0
tauxy=0

# Propriedades elásticas da lâmina
EL=148.
ET=10.5
GLT=5.61
nuLT=0.3


# ângulo de inclinação das fibras
theta=45

# Calcula as deformações epsilonx, epsilony e gammaxy
tensoesxy=np.array([sigx,sigy,tauxy])
Q=lm.calc_Q(EL,ET,GLT,nuLT)
T=lm.calc_T(theta)
Qbar=lm.calc_Qbar(Q,theta)
tensoesLT=T.dot(tensoesxy)
sigL=tensoesLT[0]
sigT=tensoesLT[1]
tauLT=tensoesLT[2]
defxy=np.linalg.inv(Qbar).dot(tensoesxy)
epsilonx=defxy[0]
epsilony=defxy[1]
gammaxy=defxy[2]	


ino=0 # Inicializa o contador de nós

# Atribui os valores aos elementos da matrz NOS e dos vetores u e v
x=xmin # inicializa x
for i in range(N):
    y=ymin # inicializa y
    for j in range(N):
        NOS[ino,0]=x
        NOS[ino,1]=y
        v[ino]=epsilony*y+gammaxy*x/2
        u[ino]=epsilonx*x+gammaxy*y/2
        y=y+deltay # atualiza y
        ino=ino+1 # atualiza o contador de nos
    x=x+deltax # atualiza x

# Cria a variável escala que será aplicada sobre os deslocamentos
distmaxima=np.sqrt((xmax-xmin)**2. +(ymax-ymin)**2)
deslocamento_total=np.sqrt(u**2+v**2)
deslocmaximo=np.max(np.abs(deslocamento_total))
escala=distmaxima/deslocmaximo*.1


# Cria os triângulos que são armazenados na matriz elem
triang = tri.Triangulation(NOS[:,0], NOS[:,1])
elem=triang.triangles


plt.figure()
# plota a malha não deformada
plt.triplot(NOS[:,0], NOS[:,1], elem, color=(0.0,0.,0.),linewidth=0.4,linestyle="--")
# plota a malha deformada
plt.triplot(NOS[:,0]+escala*u,NOS[:,1]+escala*v,color=(1.0,0.,0.),linewidth=0.8,linestyle="-")
