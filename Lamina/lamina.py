# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
import numpy as np
def  calc_T(theta):
    theta = theta*np.pi/180
    m = np.cos(theta)
    n = np.sin(theta)
    T = np.array([[m**2,n**2,2*m*n],
         [n**2,  m**2,  -2*m*n],
        [-m*n,  m*n , m**2-n**2]])
    return T
def  calc_Qbar(Q,theta):
    # Entre  com theta em graus
    T = calc_T(theta)
    Qbar = np.linalg.inv(T).dot(Q.dot(np.transpose(np.linalg.inv(T))))
    return Qbar

def calc_Q(EL,ET,GLT,nuLT):
    nuTL = nuLT*ET/EL# Coeficiente de Poisson secundário
    Q11 = EL/(1-nuLT*nuTL)
    Q22 = ET/(1-nuLT*nuTL)
    Q66 = GLT
    Q16 = 0
    Q26 = 0
    Q12 = nuTL*EL/(1-nuLT*nuTL)
    Q = np.array([[Q11, Q12, Q16],
         [Q12, Q22, Q26],
         [Q16, Q26, Q66]])
    return  Q