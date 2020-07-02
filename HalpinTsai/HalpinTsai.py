#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:28:19 2020

@author: eder
"""

def calc_Vf(wf,wm,rof,rom):
    wc=wf+wm
    Wf=wf/wc
    Wm=wm/wc
    roc=1/(Wf/rof+Wm/rom)
    Vf=roc/rof*Wf
    return Vf

def calc_EL(Vf,Ef,Em):
    Vm = 1 - Vf
    EL = Vf*Ef + Vm*Em
    return EL

def calc_EL_fibras_curtas(Vf,Ef,Em,l,d):
    Vm = 1 - Vf
    eta=(Ef/Em-1)/(Ef/Em+2*(l/d))
    EL = Em*(1+(2*l/d)*eta*Vf)/(1-eta*Vf)
    return EL

def calc_nuLT(Vf,NUf,NUm):
    Vm = 1 - Vf
    nuLT = Vf*NUf + Vm*NUm
    return nuLT

def calc_ET(Vf,Ef,Em,p):
    Vm = 1 - Vf
    if (p == 1):
        ET = 1/(Vf/Ef + Vm/Em)
    elif (p == 2):
        xi=2
        eta=(Ef/Em-1)/(Ef/Em+xi)
        ET = Em*(1+xi*eta*Vf)/(1-eta*Vf)
    return ET

def calc_GLT(Vf,Gf,Gm,p,xi=1):
    Vm = 1 - Vf
    if (p == 1):
        GLT = 1/(Vf/Gf + Vm/Gm)
    elif (p == 2):
        eta=(Gf/Gm-1)/(Gf/Gm+xi)
        GLT = Gm*(1+xi*eta*Vf)/(1-eta*Vf)
    return GLT