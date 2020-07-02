#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:58:40 2020

@author: eder
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


def calc_ABD(Material):
    n_laminae=Material.shape[0]
    thickness=np.sum(Material,0)[4]
    hant=-thickness/2
    A=np.zeros((3,3))
    B=np.zeros((3,3))
    D=np.zeros((3,3))
    for i in range(n_laminae):
        EL=Material[i,0]
        ET=Material[i,1]
        GLT=Material[i,2]
        nuLT=Material[i,3]
        h = Material[i,4]
        theta=Material[i,5]
 #       nuTL=nuLT*ET/EL
        Q=calc_Q(EL,ET,GLT,nuLT)
        Qbar=calc_Qbar(Q,theta)
        A=A+Qbar*((h+hant)-hant)
        B=B+1/2*Qbar*((h+hant)**2-hant**2)
        D=D+1/3*Qbar*((h+hant)**3-hant**3)
        hant=hant+h
    AB= np.concatenate((A,B),axis=0)
    BD= np.concatenate((B,D),axis=0)
    ABD=np.concatenate((AB,BD),axis=1)
    return ABD


def calc_ABD2(Qs,thetas,hs):
    n_laminae=len(hs)
    thickness=np.sum(hs)
    hant=-thickness/2
    A=np.zeros((3,3))
    B=np.zeros((3,3))
    D=np.zeros((3,3))
    for i in range(n_laminae):
        Q=Qs[:,:,i]
        Qbar=calc_Qbar(Q,thetas[i])
        A=A+Qbar*((hs[i]+hant)-hant)
        B=B+1/2*Qbar*((hs[i]+hant)**2-hant**2)
        D=D+1/3*Qbar*((hs[i]+hant)**3-hant**3)
        hant=hant+hs[i]
    AB= np.concatenate((A,B),axis=0)
    BD= np.concatenate((B,D),axis=0)
    ABD=np.concatenate((AB,BD),axis=1)
    return ABD


def calc_TsaiWu(sigl,sigt,tault,xt,xc,yt,yc,s):
    a=sigl**2/(xt*xc)+sigt**2/(yt*yc)+tault**2/s**2-sigl*sigt/np.sqrt(xt*xc*yt*yc)
    b=(1/xt-1/xc)*sigl+(1/yt-1/yc)*sigt
    c=-1
    delta=b**2-4*a*c
    k1=(-b+np.sqrt(delta))/(2*a)
    k2=(-b-np.sqrt(delta))/(2*a)
    if(k1>0 and k2>0):
        if(k1<k2):
            k=k1
        else:
            k=k2
    elif(k1>0):
        k=k1
    elif(k2>0):
        k=k2
    else:
        print("Critério de falha negativo: $k1 e $k2")
    return k

def calc_TsaiHill(sigl,sigt,tault,xt,xc,yt,yc,s):
    if(sigl<0):
        xt=xc
    if(sigt<0):
        yt=yc
    indice=sigl**2/xt**2+sigt**2/yt**2-sigl*sigt/xt**2+tault**2/s**2
    k=1/np.sqrt(indice)
    return k

def calc_pressao(sigl,sigt,tault,xt,xc,yt,yc,s,k):
    a=k**2*(sigl**2/(xt*xc)+sigt**2/(yt*yc)+tault**2/s**2-sigl*sigt/np.sqrt(xt*xc*yt*yc))
    b=k*((1/xt-1/xc)*sigl+(1/yt-1/yc)*sigt)
    c=-1
    delta=b*2-4*a*c
    p1=(-b+np.sqrt(delta))/(2*a)
    p2=(-b-np.sqrt(delta))/(2*a)
    if(p1>0 and p2>0):
        if(p1<p2):
            p=p1
        else:
            p=p2
    elif(p1>0):
        p=p1
    elif(p2>0):
        p=p2
    else:
        print("As duas pressões são negativas")
    return p