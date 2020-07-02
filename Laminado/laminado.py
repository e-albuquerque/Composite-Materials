# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
import mpl_toolkits.mplot3d as mp
import mpl_toolkits.mplot3d.art3d as ar

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


def plota_deslocamentos(XYZ,uvw,fator):
    # Cria a variável escala que será aplicada sobre os deslocamentos para produzir a malha deformada
    xmax=np.max(XYZ[:,0])
    ymax=np.max(XYZ[:,1])
    xmin=np.min(XYZ[:,0])
    ymin=np.min(XYZ[:,1])
    distmaxima=np.sqrt((xmax-xmin)**2+(ymax-ymin)**2) # distância máxima entre 2 nós da malha
    dt=np.sqrt(uvw[:,0]**2+uvw[:,1]**2+uvw[:,2]**2) # deslocamentos totais dos nós
    deslocmaximo=np.max(np.abs(dt)) # deslocamento total máximo
    print("deslocamento máximo = ",deslocmaximo)
    escala=fator*distmaxima/deslocmaximo
    
    # Cria os triângulos que são armazenados na matriz elem
    nostriang = tri.Triangulation(XYZ[:,0], XYZ[:,1])
    triangulos=nostriang.triangles
    pc=[]
    poligono=np.zeros((3,3))
    nelem=triangulos.shape[0]
    x=np.zeros(3)
    y=np.zeros(3)
    z=np.zeros(3)
    dc=np.zeros(3)
    zc=np.zeros(nelem)
    XYZdef=XYZ+escala*uvw # Coordenada dos nos da malha deformada
    for elem in range(nelem):
        no1=triangulos[elem,0] # no1 do elemento elem
        no2=triangulos[elem,1]
        no3=triangulos[elem,2]
        x[0]=XYZdef[no1,0]
        y[0]=XYZdef[no1,1]
        z[0]=XYZdef[no1,2]
        dc[0]=dt[no1]
        x[1]=XYZdef[no2,0]
        y[1]=XYZdef[no2,1]
        z[1]=XYZdef[no2,2]
        dc[1]=dt[no2]	
        x[2]=XYZdef[no3,0]
        y[2]=XYZdef[no3,1]
        z[2]=XYZdef[no3,2]
        dc[2]=dt[no3]	
        zc[elem]=sum(dc)/3
        poligono=np.array([[x[0], y[0], z[0]],
        [x[1], y[1], z[1]],
        [x[2], y[2], z[2]]])
        pc.append(poligono)
        # plota a malha não deformada (2D)
    plt.figure()
    plt.triplot(XYZ[:,0], XYZ[:,1], triangulos, color=(0.0,0.,0.),linewidth=0.4,linestyle="--")
    
    # plota a malha deformada (3D)
    fig = plt.figure()
    ax = mp.Axes3D(fig)
    q = ar.Poly3DCollection(pc, linewidths=1)
    ax.add_collection3d(q)
    m = cm.ScalarMappable(cmap=cm.jet)
    b=m.to_rgba(zc)
    q.set_facecolor(b)
    m.set_array([np.min(zc),np.max(zc)])
    m.set_clim(vmin=np.min(zc),vmax=np.max(zc))
    plt.colorbar(m, orientation="vertical",shrink=0.9)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(np.min(XYZ[:,0])-escala*deslocmaximo*1.1,np.max(XYZ[:,0])+escala*deslocmaximo*1.1)
    ax.set_ylim(np.min(XYZ[:,1])-escala*deslocmaximo*1.1,np.max(XYZ[:,1])+escala*deslocmaximo*1.1)
    ax.set_zlim(np.min(XYZ[:,0])-escala*deslocmaximo*1.1,np.max(XYZ[:,0])+escala*deslocmaximo*1.1)
    ax.view_init(elev=18., azim=43.)
    ax.set_title("Deslocamentos totais")
    

def gera_NOSuvw(N,epgamma):
    # gera a matriz NOS e os vetores de deslocamento u, v e w.
    epsilonx=epgamma[0]
    epsilony=epgamma[1]
    gammaxy=epgamma[2]
    kappax=epgamma[3]
    kappay=epgamma[4]
    kappaxy=epgamma[5]
    
    # inicializa a matriz NOS
    # A matriz nós terá 2 colunas. Na primeira coluna
    # armazernará as coordenadas x dos nós e na segunda coluna
    # as coordenadas y. O número de linhas é N^2
    XYZ=np.zeros((N**2,3))
    # Inicializa os vetores deslocamentos
    uvw=np.zeros((N**2,3)) # deslocamentos na direção x
    # coordenadas mínimas e máximas da lâmina
    xmin=-1
    xmax=1
    ymin=-1
    ymax=1
    # distância entre oois pontos da malha
    deltax=(xmax-xmin)/(N-1) # distância na direção x
    deltay=(ymax-ymin)/(N-1) # distância na direção y
    
    # Calcula as deformações epsilonx, epsilony e gammaxy
    ino=0 # Inicializa o contador de nós
    # Atribui os valores aos elementos da matrz NOS e dos vetores u e v
    xx=xmin # inicializa x
    for i in range(N):
        yy=ymin; # inicializa y
        for j in range(N):
            XYZ[ino,0]=xx
            XYZ[ino,1]=yy
            XYZ[ino,2]=0
            uvw[ino,1]=epsilony*yy+gammaxy*xx/2
            uvw[ino,0]=epsilonx*xx+gammaxy*yy/2
            uvw[ino,2]=kappax*xx**2/2+kappay*yy**2/2+kappaxy*xx*yy				
            yy=yy+deltay # atualiza y
            ino=ino+1 # atualiza o contador de nos
        xx=xx+deltax # atualiza x
    return XYZ,uvw    



def plota_tensoes(epgamma,Material):
    n_laminae=Material.shape[0]
    thickness=np.sum(Material,0)[4]
    hant=-thickness/2
    vetz=np.zeros(2*n_laminae)
    vettens=np.zeros((2*n_laminae,3))
    vetdef=np.zeros((2*n_laminae,3))
    for i in range(n_laminae):
        EL=Material[i,0]
        ET=Material[i,1]
        GLT=Material[i,2]
        nuLT=Material[i,3]
        h = Material[i,4]   
        theta=Material[i,5]
#        nuTL=nuLT*ET/EL
        Q=calc_Q(EL,ET,GLT,nuLT)
        Qbar=calc_Qbar(Q,theta)
        sig1=Qbar.dot(epgamma[0:3])+Qbar.dot(epgamma[3:6])*hant
        sig2=Qbar.dot(epgamma[0:3])+Qbar.dot(epgamma[3:6])*(hant+h)
        def1=epgamma[0:3]+epgamma[3:6]*hant
        def2=epgamma[0:3]+epgamma[3:6]*(hant+h)
        vettens[2*i,:]=sig1
        vettens[2*i+1,:]=sig2
        vetdef[2*i,:]=def1
        vetdef[2*i+1,:]=def2
        vetz[2*i:2*i+2]=np.array([hant, hant+h])
        hant=hant+h
    plt.figure()
    plt.plot(vettens[:,0],vetz,linestyle="-",marker="x",label=r"$\sigma_x$")
    plt.plot(vettens[:,1],vetz,linestyle="--",marker="o",label=r"$\sigma_y$")
    plt.plot(vettens[:,2],vetz,linestyle="-.",marker="p",label=r"$\tau_{xy}$")
    plt.grid("on")
    plt.legend(loc="lower right",fancybox="true") # Create a legend of all the existing plots using their labels as names
    plt.xlabel("Tensões")
    plt.ylabel("$z$")
    plt.title("Tensões ao longo da espessura do laminado")
    
    plt.figure()
    plt.plot(vetdef[:,0],vetz,linestyle="-",marker="x",label=r"$\epsilon_x$")
    plt.plot(vetdef[:,1],vetz,linestyle="--",marker="o",label=r"$\epsilon_y$")
    plt.plot(vetdef[:,2],vetz,linestyle="-.",marker="p",label=r"$\gamma_{xy}$")
    plt.grid("on")
    plt.legend(loc="lower right",fancybox="true") # Create a legend of all the existing plots using their labels as names
    plt.xlabel("Deformações")
    plt.ylabel("$z$")
    plt.title("Deformações ao longo da espessura do laminado")