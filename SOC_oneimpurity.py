#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 12:28:09 2018

@author: cristina
"""
#########
# one impurity with different SOC for spin along z and along x directions


import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2D as sc
import detect_peaks as dp

pi=np.pi

n=1.0

alpha=np.linspace(0, 5.0, 5)

##########
#out of plane spin (along z)
##########


Vpeak_z_plus = np.zeros(len(alpha))
Vpeak_z_minus = np.zeros(len(alpha))

N_atoms = 1
state = 'FM'
borde = 2
ancho = 3
k_f = 0.2
U = 0

for n_i in range(len(alpha)):
    
    #spin || z case
    (gg , N_x, N_y, N_omega , vv, Self, Go) = sc.Shiba_Chain2(n, N_atoms, state, alpha[n_i], borde, ancho, k_f, U)
    
    
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')

    for i_atom in range(N_y):
        for j_atom in range(N_x):
            I = i_atom*N_x + j_atom

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(2*pi)
             
            
    
    plt.figure(1)
    plt.style.use('seaborn-bright')
    row = int(N_y/2)
    plt.plot(vv, spectro[row, borde,:],linewidth=0.8, label = '%i' % n_i)
    ndexes = dp.detect_peaks(spectro[row, borde, :])
    peaks = vv[ndexes]
    plt.plot(peaks,spectro[row, borde, ndexes],'y*')
    plt.xlabel('mV')
    plt.ylabel('PDOS')
    plt.legend()
    plt.title('spin || z case')
    
    
    Shiba=vv[ndexes]
        
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in ndexes:
        if vv[i] >= 0:
            Shiba_plus.append(vv[i])
        elif vv[i] <= 0:
            Shiba_minus.append(vv[i])
    
    Vpeak_z_plus[n_i]= Shiba_plus[0]     
    Vpeak_z_minus[n_i]= Shiba_minus[-1]
   
 
##########
#inplane spin (along x)
##########

state = 'inplane'
Vpeak_x_plus = np.zeros(len(alpha))
Vpeak_x_minus = np.zeros(len(alpha))

for n_i in range(len(alpha)):
    
    #spin || x
    (gg , N_x, N_y, N_omega , vv, Self, Go) = sc.Shiba_Chain2(n, N_atoms, state, alpha[n_i], borde, ancho, k_f, U)
    
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')

    for i_atom in range(N_y):
        for j_atom in range(N_x):
            I = i_atom*N_x + j_atom

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(2*pi)
             
    
    plt.figure(2)
    plt.style.use('seaborn-bright')
    row = int(N_y/2)
    plt.plot(vv, spectro[row, borde,:],linewidth=0.8, label = '%i' % n_i)
    ndexes = dp.detect_peaks(spectro[row, borde,:])
    peaks = vv[ndexes]
    plt.plot(peaks,spectro[row, borde, ndexes],'y*')
    plt.xlabel('mV')
    plt.ylabel('PDOS')
    plt.title('spin || x case')
    plt.legend()
    
    
    Shiba=vv[ndexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in ndexes:
        if vv[i]>=0:
            Shiba_plus.append(vv[i])
        elif vv[i] <= 0:
            Shiba_minus.append(vv[i])

    Vpeak_x_plus[n_i]=Shiba_plus[0]
    Vpeak_x_minus[n_i]=Shiba_minus[-1]

  

#Delta= 1.0#meV
#    
plt.figure(3)
#plt.style.use('seaborn-pastel')
plt.plot(alpha,Vpeak_z_plus,'b.-',alpha,Vpeak_z_minus,'b.-',label = '|| z')
plt.plot(alpha,Vpeak_x_plus,'r.-',alpha,Vpeak_x_minus,'r.-', label = '|| x')
plt.show()
plt.legend()
plt.xlabel('alpha')
plt.ylabel('Shiba peak (meV)')
plt.title('SOC')
#
#
#plt.figure(4)
#plt.plot(alpha,Vpeak_FM1_plus/Delta,'r.-', label = 'z')
#plt.plot(alpha,Vpeak_AF_plus/Delta,'b.-', label= 'x')
#plt.xlabel('alpha')
#plt.ylabel('E/Delta')
##plt.ylim(0, 0.4)
#plt.legend()
