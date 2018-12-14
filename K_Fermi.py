#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 12:25:18 2018

@author: cristina
"""

import Shiba_Chain2D as sc2
import numpy as np
import detect_peaks as dp
import matplotlib.pyplot as plt
import time

t1 = time.time()

pi=np.pi
d = 1.0 #distance between sites
N_atoms = 20 #number of atoms
borde = 3
ancho = 5
alpha = 5.5 #SOC
state = 'FM' #spin state
k_F = 0.183
U = 5500./27211.6#%potential scatt
U = 0.0
j = 1800./27211.6 #coupling
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2001


K_Fermi = [0.19, 0.20, 0.202 ,0.205, 0.208 ,0.21, 0.215 ,0.22, 0.2205, 0.225, 0.23]
picos = []
#picos = np.full((len(K_Fermi), 2), np.inf)

for k_i in range(len(K_Fermi)):
    (gg , N_x, N_y, N_omega , vv, Go, Self, GG) = sc2.Shiba_Chain2(d, N_atoms, state, alpha, borde, ancho, 
    K_Fermi[k_i], U, j, DOS, s, delta, N_omega)
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    
    for i_atom in range(N_x):
        for j_atom in range(N_y):
            I = i_atom + (j_atom)*N_x

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                spectro[j_atom , i_atom, i_omega]= - (tr.imag)/(2*pi)
             
    
    row = int(N_y/2)
    spectro_chain = spectro[row, borde, :]#spectrum in the first atom
    ndexes = dp.detect_peaks(spectro_chain)#find the peaks
    peaks = vv[ndexes]

    minpeak = min(abs(peaks))#find the minimum
    peaks2 = abs(peaks)
    #peaks = peaks.tolist()
    #peaks2 = peaks2.tolist()
    picos_index = np.where(peaks2 == minpeak)
    picos.append( peaks[picos_index] )
    #i=peaks2.index(minpeak)#the index of the peak closest to zero
    del gg


for xe, ye in zip(K_Fermi, picos):
    plt.scatter([xe] * len(ye), ye, c = 'b', marker = 'o')

plt.xlabel('K_Fermi')
plt.ylabel('Peak position (meV)')
plt.title('Central peak position vs K_F')
plt.xlim((0.185,0.235))

t2 = time.time()

print('the program is finished after', t2 - t1)

plt.savefig('K_fermi.pdf')



