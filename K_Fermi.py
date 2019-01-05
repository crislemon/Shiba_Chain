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
N_atoms = 22 #number of atoms
borde = 2
ancho = 5
alpha = 2.5 #SOC
state = 'FM' #spin state
#k_F = 0.183
U = -5500./27211.6#%potential scatt
U = 0.0
j = -1800./27211.6 #coupling
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2001
range_omega = 1.0

N_k = 6


k_Fermi = np.linspace(0.15,0.2, N_k)
picos = []
spectro_K = np.zeros([N_omega, N_k])

for k_i in range(len(k_Fermi)):
    (gg , N_x, N_y, N_omega , vv, Go, Self, Go2, Self2) = sc2.Shiba_Chain2(d, N_atoms, state, alpha, borde, ancho, 
    k_Fermi[k_i], U, j, DOS, s, delta, N_omega, range_omega)
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    
    for i_atom in range(N_y):
        for j_atom in range(N_x):
            I = i_atom*N_x + j_atom

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(2*pi)
                
    row = int(N_y/2)         
    spectro_K[:, k_i] = spectro[row,borde,:]
    
#    spectro_chain = spectro[row, borde, :]#spectrum in the first atom
#    ndexes = dp.detect_peaks(spectro_chain)#find the peaks
#    peaks = vv[ndexes]
#
#    minpeak = min(abs(peaks))#find the minimum
#    peaks2 = abs(peaks)
#    
#    picos_index = np.where(peaks2 == minpeak)
#    picos.append( peaks[picos_index] )
#    
    del gg


#for xe, ye in zip(k_Fermi, picos):
#    plt.scatter([xe] * len(ye), ye, c = 'b', marker = 'o')
#
#plt.figure(1)
#plt.xlabel('K_Fermi')
#plt.ylabel('Peak position (meV)')
#plt.title('Central peak position vs K_F')
##plt.xlim((0.17,0.215))
#plt.savefig('results/K_fermi.pdf')

plt.figure(1)
plt.imshow(spectro_K, aspect='auto', cmap = plt.cm.gnuplot)
ticks = np.linspace(0, N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.yticks(ticks, ticklabels)

ticks2 = np.linspace(0, len(k_Fermi) - 1, 3, dtype = 'int')
ticklabels2 = k_Fermi[ticks2]
ticklabels2 = np.around(ticklabels2, decimals=3)
plt.xticks(ticks2, ticklabels2)

plt.ylabel('Energy (meV)')
plt.xlabel('k_f (a.u.)')
plt.title('E vs k_F')
plt.colorbar()
plt.savefig('results/k_fermi2_dimer.pdf')

t2 = time.time()
print('the program is finished after', t2 - t1)





