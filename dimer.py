#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 15:24:27 2019

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2D as sc
plt.rcParams.update({'font.size': 13})

pi=np.pi
d_ini = 0.4
d_final = 2.5
N = np.linspace(d_final, d_ini, 101)
#d = 1.0 #distance between sites
N_atoms = 2 #number of atoms
borde = 1
ancho = 3
#alpha = 4.0 #SOC
alpha = 3.0 #SOC
k_F = 0.4
U = -5500./27211.6#%potential scatt
U = 0.0
#U = -3500.0/27211.6#%potential scatt
j = -4000./27211.6 #coupling
#j = -3100./27211.6
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2003
range_omega = 3
row = int(ancho/2)

spectro_AF = np.zeros([len(N), N_omega ])
spectro_FM = np.zeros([len(N), N_omega ])

""""""""""""""""""""""""
"Anti-ferromagnetic spin1 = 0 spin2 = pi"
""""""""""""""""""""""""

state = 'AF' #spin state

for n_i in range(len(N)):
    
    (gg , N_x, N_y, N_omega , vv, Go, Self2) = sc.Shiba_Chain2(N[n_i], N_atoms, state, alpha, borde, ancho, 
k_F, U, j, DOS, s, delta, N_omega, range_omega)
    
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    for i_atom in range(N_y):
        for j_atom in range(N_x):
            I = i_atom*N_x + j_atom

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(pi)
                
    spectro_AF[n_i, :] = spectro[row, borde, :]
             
             
#######################
""""""""""""""""""""""""
"Ferromagnetic spin1 = 0 spin2 = 0"
""""""""""""""""""""""""

state = 'FM' #spin state

for n_i in range(len(N)):
    
    (gg , N_x, N_y, N_omega , vv, Go, Self2) = sc.Shiba_Chain2(N[n_i], N_atoms, state, alpha, borde, ancho, 
k_F, U, j, DOS, s, delta, N_omega, range_omega)
    
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    for i_atom in range(N_y):
        for j_atom in range(N_x):
            I = i_atom*N_x + j_atom

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(pi)
                
    spectro_FM[n_i, :] = spectro[row, borde, :]             

np.savetxt('results/vv.txt', vv)
             
#AF            
plt.figure(1)
plt.imshow(spectro_AF, aspect='auto', cmap = plt.cm.gnuplot)

ticks = np.linspace(0, N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.xticks(ticks, ticklabels)

ticks2 = np.linspace(0, len(N) - 1, 3, dtype = 'int')
ticklabels2 = N[ticks2]
ticklabels2 = np.around(ticklabels2, decimals=2)
plt.yticks(ticks2, ticklabels2)

plt.xlabel('Energy (meV)')
plt.ylabel('d (a)')
#plt.title('AF dimer U = 5.5 eV')
plt.colorbar()
plt.savefig('results/AF_dimer.pdf')
np.savetxt('results/AF.txt', spectro_AF)

#FM
plt.figure(2)
plt.imshow(spectro_FM, aspect='auto', cmap = plt.cm.gnuplot)

ticks = np.linspace(0,N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.xticks(ticks, ticklabels)

ticks2 = np.linspace(0, len(N) - 1, 3, dtype = 'int')
ticklabels2 = N[ticks2]
ticklabels2 = np.around(ticklabels2, decimals=2)
plt.yticks(ticks2, ticklabels2)

plt.xlabel('Energy (meV)')
plt.ylabel('d (a)')
#plt.title('FM up dimer U = 5.5 eV')
plt.colorbar()
plt.savefig('results/FM_dimer.pdf')
np.savetxt('results/FM.txt', spectro_FM)          
             
             
             
             