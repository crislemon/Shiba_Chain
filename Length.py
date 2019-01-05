#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 16:35:21 2018

@author: cristina
"""

import Shiba_Chain2D as sc2
import numpy as np
import matplotlib.pyplot as plt
import time

t1 = time.time()

pi=np.pi
d = 1.0 #distance between sites
#N_atoms = 22 #number of atoms
borde = 3
ancho = 7
alpha = 2.5 #SOC
state = 'FM' #spin state
k_F = 0.183
U = -5500./27211.6#%potential scatt
U = 0.0
j = -1800./27211.6 #coupling
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2000
range_omega = 2.0


start = 7
end = 18
N = end - start
N_atoms=np.linspace(start, end , N, dtype = 'int')


spectro_K = np.zeros([N_omega, N])

for N_i in range(N):
    
    (gg , N_x, N_y, N_omega , vv, Go, Self, Self2) = sc2.Shiba_Chain2(d, N_atoms[N_i], state, alpha, borde, ancho, 
    k_F, U, j, DOS, s, delta, N_omega, range_omega)
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    
    for i_atom in range(N_x):
        for j_atom in range(N_y):
            I = i_atom + (j_atom)*N_x

            for i_omega in range(N_omega):
             
                #tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega] + gg[I*4 + 2, I*4 + 2, N_omega - (i_omega+1)] + gg[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
             
                spectro[j_atom , i_atom, i_omega]= - (tr.imag)/(4*pi)
                
    row = int(N_y/2)         
    spectro_K[:, N_i] = spectro[row,borde,:]
    
    
    del gg
    
    
plt.figure(1)
plt.imshow(spectro_K, aspect='auto', cmap = plt.cm.gnuplot)
ticks = np.linspace(0, N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.yticks(ticks, ticklabels)

ticks2 = np.linspace(0, N - 1, 3, dtype = 'int')
ticklabels2 = N_atoms[ticks2]
plt.xticks(ticks2, ticklabels2)

plt.ylabel('Energy (meV)')
plt.xlabel('Number of atoms')

plt.colorbar()
plt.savefig('N_atoms.pdf')

t2 = time.time()
print('the program is finished after', t2 - t1)

    
    
    
    