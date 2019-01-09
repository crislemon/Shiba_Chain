#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 16:35:21 2018

@author: cristina
"""

import Shiba_Chain2D as sc2
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import time

t1 = time.time()

pi=np.pi
d = 1.0 #distance between sites
borde = 1
ancho = 3
alpha = 0.0 #SOC
state = 'FM' #spin state
k_F = 0.4
U = -5500./27211.6#%potential scatt
U = 0.0
j = -1800./27211.6 #coupling
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2003
range_omega = 2.0


start = 1
end = 9
N = end - start + 1
N_atoms=np.linspace(start, end , N, dtype = 'int')


spectro_K = np.zeros([N_omega, N])

for N_i in range(N):
    
    (gg , N_x, N_y, N_omega , vv, Go, Self2) = sc2.Shiba_Chain2(d, N_atoms[N_i], state, alpha, borde, ancho, 
    k_F, U, j, DOS, s, delta, N_omega, range_omega)
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    
    for i_atom in range(N_x):
        for j_atom in range(N_y):
            I = i_atom + (j_atom)*N_x

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                #tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
                #tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega] + gg[I*4 + 2, I*4 + 2, N_omega - (i_omega+1)] + gg[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
             
                spectro[j_atom , i_atom, i_omega]= - (tr.imag)/(2*pi)
                
    row = int(N_y/2)         
    spectro_K[:, N_i] = spectro[row,borde,:]
    
    
    del gg
    
np.savetxt('length.txt', spectro_K)
np.savetxt('vv.txt', vv)
np.savetxt('N_atoms.txt', N_atoms)
    
import plot_Length as pl
pl.plot_length(spectro_K, vv, N_atoms, N_omega, N) 

t2 = time.time()
print('the program is finished after', t2 - t1)

    
    
    
    