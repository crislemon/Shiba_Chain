#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:03:24 2019

@author: cristina
"""

import Shiba_Chain2D as sc2
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import time

t1 = time.time()

pi=np.pi
#d = 1.0 #distance between sites
borde = 2
ancho = 5
alpha = 3.5 #SOC
state = 'FM' #spin state
k_F = 0.21
U = -5500./27211.6#%potential scatt
U = 0.0
j = -1800./27211.6 #coupling
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2003
range_omega = 2.0
N_atoms = 2


start = 0.4
end = 2.5
N = 5
d = np.linspace(start, end , N)


spectro_K = np.zeros([N_omega, N])

for d_i in range(N):
    
    (gg , N_x, N_y, N_omega , vv, Go, Self2) = sc2.Shiba_Chain2(d[d_i], N_atoms, state, alpha, borde, ancho, 
    k_F, U, j, DOS, s, delta, N_omega, range_omega)
    spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    
    for i_atom in range(N_x):
        for j_atom in range(N_y):
            I = i_atom + (j_atom)*N_x

            for i_omega in range(N_omega):
             
                tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
             
                spectro[j_atom , i_atom, i_omega]= - (tr.imag)/(2*pi)
                
    row = int(N_y/2)         
    spectro_K[:, d_i] = spectro[row,borde,:]
    
    
    del gg
    
np.savetxt('spacing.txt', spectro_K)
np.savetxt('vv.txt', vv)
np.savetxt('d.txt', d)
    
import plot_spacing as pl
pl.plot_spacing(spectro_K, vv, d, N_omega, N)

t2 = time.time()
print('the program is finished after', t2 - t1)

    