#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 10:02:15 2018

@author: cristina
"""

import numpy as np
import cmath as cm
import scipy.spatial.distance

N_x = 3
N_y = 5


g = np.zeros([N_x * N_y, N_x * N_y, 4, 4], dtype = complex)
G = np.zeros([N_x * N_y * 4, N_x * N_y * 4], dtype=complex)
tupla = np.zeros([N_x * N_y, N_x * N_y], dtype = tuple)


#######calculo de distancias
i = np.arange(N_y)
j = np.arange(N_x)
I, J= np.meshgrid(j, i)
ii = np.reshape(I, ((N_x*N_y), ))
jj = np.reshape(J, ((N_x*N_y), ))
ij = zip(jj,ii)
ij = list(ij)
IJ = np.array(ij, dtype = 'double')
rr = scipy.spatial.distance.cdist(IJ, IJ, metric='euclidean')#distance between sites

for i_y in range(N_y):
    for i_x in range(N_x):
        g_i = (i_y)*(N_x) + (i_x)

        
        for j_y in range(N_y):
            for j_x in range(N_x):
                g_j = (g_j) = (j_y)*(N_x) + j_x
                
                tupla[g_i, g_j] = (g_i, g_j)
                
                
                
                g[g_i, g_j, 0, 0] = rr[g_i,g_j]
                g[g_i, g_j, 1, 1] = 2
                g[g_i, g_j, 2, 2] = 3
                g[g_i, g_j, 3, 3] = 4
                
for i in range(N_x*N_y):
        for j in range(N_x*N_y):
            for t_i in range(4):
                for t_j in range(4):
                    G[(i) * 4 + t_i, (j) * 4 + t_j] = g[i, j, t_i, t_j]
    
 
print(g)

