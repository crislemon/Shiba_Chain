#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:32:08 2018

@author: cristina
"""

import numpy as np
import cmath as cm
import scipy.spatial.distance

# Functions
pi = np.pi
sin = np.sin
cos = np.cos
sqrt = np.sqrt
exp = np.exp


def Free_Green(N_x, N_y, lomega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic):

    G = np.zeros([N_x * N_y * 4, N_x * N_y * 4], dtype=complex)
    g = np.zeros([N_x * N_y, N_x * N_y, 4, 4], dtype=complex)
    
    omega = lomega + 1j * Damping

    
    #######calculo de distancias
    i = np.arange(N_y)
    j = np.arange(N_x)
    I, J= np.meshgrid(j, i)
    ii = np.reshape(I, ((N_x*N_y), ))
    jj = np.reshape(J, ((N_x*N_y), ))
    ij = zip(jj,ii)
    ij = list(ij)
    IJ = np.array(ij, dtype = 'double')
    rr = scipy.spatial.distance.cdist(IJ, IJ, metric='euclidean')*a_interatomic#distance between sites
    rr[np.where(rr == 0)] = 100   # avoid 1 / 0 errors !!
    
    
    # Non diagonal in atom
    SS = sqrt(Delta**2 - omega**2)
    xi = Fermi_k / (mass_eff * SS)

    for i_y in range(N_y):
        for i_x in range(N_x):
            g_i = (i_y)*(N_x) + (i_x)

        
            for j_y in range(N_y):
                for j_x in range(N_x):
                    
                    g_j = (j_y)*(N_x) + j_x
                    factor = - pi * DOS_o * exp(-rr[g_i, g_j]/ xi) / (SS * Fermi_k * rr[g_i, g_j])
                    
                
                
                    g[g_i, g_j, 0, 0] = ( omega * sin(Fermi_k * rr[g_i, g_j]) + SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                    g[g_i, g_j, 1, 1] = ( omega * sin(Fermi_k * rr[g_i, g_j]) - SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                    g[g_i, g_j, 2, 2] = ( omega * sin(Fermi_k * rr[g_i, g_j]) + SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                    g[g_i, g_j, 3, 3] = ( omega * sin(Fermi_k * rr[g_i, g_j]) - SS * cos(Fermi_k * rr[g_i, g_j]) )* factor
                    
                    g[g_i, g_j, 0, 3] = - Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
                    g[g_i, g_j, 1, 2] = Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
                    g[g_i, g_j, 2, 1] = Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
                    g[g_i, g_j, 3, 0] = - Delta * sin(Fermi_k * rr[g_i, g_j]) * factor
    

    # Diagonal in atom
    omega = lomega + 1j * Damping
    SS = sqrt(Delta**2 - omega**2)
    factor_diag = - pi * DOS_o / SS

    for g_i in range(N_x * N_y):
        
            
            g[g_i, g_i, 0, 0] = omega * factor_diag
            g[g_i, g_i, 1, 1] = omega * factor_diag
            g[g_i, g_i, 2, 2] = omega * factor_diag
            g[g_i, g_i, 3, 3] = omega * factor_diag
            g[g_i, g_i, 0, 3] = -Delta * factor_diag
            g[g_i, g_i, 1, 2] = Delta * factor_diag
            g[g_i, g_i, 2, 1] = Delta * factor_diag
            g[g_i, g_i, 3, 0] = -Delta * factor_diag
            

    for i in range(N_x*N_y):
        for j in range(N_x*N_y):
            for t_i in range(4):
                for t_j in range(4):
                    G[(i) * 4 + t_i, (j) * 4 + t_j] = g[i, j, t_i, t_j]
    
    
    
    return (G, g)
