#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:13:10 2018

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
    omega = lomega + 1j * Damping

    # Non diagonal in atom
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

    t = np.arange(N_x*N_y)
    T, T2 = np.meshgrid(t, t)
    t_i = np.reshape(T, (((N_x*N_y) ** 2), ))
    t_j = np.reshape(T2, (((N_x*N_y) ** 2), ))
    
    SS = sqrt(Delta**2 - omega**2)
    xi = Fermi_k / (mass_eff * SS)
    factor = - pi * DOS_o * exp(-rr/ xi) / (SS * Fermi_k * rr)
    

    G[t_j * 4 + 0, t_i * 4 + 0] = ( omega * sin(Fermi_k * rr[t_j,t_i]) + SS * cos(Fermi_k * rr[t_j,t_i]) )* factor[t_j,t_i]
    G[t_j * 4 + 1, t_i * 4 + 1] = ( omega * sin(Fermi_k * rr[t_j,t_i]) - SS * cos(Fermi_k * rr[t_j,t_i]) )* factor[t_j,t_i]
    G[t_j * 4 + 2, t_i * 4 + 2] = ( omega * sin(Fermi_k * rr[t_j,t_i]) + SS * cos(Fermi_k * rr[t_j,t_i]) )* factor[t_j,t_i]
    G[t_j * 4 + 3, t_i * 4 + 3] = ( omega * sin(Fermi_k * rr[t_j,t_i]) - SS * cos(Fermi_k * rr[t_j,t_i]) )* factor[t_j,t_i]

    G[t_j * 4 + 0, t_i * 4 + 3] = - Delta * sin(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]
    G[t_j * 4 + 1, t_i * 4 + 2] = Delta * sin(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]
    G[t_j * 4 + 2, t_i * 4 + 1] = Delta * sin(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]
    G[t_j * 4 + 3, t_i * 4 + 0] = - Delta * sin(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]

    # Diagonal in atom
    omega = lomega + 1j * Damping
    SS = sqrt(Delta**2 - omega**2)
    factor_diag = - pi * DOS_o / SS

    G[t * 4 + 0, t * 4 + 0] = omega * factor_diag
    G[t * 4 + 1, t * 4 + 1] = omega * factor_diag
    G[t * 4 + 2, t * 4 + 2] = omega * factor_diag
    G[t * 4 + 3, t * 4 + 3] = omega * factor_diag
    G[t * 4 + 0, t * 4 + 3] = -Delta * factor_diag
    G[t * 4 + 1, t * 4 + 2] = Delta * factor_diag
    G[t * 4 + 2, t * 4 + 1] = Delta * factor_diag
    G[t * 4 + 3, t * 4 + 0] = -Delta * factor_diag
    
    
    return (G)

