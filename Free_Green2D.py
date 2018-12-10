#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import cmath as cm

# Functions
pi = np.pi
sin = np.sin
cos = np.cos
sqrt = np.sqrt
exp = np.exp

complex_v = np.vectorize(complex)


def Free_Green(N_x, N_y, lomega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic):

    G = np.zeros([N_x * N_y * 4, N_x * N_y * 4], dtype=complex)
    omega = lomega + 1j * Damping

    # Non diagonal in atom
    i = np.arange(N_y)
    j = np.arange(N_x)
    
    I, J, N, M = np.meshgrid(i, j, i, j)
    ii = np.reshape(I, ((N_x*N_y) ** 2, 1))
    jj = np.reshape(J, ((N_x*N_y) ** 2, 1))
    nn = np.reshape(N, ((N_x*N_y) ** 2, 1))
    mm = np.reshape(M, ((N_x*N_y) ** 2, 1))
    

    rr = sqrt((ii - nn) ** 2 + (jj - mm) ** 2) * a_interatomic
    rr[np.where(rr == 0)] = 1   # avoid 1 / 0 errors !!

    
    SS = sqrt(Delta**2 - omega**2)
    #SS = complex_v(sqrt(SS.real ** 2), SS.imag)#impose Re(SS) > 0
    xi = Fermi_k / (mass_eff * SS)
    factor = - pi * DOS_o * exp(-rr / xi) / (SS * Fermi_k * rr)

    g_i = (ii - 1) * N_x + jj
    g_j = (nn - 1) * N_x + mm

    G[g_i * 4 + 0, g_j * 4 + 0] = omega * sin(Fermi_k * rr) +\
     SS * cos(Fermi_k * rr) * factor
    G[g_i * 4 + 1, g_j * 4 + 1] = omega * sin(Fermi_k * rr) -\
     SS * cos(Fermi_k * rr) * factor
    G[g_i * 4 + 2, g_j * 4 + 2] = omega * sin(Fermi_k * rr) +\
     SS * cos(Fermi_k * rr) * factor
    G[g_i * 4 + 3, g_j * 4 + 3] = omega * sin(Fermi_k * rr) -\
     SS * cos(Fermi_k * rr) * factor

    G[g_i * 4 + 0, g_j * 4 + 3] = - Delta * sin(Fermi_k * rr) * factor
    G[g_i * 4 + 1, g_j * 4 + 2] = Delta * sin(Fermi_k * rr) * factor
    G[g_i * 4 + 2, g_j * 4 + 1] = Delta * sin(Fermi_k * rr) * factor
    G[g_i * 4 + 3, g_j * 4 + 0] = - Delta * sin(Fermi_k * rr) * factor

    # Diagonal in atom
    omega = lomega + 1j * Damping
    SS = sqrt(Delta**2 - omega**2)
    factor = - pi * DOS_o / SS

    G[g_i * 4 + 0, g_i * 4 + 0] = omega * factor
    G[g_i * 4 + 1, g_i * 4 + 1] = omega * factor
    G[g_i * 4 + 2, g_i * 4 + 2] = omega * factor
    G[g_i * 4 + 3, g_i * 4 + 3] = omega * factor
    G[g_i * 4 + 0, g_i * 4 + 3] = -Delta * factor
    G[g_i * 4 + 1, g_i * 4 + 2] = Delta * factor
    G[g_i * 4 + 2, g_i * 4 + 1] = Delta * factor
    G[g_i * 4 + 3, g_i * 4 + 0] = -Delta * factor
    
    
    return (G)

