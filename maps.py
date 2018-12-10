#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 17:11:39 2018

@author: cristina
"""
import numpy as np

def maps(N_y, N_x, spectro, ndexes, i, N_omega, row, borde, vv):
    #N_x X N_y array
    z = np.zeros([N_y, N_x], dtype = float)
    spectra = spectro[:,:,ndexes[i]]


    for j_atom in range(N_y):
        for i_atom in range(N_x):
            z[j_atom,i_atom] = spectra[j_atom,i_atom]

    #array for omega
    e = np.zeros([N_x, N_omega], dtype = float)
    spectra2 = spectro[row,:,:]

    for i_omega in range(N_omega):
        for i_atom in range(N_x):
            e[i_atom,i_omega] = spectra2[i_atom,i_omega]

    ###
    "Save data"
    data_spectro = np.array([vv,spectro[row, borde, :]])
    data_3D = z
    np.savetxt('results/data_spectro.txt', data_spectro)
    np.savetxt('results/data_3D.txt', data_3D)
    
    return(e, z)