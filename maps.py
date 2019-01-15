#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 17:11:39 2018

@author: cristina
"""
import numpy as np

def maps(N_y, N_x, spectro, ndexes, i, N_omega, row, borde, vv, spectro_spinup, spectro_spindown, spectro_spinx):
    #N_x X N_y array
    z = np.zeros([N_y, N_x], dtype = float)
    z_up = np.zeros([N_y, N_x], dtype = float)
    z_down = np.zeros([N_y, N_x], dtype = float)
    z_x = np.zeros([N_y, N_x], dtype = float)
                      
    spectra = spectro[:,:,ndexes[i]]
    spectra_up = spectro_spinup[:,:,ndexes[i]]
    spectra_down = spectro_spindown[:,:,ndexes[i]]
    spectra_x = spectro_spinx[:,:,ndexes[i]]
    

    for j_atom in range(N_y):
        for i_atom in range(N_x):
            z[j_atom,i_atom] = spectra[j_atom,i_atom]
            z_up[j_atom,i_atom] = spectra_up[j_atom,i_atom]
            z_down[j_atom,i_atom] = spectra_down[j_atom,i_atom]
            z_x[j_atom,i_atom] = spectra_x[j_atom,i_atom]

    #array for omega
    e = np.zeros([N_x, N_omega], dtype = float)
    e_up = np.zeros([N_x, N_omega], dtype = float)
    e_down = np.zeros([N_x, N_omega], dtype = float)
    e_x = np.zeros([N_x, N_omega], dtype = float)
    
    spectra2 = spectro[row,:,:]
    spectra2_up = spectro_spinup[row,:,:]
    spectra2_down = spectro_spindown[row,:,:]
    spectra2_x = spectro_spinx[row,:,:]
    

    for i_omega in range(N_omega):
        for i_atom in range(N_x):
            e[i_atom,i_omega] = spectra2[i_atom,i_omega]
            e_up[i_atom,i_omega] = spectra2_up[i_atom,i_omega]
            e_down[i_atom,i_omega] = spectra2_down[i_atom,i_omega]
            e_x[i_atom,i_omega] = spectra2_x[i_atom,i_omega]
    ###
    "Save data"
    data_spectro = np.array([vv,spectro[row, borde, :]])
    data_3D = z
    np.savetxt('results/data_spectro.txt', data_spectro)
    np.savetxt('results/data_3D.txt', data_3D)
    
    return(e, e_up, e_down, e_x, z, z_up, z_down, z_x)