#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:34:02 2019

@author: cristina
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_length(spectro_K, vv, N_atoms, N_omega, N):
    
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