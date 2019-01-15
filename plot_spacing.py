#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:05:56 2019

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})

def plot_spacing(spectro_K, vv, d, N_omega, N):
    
    plt.figure(1)
    plt.imshow(spectro_K, aspect='auto', cmap = plt.cm.gnuplot)
    ticks = np.linspace(0, N_omega - 1, 3, dtype = 'int')
    ticklabels = vv[ticks]
    ticklabels = np.around(ticklabels, decimals=2)
    plt.yticks(ticks, ticklabels)
    
    ticks2 = np.linspace(0, N - 1, 3, dtype = 'int')
    ticklabels2 = d[ticks2]
    ticklabels2 = np.around(ticklabels2, decimals=2)
    plt.xticks(ticks2, ticklabels2)

    plt.ylabel('Energy (meV)')
    plt.xlabel('Distance between sites (a)')

    plt.colorbar()
    plt.savefig('spacing.pdf')