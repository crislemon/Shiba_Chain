#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 17:24:47 2018

@author: cristina
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

def map2D_3D(e, e_up, e_down, z, z_up, z_down, titulo, N_x, N_y, N_omega, vv):
    
    #2D plots
    plt.figure(5)
    plt.imshow(z, cmap = plt.cm.jet)
    plt.colorbar()
    plt.title('FM E = %f meV' %titulo)
    plt.savefig('results/2D.pdf')
    
    plt.figure(6)
    plt.imshow(z_up, cmap = plt.cm.jet)
    plt.colorbar()
    plt.title('Spin up' %titulo)
    plt.savefig('results/2D_N1.pdf')
    
    plt.figure(7)
    plt.imshow(z_down, cmap = plt.cm.jet)
    plt.colorbar()
    plt.title('Spin down' %titulo)
    plt.savefig('results/2D_N2.pdf')

    plt.figure(8)
    plt.imshow(e, aspect='auto', cmap = plt.cm.jet)
    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
    ticks = np.linspace(0,N_omega-1,3, dtype = 'int')
    ticklabels = vv[ticks]
    
    for i in range(len(ticklabels)):
        ticklabels[i] = format(ticklabels[i], ".3g")
        
    plt.xticks(ticks, ticklabels)
    plt.yticks(ticks2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('atom index')
    plt.title('Atom index vs spectro')
    plt.colorbar()
    plt.savefig('results/map.pdf')

    
    plt.figure(9)
    plt.imshow(e_up, aspect='auto', cmap = plt.cm.jet)
    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
    plt.xticks(ticks, ticklabels)
    plt.yticks(ticks2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('atom index')
    plt.title('Spin up')
    plt.savefig('results/map_up.pdf')
    plt.colorbar()
    
    plt.figure(10)
    plt.imshow(e_down, aspect='auto', cmap = plt.cm.jet)
    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
    plt.xticks(ticks, ticklabels)
    plt.yticks(ticks2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('atom index')
    plt.title('Spin down')
    plt.savefig('results/map_down.pdf')
    plt.colorbar()

    
    #3D plot
    Y = list(range(N_y))
    X = list(range (N_x))
    X, Y = np.meshgrid(X, Y)
    fig2 = plt.figure(11)
    ax = fig2.add_subplot((111), projection='3d')
    ax.plot_wireframe(X, Y, z)
    plt.title('FM E = %f meV' %titulo)
    ax.set_zlabel('PDOS')
    plt.savefig('results/3D.pdf')
