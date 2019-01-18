#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 17:24:47 2018

@author: cristina
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
plt.rcParams.update({'font.size': 13})

def map2D_3D(e, e_up, e_down, e_x, e_13, e_24, z, z_up, z_down, z_x, z_13, z_24, titulo, N_x, N_y, N_omega, vv):
    
    #2D plots
    plt.figure(7)
    plt.imshow(z, cmap = plt.cm.jet)
    plt.colorbar(orientation='horizontal')
    plt.title('FM E = %f meV')
    plt.savefig('results/2D.pdf')
    
#    plt.figure(6)
#    plt.imshow(z_up, cmap = plt.cm.jet)
#    plt.colorbar(orientation='horizontal')
#    plt.title('Spin up' %titulo)
#    plt.savefig('results/2D_up.pdf')
#    
#    plt.figure(7)
#    plt.imshow(z_down, cmap = plt.cm.jet)
#    plt.colorbar(orientation='horizontal')
#    plt.title('Spin down' %titulo)
#    plt.savefig('results/2D_down.pdf')
#    
    plt.figure(8)
    plt.imshow(z_x, cmap = plt.cm.jet)
    plt.colorbar(orientation='horizontal')
    plt.title('Spin x')
    plt.savefig('results/2D_x.pdf')
    
    plt.figure(9)
    plt.imshow(z_up - z_down, cmap = plt.cm.jet)
    plt.colorbar(orientation='horizontal')
    plt.title('Spin z')
    plt.savefig('results/2D_z.pdf')
    
    plt.figure(10)
    plt.imshow(z_13 - z_24, cmap = plt.cm.jet)
    plt.colorbar(orientation='horizontal')
    plt.title('Spin 13 - 24')
    plt.savefig('results/2D_13.pdf')
    
    
    #################
    #################
    
    plt.figure(11)
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
    #plt.title('Atom index vs spectro')
    plt.colorbar()
    plt.savefig('results/map.pdf')

#    plt.figure(9)
#    plt.imshow(e_up, aspect='auto', cmap = plt.cm.jet)
#    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
#    plt.xticks(ticks, ticklabels)
#    plt.yticks(ticks2)
#    plt.xlabel('Energy (meV)')
#    plt.ylabel('atom index')
#    plt.title('Spin up')
#    plt.savefig('results/map_up.pdf')
#    plt.colorbar()
#    
#    plt.figure(10)
#    plt.imshow(e_down, aspect='auto', cmap = plt.cm.jet)
#    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
#    plt.xticks(ticks, ticklabels)
#    plt.yticks(ticks2)
#    plt.xlabel('Energy (meV)')
#    plt.ylabel('atom index')
#    plt.title('Spin down')
#    plt.savefig('results/map_down.pdf')
#    plt.colorbar()
    
    plt.figure(12)
    plt.imshow(e_up - e_down, aspect='auto', cmap = plt.cm.jet)
    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
    plt.xticks(ticks, ticklabels)
    plt.yticks(ticks2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('atom index')
    plt.title('Spin z')
    plt.savefig('results/map_z.pdf')
    plt.colorbar()
    
    plt.figure(13)
    plt.imshow(e_13 - e_24, aspect='auto', cmap = plt.cm.jet)
    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
    plt.xticks(ticks, ticklabels)
    plt.yticks(ticks2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('atom index')
    plt.title('Spin 13 - 24')
    plt.savefig('results/map_13.pdf')
    plt.colorbar()
    
    plt.figure(14)
    plt.imshow(e_x, aspect='auto', cmap = plt.cm.jet)
    ticks2 = np.linspace(0,N_x-1,5, dtype = 'int')
    plt.xticks(ticks, ticklabels)
    plt.yticks(ticks2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('atom index')
    plt.title('Spin x')
    plt.savefig('results/map_x.pdf')
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
