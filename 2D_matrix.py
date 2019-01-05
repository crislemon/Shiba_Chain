#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 17:53:46 2018

@author: cristina
"""

"this program makes all the plots"


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import time

pi=np.pi
d = 1.0 #distance between sites
N_atoms = 1 #number of atoms
borde = 2
ancho = 5
alpha = 2.5 #SOC
state = 'FM' #spin state
k_F = 0.183
U = -5500./27211.6#%potential scatt
#U = 0.0
#U = -1000.0/27211.6#%potential scatt
j = -1800./27211.6 #coupling
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2001
range_omega = 4

################################################# We solve Dyson's equation

import Shiba_Chain2D as sc2
t1=time.time()
(gg , N_x, N_y, N_omega , vv, Go, Self, Self2) = sc2.Shiba_Chain2(d, N_atoms, state, alpha, borde, ancho, 
k_F, U, j, DOS, s, delta, N_omega, range_omega)
t2 = time.time()
 
print('The program is finished after', t2 - t1)

##################################################

#####
"The spectrum is obtained all Nambu components"
spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')

spectro_up = np.zeros([N_y, N_x, N_omega], dtype= 'float')#Nambu 1 spectrum
spectro_down = np.zeros([N_y, N_x, N_omega], dtype= 'float')#Nambu 2 spectrum
spectro_uphole = np.zeros([N_y, N_x, N_omega], dtype= 'float')#Nambu 3 spectrum
spectro_downhole = np.zeros([N_y, N_x, N_omega], dtype= 'float')#Nambu 4 spectrum

spectro_spinup = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin up spectrum
spectro_spindown = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin down spectru

for i_atom in range(N_y):
    for j_atom in range(N_x):
         I = i_atom*N_x + j_atom

         for i_omega in range(N_omega):
             
             tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 2, I*4 + 2, N_omega - (i_omega+1)]
             spectro_spinup[i_atom , j_atom, i_omega]= - (tr.imag)/(2*pi)
             
             tr2 = gg[I*4 + 1, I*4 + 1, i_omega] + gg[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
             spectro_spindown[i_atom , j_atom, i_omega]= - (tr2.imag)/(2*pi)
             
             tr3 = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega] + gg[I*4 + 2, I*4 + 2, N_omega - (i_omega+1)] + gg[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
             spectro[i_atom , j_atom, i_omega]= - (tr3.imag)/(2*pi)
             
             trup = gg[I*4 + 0, I*4 + 0, i_omega]
             spectro_up[i_atom , j_atom, i_omega]= - (trup.imag)/(2*pi)
             
             trdown = gg[I*4 + 1, I*4 + 1, i_omega]
             spectro_down[i_atom , j_atom, i_omega]= - (trdown.imag)/(2*pi)
             
             trup_hole = gg[I*4 + 2, I*4 + 2, i_omega]
             spectro_uphole[i_atom , j_atom, i_omega]= - (trup_hole.imag)/(2*pi)
             
             trdown_hole = gg[I*4 + 3, I*4 + 3, i_omega]
             spectro_downhole[i_atom , j_atom, i_omega]= - (trdown_hole.imag)/(2*pi)
#####
"Plot the spectrum in the first atom"
row = int(N_y/2)
medio=int(N_x/2)
import plot_espectro as spect
#(titulo, ndexes, i) = spect.espectro(spectro, spectro_2, spectro_3 ,row, vv, borde)
(titulo, ndexes, i) = spect.espectro(spectro, spectro_spinup, spectro_spindown ,row, vv, borde)


#####
"Plot the chain profile"
import plot_profile as pf
pf.profile(N_x, titulo, spectro, row, ndexes, i)


######
"Plot the spectrum for all Nambu operators"
import plot_nambu as nb
nb.Nambu(spectro_up, spectro_down, spectro_uphole, spectro_downhole, vv, row, borde)


###
"creates 2D maps"

import maps as mp
(e, e_up, e_down, z, z_up, z_down) = mp.maps(N_y, N_x, 
spectro, ndexes, i, N_omega, row, borde, vv, spectro_spinup, spectro_spindown)

#z is the PDOS every where in the array for the energy 
#corresponding to the closest peaks to zero in the first atom

#e is the spectrum along the atomic chain


###
"2D and 3D plots"
import plot_2D3D as map2D
map2D.map2D_3D(e, e_up, e_down, z, z_up, z_down, titulo, N_x, N_y, N_omega, vv)


"plot Green's function"
import plot_Green as pG
#pG.plot_Green(gg, N_x, N_y, N_omega, row, borde, N_atoms, vv)





