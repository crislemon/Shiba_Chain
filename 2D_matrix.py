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
N_atoms = 24 #number of atoms
borde = 2
ancho = 5
#alpha = 4.0 #SOC
alpha = 3.0 #SOC
state = 'FM' #spin state
k_F = 0.55
U = -5500./27211.6#%potential scatt
U = 0.0
#U = -3500.0/27211.6#%potential scatt
j = -1800./27211.6 #coupling
j = -3100./27211.6
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
N_omega = 2003
range_omega = 4

################################################# We solve Dyson's equation

import Shiba_Chain2D as sc2
t1=time.time()
(gg , N_x, N_y, N_omega , vv, Go, Self2) = sc2.Shiba_Chain2(d, N_atoms, state, alpha, borde, ancho, 
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

spectro_spinx = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin x spectrum

spectro_13 = np.zeros([N_y, N_x, N_omega], dtype= 'float')#13
spectro_31 = np.zeros([N_y, N_x, N_omega], dtype= 'float')#31

spectro_24 = np.zeros([N_y, N_x, N_omega], dtype= 'float')#13
spectro_42 = np.zeros([N_y, N_x, N_omega], dtype= 'float')#31

for i_atom in range(N_y):
    for j_atom in range(N_x):
         I = i_atom*N_x + j_atom

         for i_omega in range(N_omega):
             
             tr = gg[I*4 + 0, I*4 + 0, i_omega]# + gg[I*4 + 2, I*4 + 2, N_omega - (i_omega+1)]
             spectro_spinup[i_atom , j_atom, i_omega]= - (tr.imag)/(pi)
             
             tr2 = gg[I*4 + 1, I*4 + 1, i_omega]# + gg[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
             spectro_spindown[i_atom , j_atom, i_omega]= - (tr2.imag)/(pi)
             
             #tr3 = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega] + gg[I*4 + 2, I*4 + 2, N_omega - (i_omega+1)] + gg[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
             tr3 = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
             spectro[i_atom , j_atom, i_omega] = - (tr3.imag)/(pi)
             
             trup = gg[I*4 + 0, I*4 + 0, i_omega]
             spectro_up[i_atom , j_atom, i_omega]= - (trup.imag)/(pi)
             
             trdown = gg[I*4 + 1, I*4 + 1, i_omega]
             spectro_down[i_atom , j_atom, i_omega]= - (trdown.imag)/(pi)
             
             trup_hole = gg[I*4 + 2, I*4 + 2, i_omega]
             spectro_uphole[i_atom , j_atom, i_omega]= - (trup_hole.imag)/(pi)
             
             trdown_hole = gg[I*4 + 3, I*4 + 3, i_omega]
             spectro_downhole[i_atom , j_atom, i_omega]= - (trdown_hole.imag)/(pi)
             
             trx = gg[I*4 + 2, I*4 + 3, i_omega] + gg[I*4 + 3, I*4 + 2, i_omega]
             spectro_spinx[i_atom , j_atom, i_omega] = - (trx.imag)/(2*pi)
             
             tr_13 = gg[I*4 + 0, I*4 + 2, i_omega]
             spectro_13[i_atom , j_atom, i_omega] = - (tr_13.imag)/(pi)
             tr_31 = gg[I*4 + 2, I*4 + 0, i_omega]
             spectro_31[i_atom , j_atom, i_omega] = - (tr_31.imag)/(pi)
             
             tr_24 = gg[I*4 + 1, I*4 + 3, i_omega]
             spectro_24[i_atom , j_atom, i_omega] = - (tr_13.imag)/(pi)
             tr_42 = gg[I*4 + 3, I*4 + 1, i_omega]
             spectro_42[i_atom , j_atom, i_omega] = - (tr_31.imag)/(pi)
             
#####
"Plot the spectrum in the first atom"
row = int(N_y/2)
medio=int(N_x/2)
import plot_espectro as spect

(titulo, ndexes, i) = spect.espectro(spectro, spectro_spinup, spectro_spindown, spectro_13, spectro_31, spectro_24, spectro_42, row, vv, borde)


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
(e, e_up, e_down, e_x, e_13, e_24, z, z_up, z_down, z_x, z_13, z_24) = mp.maps(N_y, N_x, 
spectro, ndexes, i, N_omega, row, borde, vv, spectro_spinup, spectro_spindown, spectro_spinx, spectro_13, spectro_24)

#z is the PDOS every where in the array for the energy 
#corresponding to the closest peaks to zero in the first atom

#e is the spectrum along the atomic chain


###
"2D and 3D plots"
import plot_2D3D as map2D
map2D.map2D_3D(e, e_up, e_down, e_x, e_13, e_24, z, z_up, z_down, z_x, z_13, z_24, titulo, N_x, N_y, N_omega, vv)


"plot Green's function"
import plot_Green as pG
#pG.plot_Green(gg, N_x, N_y, N_omega, row, borde, N_atoms, vv)





