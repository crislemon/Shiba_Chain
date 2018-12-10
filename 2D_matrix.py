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
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


import time

pi=np.pi
d = 1.0 #distance between sites
N_atoms = 3 #number of atoms
borde = 2
ancho = 3
alpha = 3.0 #SOC
state = 'FM' #spin state
k_F = 0.5
U = 5500./27211.6#%potential scatt
U = 0
j = 1800./27211.6
DOS = 1.0
s = 5.0/2.0
delta = 0.75/27211.6 #SC gap

################################################# We solve Dyson's equation

import Shiba_Chain2D as sc2
t1=time.time()
(gg , N_x, N_y, N_omega , vv, Self, Go) = sc2.Shiba_Chain2(d, N_atoms, state, alpha, borde, ancho, k_F, U,
j, DOS, s, delta)
t2 = time.time()
 
print('The program is finished after', t2 - t1)

##################################################

spectro = np.zeros([N_y, N_x, N_omega], dtype= 'float')
spectro_up = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin up spectrum
spectro_down = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin down spectrum
spectro_uphole = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin up_hole spectrum
spectro_downhole = np.zeros([N_y, N_x, N_omega], dtype= 'float')#spin down_hole spectrum

for i_atom in range(N_y):
    for j_atom in range(N_x):
         I = i_atom*N_x + j_atom

         for i_omega in range(N_omega):
             
             tr = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
             spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(2*pi)
             
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
import plot_espectro as spect
(titulo, ndexes, i) = spect.espectro(spectro, row, vv, borde)


#####
"Plot the chain profile"
import plot_profile as pf
pf.profile(N_x, titulo, spectro, row, ndexes, i)


######
"Plot the spectrum for all Nambu operators"
plt.figure(3)
###cambiado!!!!!!
spectro_up_1 = spectro_up[row, borde, :]
spectro_down_1 = spectro_down[row, borde, :]
spectro_uphole_1 = spectro_uphole[row, borde, :]
spectro_downhole_1 = spectro_downhole[row, borde, :]

plt.plot(vv, spectro_up_1,label='Nambu 1',linewidth=0.8)
plt.plot(vv, spectro_down_1,label='Nambu 2',linewidth=0.8)
plt.plot(vv, spectro_uphole_1,label='Nambu 3',linewidth=0.8)
plt.plot(vv, spectro_downhole_1,label='Nambu 4',linewidth=0.8)
plt.legend()
plt.xlabel('mV')
plt.ylabel('PDOS')
plt.title('Nambu components')
plt.savefig('results/Nambu.pdf')


###
"2D and 3D maps"
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

###
"2D plots"
plt.figure(4)
plt.imshow(z, cmap = plt.cm.jet)
plt.colorbar()
plt.title('FM E = %f meV' %titulo)
plt.savefig('results/2D.pdf')

plt.figure(5)
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


###
"3D plot"
Y = list(range(N_y))
X = list(range (N_x))
X, Y = np.meshgrid(X, Y)
fig2 = plt.figure(6)
ax = fig2.add_subplot((111), projection='3d')
ax.plot_wireframe(X, Y, z)
plt.title('FM E = %f meV' %titulo)
ax.set_zlabel('PDOS')
plt.savefig('results/3D.pdf')



"plot Green's function"
import plot_Green as pG
#pG.plot_Green(gg, N_x, N_y, N_omega, row, borde, N_atoms, vv)





