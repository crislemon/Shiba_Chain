#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:16:14 2018

@author: cristina
"""

import numpy as np
import scipy.spatial.distance


N_x = 5
N_y = 3
nstep = 1.0

#factor = 33
Damping = 0.018/27211.6 #Dynes damping
omega = 0.35
omega = omega + 1j * Damping
a_interatomic=nstep*3.36/0.529

Delta=0.75/27211.6 #SC gap
Delta = 1
DOS_o=1 #Normal phase DOS
Fermi_k = 0.2
mass_eff=1 #SC Band effective mass



i = range(N_y)
j = range(N_x)
   
I, J = np.meshgrid(j, i)
ii = np.reshape(I, ((N_x*N_y), ))
jj = np.reshape(J, ((N_x*N_y), ))


ij = zip(jj,ii)
ij = list(ij)
IJ = np.array(ij, dtype = 'double')
rr = scipy.spatial.distance.cdist(IJ, IJ, metric='euclidean')*a_interatomic
rr[np.where(rr == 0)] = 100   # avoid 1 / 0 errors !!

t = np.arange(N_x*N_y)
T, T2 = np.meshgrid(t, t)
t_i = np.reshape(T, (((N_x*N_y) ** 2), ))
t_j = np.reshape(T2, (((N_x*N_y) ** 2), ))


SS = np.sqrt(Delta**2 - omega**2)
xi = Fermi_k / (mass_eff * SS)
factor = - np.pi * DOS_o * np.exp(-rr / xi) / (SS * Fermi_k * rr)

Free_Green = np.zeros([N_y * N_x, N_y * N_x, 4, 4], dtype= 'complex')



Free_Green[t_j, t_i, 0, 0] = omega * np.sin(Fermi_k * rr[t_j,t_i]) +\
SS * np.cos(Fermi_k * rr[t_j,t_i])# * factor[t_j,t_i]    
Free_Green[t_j, t_i, 1, 1] = omega * np.sin(Fermi_k * rr[t_j,t_i]) -\
     SS * np.cos(Fermi_k * rr[t_j,t_i])# * factor[t_j,t_i]
Free_Green[t_j, t_i, 2, 2] = omega * np.sin(Fermi_k * rr[t_j,t_i]) +\
     SS * np.cos(Fermi_k * rr[t_j,t_i])# * factor[t_j,t_i]
Free_Green[t_j, t_i, 3, 3] = omega * np.sin(Fermi_k * rr[t_j,t_i]) -\
     SS * np.cos(Fermi_k * rr[t_j,t_i])# * factor[t_j,t_i]
     
#for i in range(4):
#    for j in range(4):
#        Free_Green[t_j, t_i, i, j] = factor[t_j,t_i] * Free_Green[t_j, t_i, i, j]
     
     
SS = np.sqrt(Delta**2 - omega**2)
factor_diag = - np.pi * DOS_o / SS
Free_Green[t, t, 0, 0] = omega# * factor_diag #omega * factor
Free_Green[t, t, 1, 1] = omega# * factor_diag #omega * factor
Free_Green[t, t, 2, 2] = omega# * factor_diag#omega * factor
Free_Green[t, t, 3, 3] = omega# * factor_diag#omega * factor
Free_Green[t, t, 0, 3] = -Delta# * factor_diag
Free_Green[t, t, 1, 2] = Delta# * factor_diag
Free_Green[t, t, 2, 1] = Delta# * factor_diag
Free_Green[t, t, 3, 0] = -Delta# * factor_diag

#for i in range(4):
#    for j in range(4):
#        Free_Green[t_j, t_i, i, j] = factor_diag * Free_Green[t_j, t_i, i, j]



print(Free_Green)           