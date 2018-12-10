#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 13:00:59 2018

@author: cristina
"""

import numpy as np
import scipy.spatial.distance


N_x = 4
N_y = 3
nstep = 1.0

#factor = 33
Damping = 0.018/27211.6 #Dynes damping
omega = 0.76/27211.6
omega = omega + 1j * Damping
a_interatomic=nstep*3.36/0.529

Delta= 0.75/27211.6 #SC gap
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

Free_Green = np.zeros([4 * N_y * N_x, 4* N_y * N_x], dtype= 'complex')



Free_Green[t_j * 4 + 0, t_i * 4 + 0] = omega * np.sin(Fermi_k * rr[t_j,t_i]) +\
     SS * np.cos(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]    
Free_Green[t_j * 4 + 1, t_i * 4 + 1] = omega * np.sin(Fermi_k * rr[t_j,t_i]) -\
     SS * np.cos(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]
Free_Green[t_j * 4 + 2, t_i * 4 + 2] = omega * np.sin(Fermi_k * rr[t_j,t_i]) +\
     SS * np.cos(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]
Free_Green[t_j * 4 + 3, t_i * 4 + 3] = omega * np.sin(Fermi_k * rr[t_j,t_i]) -\
     SS * np.cos(Fermi_k * rr[t_j,t_i]) * factor[t_j,t_i]
     
     
SS = np.sqrt(Delta**2 - omega**2)
factor_diag = - np.pi * DOS_o / SS

Free_Green[t * 4 + 0, t * 4 + 0] = omega * factor_diag 
Free_Green[t * 4 + 1, t * 4 + 1] = omega * factor_diag 
Free_Green[t * 4 + 2, t * 4 + 2] = omega * factor_diag
Free_Green[t * 4 + 3, t * 4 + 3] = omega * factor_diag
Free_Green[t * 4 + 0, t * 4 + 3] = -Delta * factor_diag
Free_Green[t * 4 + 1, t * 4 + 2] = Delta * factor_diag
Free_Green[t * 4 + 2, t * 4 + 1] = Delta * factor_diag
Free_Green[t * 4 + 3, t * 4 + 0] = -Delta * factor_diag



print(Free_Green)           