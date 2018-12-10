#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 11:41:27 2018

@author: cristina
"""

import numpy as np
sin = np.sin
cos = np.cos
exp=np.exp
ui = complex(0.0, 1.0)


N_atoms = 2
N_x = 4
N_y = 3
lamda = 3
borde = 1
J = 1800./27211.6
S = 5.0/2.0
U = 0
theta_i = 0

Self = np.zeros([N_y * N_x * 4, N_y * N_x * 4], dtype= complex)


i = np.arange(N_atoms)
I = np.meshgrid(i)
ii = np.reshape(I, (N_atoms, ))
    
i_matrix = int(N_y/2)
    
t_i = (i_matrix) * N_x + (ii + borde)
theta_i = np.pi/3.0
phi_i = 0

Self [t_i * 4 + 0, t_i * 4 + 0]=J*S*cos(theta_i)-U
Self [t_i * 4 + 1, t_i * 4 + 1]=-J*S*cos(theta_i)-U
Self [t_i * 4 + 0, t_i * 4 + 1]=J*S*sin(theta_i)*exp(-ui*phi_i)
Self [t_i * 4 + 1, t_i * 4 + 0]=J*S*sin(theta_i)*exp(ui*phi_i)
Self [t_i * 4 + 2, t_i * 4 + 2]=-J*S*cos(theta_i)+U
Self [t_i * 4 + 3, t_i * 4 + 3]=J*S*cos(theta_i)+U
Self [t_i * 4 + 2, t_i * 4 + 3]=-J*S*sin(theta_i)*exp(ui*phi_i)
Self [t_i * 4 + 3, t_i * 4 + 2]=-J*S*sin(theta_i)*exp(-ui*phi_i)

#for i_matrix in range(N_y ):
#        for j_matrix in range(N_x -1):
#            
#            g_i = (i_matrix)*N_x + j_matrix 
#            g_j = (i_matrix)*N_x + j_matrix + 1
#            
#            Self [g_i, g_j]= lamda
#            Self [g_j, g_i]= -lamda
#            
##i_matrix = N_y - 1 
##for j_matrix in range(N_x - 1):
##            
##    g_i = (i_matrix)*N_x + j_matrix
##    g_j = (i_matrix)*N_x +( j_matrix + 1 )
##            
##    Self [g_i, g_j]= lamda
##    Self [g_j, g_i]= -lamda
#
#
#for i_matrix in range(N_y - 1):
#        for j_matrix in range(N_x):
#            
#            g_i = (i_matrix)*N_x + j_matrix
#            g_j = (i_matrix + 1)*N_x + j_matrix
#            
#            Self [g_i, g_j]= 2*lamda
#            Self [g_j, g_i]= -2*lamda
            
print(Self)           