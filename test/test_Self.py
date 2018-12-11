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


N_atoms = 1
N_x = 3
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

Self2 = np.zeros([N_y * N_x * 4, N_y * N_x * 4], dtype=complex)
se2 = np.zeros([N_y * N_x, N_y * N_x, 4, 4], dtype=complex)
     
for i_matrix in range(N_y):
    for j_matrix in range(N_x - 1):
            
        g_i = (i_matrix)*N_x + j_matrix
        g_j = (i_matrix)*N_x + (j_matrix + 1)
                
        se2[g_i, g_j, 0, 1]= lamda
        se2[g_i, g_j, 1, 0]= - lamda
        se2[g_i, g_j, 2, 3]= - lamda
        se2[g_i, g_j, 3, 2]= lamda
            
            
        se2[g_j, g_i, 0, 1]= -lamda
        se2[g_j, g_i, 1, 0]=  lamda
        se2[g_j, g_i, 2, 3]=  lamda
        se2[g_j, g_i, 3, 2]= -lamda
        ###cambiar signo o no????

    
#the coupling along y direction (sigma_x)
    
for i_matrix in range(N_y - 1):    
    for j_matrix in range(N_x):
            
        g_i = (i_matrix)*N_x + j_matrix
        g_j = (i_matrix + 1)*N_x + j_matrix
            
        se2 [g_i, g_j, 0, 1]= -ui * lamda
        se2 [g_i, g_j, 1, 0]= -ui * lamda
        se2 [g_i, g_j, 2, 3]= -ui * lamda
        se2 [g_i, g_j, 3, 2]= -ui * lamda

        se2 [g_j, g_i, 0, 1]= ui * lamda
        se2 [g_j, g_i, 1, 0]= ui * lamda
        se2 [g_j, g_i, 2, 3]= ui * lamda
        se2 [g_j, g_i, 3, 2]= ui * lamda
            
    
    
    
    
for i_matrix in range(N_y * N_x):
    for j_matrix in range(N_y * N_x):
        for i in range (4):
            for t in range (4):
                Self2[i_matrix*4+i,j_matrix*4+t] = se2[i_matrix,j_matrix,i,t]


Self2 = np.matrix(Self2)            
print(Self2)   
test = np.allclose(Self2, Self2.getH(), atol=1e-8)
print(test)
print(Self)
test2 = np.allclose(Self, np.transpose(Self), atol=1e-8)
print(test2)


        