#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 10:03:36 2018

@author: cristina
"""

"We plot Green's functions"

def plot_Green(G, N_x, N_y, N_omega, row, borde, N_atoms, vv):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    #Real and Imaginary Green's functions
    Green_Re_0 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    Green_Re_1 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    Green_Re_2 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    Green_Re_3 = np.zeros([N_y, N_x, N_omega], dtype= 'float')

    Green_Im_0 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    Green_Im_1 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    Green_Im_2 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    Green_Im_3 = np.zeros([N_y, N_x, N_omega], dtype= 'float')
    
    for i_atom in range(N_x):
        for j_atom in range(N_y):
            I = i_atom + (j_atom)*N_x

            for i_omega in range(N_omega):
             
                trup = G[I*4 + 0, I*4 + 0, i_omega]
                trdown = G[I*4 + 1, I*4 + 1, i_omega]
                trup_hole = G[I*4 + 2, I*4 + 2, i_omega]
                trdown_hole = G[I*4 + 3, I*4 + 3, i_omega]
             
                Green_Re_0[j_atom , i_atom, i_omega] = trup.real
                Green_Re_1[j_atom , i_atom, i_omega] = trdown.real
                Green_Re_2[j_atom , i_atom, i_omega] = trup_hole.real
                Green_Re_3[j_atom , i_atom, i_omega] = trdown_hole.real
             
                Green_Im_0[j_atom , i_atom, i_omega] = trup.imag
                Green_Im_1[j_atom , i_atom, i_omega] = trdown.imag
                Green_Im_2[j_atom , i_atom, i_omega] = trup_hole.imag
                Green_Im_3[j_atom , i_atom, i_omega] = trdown_hole.imag
                
    #plot Greens functions Real and Imag for all 4 Nambus
    Green_Re_0_first = Green_Re_0[row, borde, :]#first atom
    Green_Re_0_last = Green_Re_0[row, borde + N_atoms - 1, :]#last atom
    Green_Re_0_middle = Green_Re_0[row, borde + int(N_atoms/2.0), :]#middle
    Green_Re_0_out = Green_Re_0[0, 0, :]#outside

    Green_Re_1_first = Green_Re_1[row, borde, :]
    Green_Re_1_last = Green_Re_1[row, borde + N_atoms - 1, :]
    Green_Re_1_middle = Green_Re_1[row, borde + int(N_atoms/2.0), :]
    Green_Re_1_out = Green_Re_1[0, 0, :]

    Green_Re_2_first = Green_Re_2[row, borde, :]
    Green_Re_2_last = Green_Re_2[row, borde + N_atoms - 1, :]
    Green_Re_2_middle = Green_Re_2[row, borde + int(N_atoms/2.0), :]
    Green_Re_2_out = Green_Re_2[0, 0, :]

    Green_Re_3_first = Green_Re_3[row, borde, :]
    Green_Re_3_last = Green_Re_3[row, borde + N_atoms - 1, :]
    Green_Re_3_middle = Green_Re_3[row, borde + int(N_atoms/2.0), :]
    Green_Re_3_out = Green_Re_3[0, 0, :]

    Green_Im_0_first = Green_Im_0[row, borde, :]
    Green_Im_0_last = Green_Im_0[row, borde + N_atoms - 1, :]
    Green_Im_0_middle = Green_Im_0[row, borde + int(N_atoms/2.0), :]
    Green_Im_0_out = Green_Im_0[0, 0, :]

    Green_Im_1_first = Green_Im_1[row, borde, :]
    Green_Im_1_last = Green_Im_1[row, borde + N_atoms - 1, :]
    Green_Im_1_middle = Green_Im_1[row, borde + int(N_atoms/2.0), :]
    Green_Im_1_out = Green_Im_1[0, 0, :]

    Green_Im_2_first = Green_Im_2[row, borde, :]
    Green_Im_2_last = Green_Im_2[row, borde + N_atoms - 1, :]
    Green_Im_2_middle = Green_Im_2[row, borde + int(N_atoms/2.0), :]
    Green_Im_2_out = Green_Im_2[0, 0, :]

    Green_Im_3_first = Green_Im_3[row, borde, :]
    Green_Im_3_last = Green_Im_3[row, borde + N_atoms - 1, :]
    Green_Im_3_middle = Green_Im_3[row, borde + int(N_atoms/2.0), :]
    Green_Im_3_out = Green_Im_3[0, 0, :]
    
    
    plt.figure(7)
    plt.title('Real Greens function first atom')
    plt.plot(vv, Green_Re_0_first, label='up',linewidth=0.8)
    plt.plot(vv, Green_Re_1_first, label='down',linewidth=0.8)
    plt.plot(vv, Green_Re_2_first, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Re_3_first, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Real_first.pdf')

    plt.figure(8)
    plt.title('Imag Greens function first atom')
    plt.plot(vv, Green_Im_0_first, label='up',linewidth=0.8)
    plt.plot(vv, Green_Im_1_first, label='down',linewidth=0.8)
    plt.plot(vv, Green_Im_2_first, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Im_3_first, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Imag_first.pdf')

    plt.figure(9)
    plt.title('Real Greens function last atom')
    plt.plot(vv, Green_Re_0_last, label='up',linewidth=0.8)
    plt.plot(vv, Green_Re_1_last, label='down',linewidth=0.8)
    plt.plot(vv, Green_Re_2_last, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Re_3_last, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Real_last.pdf')

    plt.figure(10)
    plt.title('Imag Greens function last atom')
    plt.plot(vv, Green_Im_0_last, label='up',linewidth=0.8)
    plt.plot(vv, Green_Im_1_last, label='down',linewidth=0.8)
    plt.plot(vv, Green_Im_2_last, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Im_3_last, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Imag_last.pdf')

    plt.figure(11)
    plt.title('Real Greens function middle atom')
    plt.plot(vv, Green_Re_0_middle, label='up',linewidth=0.8)
    plt.plot(vv, Green_Re_1_middle, label='down',linewidth=0.8)
    plt.plot(vv, Green_Re_2_middle, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Re_3_middle, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Real_middle.pdf')

    plt.figure(12)
    plt.title('Imag Greens function middle atom')
    plt.plot(vv, Green_Im_0_middle, label='up',linewidth=0.8)
    plt.plot(vv, Green_Im_1_middle, label='down',linewidth=0.8)
    plt.plot(vv, Green_Im_2_middle, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Im_3_middle, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Imag_middle.pdf')

    plt.figure(13)
    plt.title('Real Greens function on SC')
    plt.plot(vv, Green_Re_0_out, label='up',linewidth=0.8)
    plt.plot(vv, Green_Re_1_out, label='down',linewidth=0.8)
    plt.plot(vv, Green_Re_2_out, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Re_3_out, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Real_out.pdf')

    plt.figure(14)
    plt.title('Imag Greens function on SC')
    plt.plot(vv, Green_Im_0_out, label='up',linewidth=0.8)
    plt.plot(vv, Green_Im_1_out, label='down',linewidth=0.8)
    plt.plot(vv, Green_Im_2_out, label='up hole',linewidth=0.8)
    plt.plot(vv, Green_Im_3_out, label='down_hole',linewidth=0.8)
    plt.legend()
    plt.savefig('Imag_out.pdf')

    
    
    