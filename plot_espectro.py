#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 16:26:57 2018

@author: cristina
"""
import matplotlib.pyplot as plt
import detect_peaks as dp
plt.rcParams.update({'font.size': 13})

def espectro(spectro, spectro_spinup, spectro_spindown, spectro_13, spectro_31, spectro_24, spectro_42, row, vv, borde):
    

    spectro_chain = spectro[row, borde, :]#spectrum in the first atom
    spectro_13_one = spectro_13[row, borde, :]
    spectro_31_one = spectro_31[row, borde, :]
    
    spectro_24_one = spectro_24[row, borde, :]
    spectro_42_one = spectro_42[row, borde, :]
    
    spectro_chain_up = spectro_spinup[row, borde, :]#spectrum in the first atom
    spectro_chain_down = spectro_spindown[row, borde, :]#spectrum in the first atom
    
    ndexes = dp.detect_peaks(spectro_chain)#find the peaks
    peaks = vv[ndexes]

    minpeak = min(abs(peaks))#find the minimum
    peaks2 = abs(peaks)
    peaks = peaks.tolist()
    peaks2 = peaks2.tolist()
    i=peaks2.index(minpeak)#the index of the peak closest to zero
    titulo = vv[ndexes[i]]

    plt.figure(1)
    plt.style.use('seaborn-bright')
    plt.plot(vv, spectro_chain, linewidth=1.0, label = 'total')
    plt.plot(peaks,spectro_chain[ndexes],'y*')
    plt.xlabel('meV')
    plt.ylabel('PDOS')
    #plt.title('We use peak # %i ' %i)
    plt.savefig('results/spectro.pdf')
    

    
    plt.figure(2)
    plt.style.use('seaborn-bright')
    plt.plot(vv, spectro_chain,linewidth=1.0, label = 'Total')
    plt.plot(vv, spectro_chain_up,linewidth=0.8, label = 'Spin up')
    plt.plot(vv, spectro_chain_down,linewidth=0.8, label = 'Spin down')
    plt.xlabel('meV')
    plt.ylabel('PDOS')
    plt.legend()
    plt.savefig('results/spectro2.pdf')
    
    
    plt.figure(3)
    plt.style.use('seaborn-bright')
    #plt.plot(vv, spectro_chain, linewidth=1.0, label = 'total')
    plt.plot(vv, spectro_13_one, linewidth=0.8, label = 'G_13')
    plt.plot(vv, spectro_31_one, linewidth=0.8, label = 'G_31')
    plt.plot(vv, spectro_chain_up, linewidth=0.8, label = 'G_11')
    plt.legend()
    plt.xlabel('meV')
    plt.ylabel('PDOS')
    plt.title('Spin up')
    plt.savefig('results/spectro_13')
    
    plt.figure(4)
    plt.style.use('seaborn-bright')
    #plt.plot(vv, spectro_chain, linewidth=1.0, label = 'total')
    plt.plot(vv, spectro_24_one, linewidth=0.8, label = 'G_24')
    plt.plot(vv, spectro_42_one, linewidth=0.8, label = 'G_42')
    plt.plot(vv, spectro_chain_down, linewidth=0.8, label = 'G_22')
    plt.legend()
    plt.xlabel('meV')
    plt.ylabel('PDOS')
    plt.title('Spin down')
    plt.savefig('results/spectro_24')
    
    
    
    return(titulo, ndexes, i)
    
    
    
    
