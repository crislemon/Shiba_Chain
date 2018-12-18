#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 16:26:57 2018

@author: cristina
"""
import matplotlib.pyplot as plt
import detect_peaks as dp

def espectro(spectro, row, vv, borde):
    

    spectro_chain = spectro[row, borde, :]#spectrum in the first atom
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
    plt.plot(vv, spectro_chain,linewidth=0.8)
    plt.plot(peaks,spectro_chain[ndexes],'y*')
    plt.xlabel('meV')
    plt.ylabel('PDOS')
    plt.title('We use peak # %i ' %i)
    plt.savefig('results/spectro.pdf')
    
    
    return(titulo, ndexes, i)
