#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 16:48:30 2018

@author: cristina
"""

import matplotlib.pyplot as plt
import numpy as np

def profile(N_x, titulo, spectro, row, ndexes, i):
    plt.figure(3)
    profile = np.zeros(N_x, dtype= 'float')
    profile = spectro[row,:,ndexes[i]]
    plt.plot(profile,'--bo', linewidth=0.8)
    plt.title('Profile at E = %f meV' %titulo)
    plt.xlabel('Atom index')
    plt.ylabel('PDOS')
    plt.savefig('results/profile.pdf')
