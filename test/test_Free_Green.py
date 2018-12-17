#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 10:02:15 2018

@author: cristina
"""

import numpy as np
import cmath as cm
import scipy.spatial.distance

N_x = 3
N_y = 2


g = np.zeros([N_x, N_y, 4, 4], dtype = complex)