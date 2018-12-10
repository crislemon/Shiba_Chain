# -*- coding: utf-8 -*-

import numpy as np
from numpy.linalg import inv

#we solve Dyson's equation
def Dyson_eq(Go , Self , N_x, N_y):
    
    
    Id = np.identity(4 * N_y * N_x)
   
    matrx_inv = inv(Id -np.dot(Go, Self))
    gg = np.dot(matrx_inv , Go)
    
    
    
    
    
    
    return(gg)
