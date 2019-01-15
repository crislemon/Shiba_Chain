#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 10:43:11 2019

@author: cristina
"""

import theano
a = theano.tensor.vector() # declare variable
b = theano.tensor.vector() # declare variable
out = a ** 2 + b ** 2 + 2 * a * b              # build symbolic expression
f = theano.function([a, b], out)   # compile function
print(f([0, 1, 2],[0, 1, 2]))
