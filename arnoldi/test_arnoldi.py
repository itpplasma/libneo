#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 12:51:49 2016

@author: Christopher Albert
"""

import numpy as np
import numpy.linalg as npl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

Q = np.array([[1,1,0,0,0,0,0,0,0,0],[0,1,1j,0,0,0,0,0,0,0],\
              [0,0,1,1j,0,0,0,0,0,0],[0,0,0,1,1,0,0,0,0,0],\
              [0,0,0,0,1,1,0,0,0,0],[0,0,0,0,0,1,1,0,0,0],\
              [0,0,0,0,0,0,1,1j,0,0],[0,0,0,0,0,0,0,-1j,1,0],\
              [0,0,0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0,0,0]])
    
# case 1: spectral radius > 1, divergence          
L = np.diag([10+5j,.9-.01j,-.4+.2j,.3-.01j,3+2j,.2,.15j,.1,.123,.05j])

# case 2: small spectral radius, fast convergence
# L = 0.01*L

# case 3: spectral radius close to 1, slow convergence
# L = 0.08*L

M = np.dot(np.dot(Q,L),npl.inv(Q))
y = [1+0.5j,0,0,0,0,0,0,0,2j,1]

A = np.eye(10,dtype=complex) - M
b = np.dot(M,y)

xsol = np.linalg.solve(A,b)

text = ['({},{})'.format(np.real(n), np.imag(n)) for n in M.T.flatten()]
text = ' '.join(text)
with open('../MC/RUN/mmat.dat','w') as f:
    f.write(text)

text = ['({},{})'.format(np.real(n), np.imag(n)) for n in A.T.flatten()]
text = ' '.join(text)
with open('../MC/RUN/amat.dat','w') as f:
    f.write(text)
    
text = ['({},{})'.format(np.real(n), np.imag(n)) for n in y]
text = ' '.join(text)
with open('../MC/RUN/yvec.dat','w') as f:
    f.write(text)
    
text = ['({},{})'.format(np.real(n), np.imag(n)) for n in b]
text = ' '.join(text)
with open('../MC/RUN/bvec.dat','w') as f:
    f.write(text)
    
    
text = ['({},{})'.format(np.real(n), np.imag(n)) for n in xsol]
text = ' '.join(text)
with open('../MC/RUN/xsol.dat','w') as f:
    f.write(text)