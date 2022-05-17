#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np

## CONSTANTS ##


## FUNCTIONS ##

def I(n):
    '''
    Build identity matrix
    '''
    return np.identity(n)

def diag_definite_matrix(N, M0):
    '''
    Generate a diagonally dominant positive definite matrix from matrix M0
    '''
    #get eigenvalues
    eig = np.linalg.eigvals(M0)
    if np.any(eig<0):
        #make diagonally dominant
        M = M0 - I(N)*min(eig.real)
    else:
        M = M0
    return M

def generate_doubly_stochastic(N, M0, epsilon=1e-9):
    '''
    Generate a doubly stochastic matrix of dimension N, starting form matrix 
    M0
    '''
    converged = False
    #repeat the following until convergence
    D1 = np.eye(N)
    D2 = np.eye(N)
    while not converged:
        #make the matrix row stochastic
        Dv = 1/np.sum(M0, axis=1)*I(N)
        M_row = Dv @ M0
        Dw = 1/np.sum(M_row, axis=0)*I(N)
        #make the matrix column stochastic
        M_col = M_row @ Dw
        #Record matrices
        D1 *= Dv
        D2 *= Dw
        #get sum of columns
        col_sum = np.sum(M_col, axis = 0)
        row_sum = np.sum(M_col, axis = 1)
        #check for convergence
        converge_col = np.all((col_sum > 1-epsilon) & (col_sum < 1 + epsilon))
        converge_row = np.all((row_sum > 1-epsilon) & (row_sum < 1 + epsilon))
        converged = np.all((converge_col) & (converge_row))
        M0 = M_col
    #Find the constant relating them both
    w = np.diag(D2)/np.diag(D1)
    #Find unique D
    D = w**(1/2)*D1
    return D

def build_perm_matrix(permutation_vector):
    '''
    Get corresponding permutation matrix from permutation vector
    '''
    #get size of matrix
    n = len(permutation_vector)
    #preallocate permutation matrix
    P = np.zeros(shape=(n, n))
    #loop over rows
    for i in range(n): 
        P[i,permutation_vector[i]]=1
    return P

def pos_def_permutation(n):
    '''
    sample a permutation with a positive definite permutation matrix
    '''
    negative = True
    while negative:
        #sample permutation
        p_vec = np.random.permutation(n)
        #transform into matrix
        P = build_perm_matrix(p_vec)
        #check eigenvalues
        eig = np.linalg.eigvals(P)
        #repeat if any is negative
        if np.all(eig>0):
            negative = False
    return P

def main(argv):
    '''Main function'''
    #leakage
    l = 0.5
    N = 5
    #create random matrix
    M = np.random.random((N, N))
    #Make symmetric and positive definite
    A = M@M.T
    S = generate_doubly_stochastic(N, A)
    #create sign fix matrix
    Sfix = -1/(N-1)*np.ones(shape=(N, N))
    np.fill_diagonal(Sfix, np.ones(N))
    B = (1-l)*(S@A@S + Sfix) 
    D = 1/l*(np.eye(N)-B)
    import ipdb; ipdb.set_trace(context = 20)
    #sample positive definite permutation matrix
    A = pos_def_permutation(N)
    D = 1/2*(A + A.T)
    import ipdb; ipdb.set_trace(context = 20)
    #symmetric = False
    #while not symmetric:
    #    #Make positive definite
    #    D_pos = diag_definite_matrix(N, A)
    #    #Make doubly stochastic
    #    DDS = generate_doubly_stochastic(N, D_pos)
    #    #check if symmetric
    #    symm = DDS - DDS.T
    #    A = DDS
    #    print(np.max(symm))
    #    if np.all(symm < 1e-3):
    #        symmetric = True
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

