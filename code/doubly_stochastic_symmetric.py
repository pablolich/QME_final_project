#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt

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

def sample_laplacian_matrix(S):
    '''
    Given the matrix S, sample a Laplacian matrix that fullfils the constraints
    1. Lii < 1 - Sii
    2. sum_j L_ij = -Lii
    3. Lij > Sij
    '''
    #size of S
    m = len(S)
    #sample elements of the diagonal 
    Lii = np.array([np.random.uniform(0, 1-S[i, i]) for i in range(m)])
    #sample off diagonal elements fullfilling constraint 3
    Lij = np.zeros(m, m)
    for i in range(m):
        for j in range(m):
            if i == j:
                L[i,j] = Lii[i]
            else:
                Li_vec = np.array([np.random.uniform(-S[i, j]-1, -S[i, j]) \
                                   for j in range(m-1)])
                #make sure it sums to the diagonal
                Li_vec = Lii[i]*Li_vec/sum(Li_vec)
                #here I have to run the algorithm that "adds the missing part 
                #to those Lij < Sij by substracting from those Lij > Sij in a 
                #ranked fashion until all Lij > Sij."
                #insert the diagonal elements
    return L

def stochastic_vector(length):
    '''
    sample a stochastic vector
    '''
    return np.random.dirichlet(np.ones(length))

def bound_components(vector, lower_bounds):
    '''
    Given a vector which components add up to a constant, modify them such that
    each is bounded from below, while still satisfying the sum constraint.

    Parameters:

        vector (1xn array): vector to be fine tunned
        bounds (1xn array): vector of bounds for each component
    '''
    #identify components that do and don't satisfy the constraint
    diff = vector - lower_bounds
    #check if problem is feasible
    surplus = sum([i for i in diff if diff > 0])
    shortage = sum([i for i in diff if diff < 0])
    if surplus < shortage: 
        print("unfeasible problem")
        raise ValueError
    #rank differences, and order both bounds and vector according to that
    rank_ind = np.argsort(diff)
    vector_rnk = vector[rank_ind]
    bounds_rnk = lower_bounds[rank_ind]


    #plot for sanity check
    plt.plot(np.arange(len(vector)), vector_rnk - bounds_rnk)
    plt.show()
    import ipdb; ipdb.set_trace(context = 20)


    #CHANGE THE SIGN OF ALL THE COMPONENTS TO BE NEGATIVE
    return 0

def build_L(size, S):
    '''
    build symmetric L satisfying all constraints
    1. Lii < 1 - Sii
    2. sum_j L_ij = -Lii
    3. Lij > Sij
    '''
    #sample diagonal satisfying (1)
    diag = np.array([np.random.uniform(0, 1-S[i, i]) for i in range(m)])
    L = np.zeros(shape=(size, size))
    #sample rows
    for i in range(size)+1:
        #add diagonal element
        L[i-1, i-1] = diag[i]
        #get stochastic vector
        vec = stochastic_vector(size - i)
        vec_sum = sum(L[i-1, 0:i])
        #make rows sum to vec_sum
        vec = vec_sum*vec
        #bound vector
        bound_vec = bound_components(vec, S[i-1, i:])
        #add to matrix L on both the rows and columns
        L[i-1, i:] = bound_vec
        L[i:, i-1] = bound_vec
    return L

def main(argv):
    '''Main function'''
    #leakage
    l = 0.99
    N = 20
    #create random matrix
    M = np.random.random((N, N))
    #Make symmetric and positive definite
    A = M@M.T
    P = generate_doubly_stochastic(N, A)
    S = P@A@P
    #create sign fix matrix
    Q = -1/(N-1)*np.ones(shape=(N, N))
    np.fill_diagonal(Q, np.ones(N))
    B = (1-l)*(S + Q) 
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
     

