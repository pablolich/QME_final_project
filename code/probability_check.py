#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import scipy.special as sc

## CONSTANTS ##


## FUNCTIONS ##

def nc_n_N_distribution(n_c, n, N):
    '''
    Sample two subsets with length n with replacement from a set of 
    cardinality N. This function computes the probability that this two 
    subsets have n_c elements in common.
    '''
    return sc.binom(n, n_c) * sc.factorial(n)/sc.factorial(n-n_c) * \
           (1/N)**n_c * (1-1/N)**(n-n_c)

def main(argv):
    '''Main function'''
    N = 1000
    n = 10
    n_c = 3
    p = nc_n_N_distribution(n_c, n, N)
    #check with simulations
    S = np.random.uniform(0, 1, N)
    commons = 0
    runs = 0
    p_approx = np.inf
    diff = abs(p_approx - p)
    tol = 1e-9
    while diff > tol:
        #sample first set
        sample1 = set(np.random.choice(S, n))
        #sample second set
        sample2 = set(np.random.choice(S, n))
        #check number of commonn elements
        n_common = len(sample1.intersection(sample2))
        #if its n_c, add to count
        runs += 1
        if n_common == n_c:
            commons += 1
            p_approx = commons/runs
        diff = abs(p-p_approx)
        print('Theory: ', p, '\nSimulation: ', p_approx, 
              '\nError: ', diff, end = '\r')
    print(p_approx)
    print(p)
    import ipdb; ipdb.set_trace(context = 20)
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

