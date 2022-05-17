#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from essential_tools import *
import numpy as np
import pandas as pd
from scipy.stats import beta, dirichlet
from itertools import combinations
import math
from scipy import spatial
import matplotlib.pylab as plt

## CONSTANTS ##


## FUNCTIONS ##

def abund_approx(ab_opt_comp, C, l, D, fun):
    '''
    Approximate abundances of subcommunity as a function of the optimal
    abundances plus a leakage term
    '''
    I = np.identity(D.shape[0])
    B = I - l*D
    return np.linalg.inv(C @ C.T) @ C @ (l*D@C.T@np.linalg.inv(C@B@C.T)@C +\
                                         B) @ fun

def optimal_abundances(C_k, C, x):
    x = x.reshape(len(x), 1)
    return (np.linalg.inv(C_k @ C_k.T) @ C_k @ C.T @ x).T[0]

def comm_function(C, abundances):
    x = abundances.reshape(len(abundances), 1)
    return (C.T@x).T[0]

def randbin(n, m, count):
    '''
    Create a nxm binary random matrix conditioned to the constrained that there
    need be "count" ones in each row.
    '''
    #output to assign ones into
    result = np.zeros((n, m))
    #simulate sampling with replacement in one axis
    col_ind = np.argsort(np.random.random(size=(n, m)), axis=1)
    #turn it into a mask over col_ind using a clever broadcast
    try:
        mask = np.arange(m) < np.round(count).reshape(n, 1)
    except:
        mask = np.arange(m) < np.round(count)
    #apply the mask not only to col_ind, but also the corresponding row_ind
    col_ind = col_ind[mask]
    row_ind = np.broadcast_to(np.arange(n).reshape(-1, 1), (n, m))[mask]
    #Set the corresponding elements to 1
    result[row_ind, col_ind] = 1
    return result

def all_comms(n):
    '''
    Get all 2**n -1 (leave the bare ground out) possible combination of n 
    species in binary
    '''
    n_rows = 2**n - 1
    vec = np.zeros(shape=(n_rows, n), dtype = int)
    binary_string = np.zeros(shape=(n_rows, 1), dtype = '<U'+str(n))
    i = 0
    while i < n_rows:
        #transform to binary (start at 1 to leave out the bare ground)
        binary = dec2bin(int(i+1))
        #complete with zeros
        vec[i, :] = np.hstack((np.zeros(n - len(binary)), binary))
        i += 1
    return vec

def dec2bin(x):
    '''
    get binary number representation of a number
    '''
    return list(map(int, bin(x)[2:]))

def build_positive_definite_D(m):
    #build crossfeeding matrix (row stochastic)
    D = dirichlet.rvs(alpha=np.ones(m), size=m)
    #get eigenvalues
    lam = np.linalg.eigvals(D)
    #shift matrix such that the minimum eigenvalue becomes positive
    min_val = min(lam)
    shift = min_val.real*np.identity(m)
    shift_D = D - shift 
    #make row stochastic
    D_final = 1/np.sum(shift_D, axis=1)*np.identity(m) @ shift_D
    return D_final

def main(argv):
    '''Main function'''
    #set parameters
    n, m = (5, 5)
    build_positive_definite_D(m)
    l_vec = np.arange(0.01, 0.2, 0.01)
    #Set growth and death rates for each simulation
    d = min(1-l_vec)*np.repeat(np.random.uniform(0, 1), n)
    r = 1+d+np.random.uniform(0, 1, m)
    i = 0
    error_approx_vec = np.zeros(len(l_vec))
    error_opt_vec = np.zeros(len(l_vec))
    for l in l_vec:
        #select a particular ROS
        a = 1
        b = 1
        #sample number of preferences of each species from beta 
        #distribution
        n_pref = beta.rvs(a, b, loc=1, scale=m-1, size=n)
        #build preference matrix
        C = randbin(n, m, n_pref)    
        #normalize it
        C = (C.T * 1/np.round(n_pref)).T
        #build crossfeeding matrix (row stochastic)
        D = build_positive_definite_D
        #get parameters of equivalent lotka volterra
        #A = C@C.T
        #rho = C@r - d
        I = np.identity(m)
        A = (1-l)*(C@(I-l*D)@C.T)
        rho = (1-l)*C@r-d
        #initialize community object
        glv_community = Community(np.ones(n), GLV, A=A, r=rho)
        #assembly community and forget its history
        glv_community.assembly()
        #remove preferences from extinct species
        C_ext = C[glv_community.presence, :]
        #forget first assembly history
        glv_community = glv_community.delete_history()
        #Measure Cxstar
        f = comm_function(C_ext, glv_community.n)
        #get all subcommunities of current one
        vec = all_comms(glv_community.richness)
        #get number of species in each subcommunity
        n_spp = np.sum(vec, axis=1)
        #select subcommunities with j species
        j = glv_community.richness//2
        #select all subcommunities with j species
        ind = np.where(n_spp == j)[0]
        finished = False
        while not finished:
            #sample one of these subcommunities
            ind_sample = np.random.choice(ind)
            sub_comm = vec[ind_sample, :]
            #get indices of species to be removed
            ind_rem = np.where(sub_comm == 0)[0]
            #form subcommunity
            glv_sub_comm = glv_community.remove_spp(ind_rem)
            #assembly the subcommunity
            glv_sub_comm.assembly()
            #only work with those communities with j species after 
            #assembly
            if glv_sub_comm.richness == j: 
                #get abundances
                abundances = np.array(glv_sub_comm.n)
                #remove preferences from extinct species
                C_sub = C_ext[glv_sub_comm.presence, :]
                #get optimal abundances
                abund_opt = optimal_abundances(C_sub, C_ext, glv_community.n) 
                #get abundances approximately
                approx_ab = abund_approx(abund_opt, C_sub, l, D, f)
                #get errors
                error_vec = abundances-approx_ab
                error_approx_vec[i] = np.dot(error_vec, error_vec)
                error_vec_opt = abundances-abund_opt
                error_opt_vec[i] = np.dot(error_vec_opt, error_vec_opt)
                finished = True
                i+=1

    plt.plot(l_vec, np.log(1+error_approx_vec))
    plt.plot(l_vec, np.log(1+error_opt_vec))
    plt.show()
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

