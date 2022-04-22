import os
from models import *
import numpy as np
from scipy.integrate import solve_ivp

def lemke_howson_wrapper(A, r):
    np.savetxt('../data/A.csv', A, delimiter=',')
    np.savetxt('../data/r.csv', r, delimiter=',')
    os.system('Rscript call_lr.r')
    x = np.loadtxt('../data/equilibrium.csv', delimiter=',')
    return x

def single_extinction(t, n, A, r, tol):
    n = n[n!=0]
    return np.any(abs(n) < tol) -1

def check_constant(sol_mat, tol):
    '''
    Check if all the solutions have reached steady state (constant)
    '''
    #Get differences between solutions
    diff_sol = sol_mat[:, 1:] - sol_mat[:, 0:-1]
    #Get last 3 timepoints
    last_3 = diff_sol[:, -1:-3:-1]
    #Note that we only impose three because there are no oscillations here. 
    const = np.all(abs(last_3) < tol)
    return const

def prune_community(fun, x0, tol, args, events=single_extinction):
    '''
    Function to prune community. Every time a species goes extinct, integration
    restarts with the pruned system
    '''
    single_extinction.terminal = True
    t_span = [0, 1e6]
    #add tolerance to tuple of arguments
    args += (tol, )
    #get initial number of species
    n_sp = len(x0)
    constant = False
    while n_sp > 1 and not constant :
        sol = solve_ivp(fun, t_span, x0, events=events, args=args, 
                        method='BDF') 
        #set species below threshold to 0
        end_point = sol.y[:, -1]
        ind_ext = np.where(end_point < tol)[0]
        end_point[ind_ext] = int(0)
        n_sp = len(end_point) - len(ind_ext)
        #initial condition of next integration is end point of previous one
        x0 = end_point
        #check if solution is constant
        constant = check_constant(sol.y, tol)
    return sol
