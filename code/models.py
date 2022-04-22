import numpy as np

def GLV(t, x, A, r, tol):
    return (x*(r + A@x)).T

