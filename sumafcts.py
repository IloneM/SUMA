import numpy as np
from cma import Rotation

rotate = Rotation()

def translated_sphere(X, args=()):
    c = [1] * len(X);
    if len(args) > 1 and len(X) == len(args):
        c = args;
    trX = X - c;
    return np.dot(trX,trX);

def weightednorm(X, args=()):
    D = len(X)
    Z = rotate(X)
    res = 0
    basis = 1. + 7./D
    for i in range(D):
        #res += 10 ** (-i) * Z[i] ** 2
        res += basis ** (D/2-i) * Z[i] ** 2
        #res += Z[i] ** (2 * (i+1))
    return res

