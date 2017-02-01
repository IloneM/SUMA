import numpy as np
from owncma import Rotation

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
    for i in range(D):
        res += 10 ** (-i) * Z[i] ** 2
    return res

