import numpy as np

def translated_sphere(X, args=()):
    c = [1] * len(X);
    if len(args) > 1 and len(X) == len(args):
        c = args;
    trX = X - c;
    return np.dot(trX,trX);

