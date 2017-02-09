import sumatricheur as st
import matplotlib.pyplot as plt
import numpy as np
import threading as th
from math import ceil,floor
from importlib import reload
import stdrun

def plot(path):
    dat = np.array(st.cma._fileToMatrix(path))
    x = list(range(dat.shape[0]))
    if dat.shape[1] > 1:
        nbdat = dat.shape[1]
        c = int(ceil(nbdat ** 0.5))
        l = nbdat // c + (1 if nbdat % c > 0 else 0)
        f, ax = plt.subplots(l,c)
        for i in range(nbdat):
            ax[i // c, i % c].scatter(x, dat[:,i], marker='.')
    else:
        plt.figure()
        plt.scatter(x, dat, marker='.')

def genpath(d,D):
    return 'results/NDCG/D=%d_d=%d' % (D,d)

def run(tdD):
    d,D = tdD
    es = st.Suma([1] * D, 0.5, {'SUMA_d': d, 'verb_log': 0, 'SUMA_nrp': 25, 'verb_filenameprefix': genpath(d,D)})
    es.optimize(st.cma.fcts.ellirot)

#if __name__ == 'main':
stdrun.launch({'tdDs': [(2,100),(5,100),(10,100),(20,100),(30,100),(2,200),(5,200),(10,200),(20,200),(30,200),(40,200)]}, run)

#thqueue = []
#
#for d,D in [(2,100),(5,100),(10,100),(20,100),(30,100),(2,200),(5,200),(10,200),(20,200),(30,200),(40,200)]:
#    newth = th.Thread(target=run, args=(d,D))
#    newth.start()
#    thqueue.append(newth)
#
#for th in thqueue:
#    th.join()
    

