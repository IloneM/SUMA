import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from os import walk
from os import listdir as ls
from os.path import isfile,join
import re
from cma import _fileToMatrix
from six import itervalues
from math import log
import sys
import getopt

if __name__ == "__main__":
    opts, args = getopt.getopt(sys.argv[1:], "f:t:", ["func=", 'type='])
    for o, a in opts:
        print(o)
        if o in ("-f", "--func"):
            func = a
        elif o in ("-t", '--type'):
            plttype = a

inputpath = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/ert/"

#launch regex representing result file on all files matched in result dir..
funcs = [re.search('result_(.+).dat', f) for f in ls(inputpath) if isfile(join(inputpath, f))]
#..and extract func name: data (as dict) on matched files thanks to group match in regex
funcs = {f.group(1): _fileToMatrix(inputpath + f.group()) for f in funcs if f}

#for (dirpath, dirnames, filenames) in walk(inputpath):
#    f.extend(filenames)

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

plt.figure("ERTs for function %s (%s flavour)" % (func, plttype))

if plttype == "marc":
    splt=1

    data = np.array(sorted([row for row in funcs[func]]))
    X = sorted([-log(prec,10) for prec in set(data[:,2])])
    for d in sorted(list(set(data[:,1]))):
        plt.subplot(2,2,splt)
        splt += 1

        for k in sorted(list(set(data[:,0]))):
            #Z = np.ndarray(shape = (len(Y), len(X)), dtype=float)
            plt.plot(X, [row[3] for row in sorted(
                [row for row in data if row[1] == d and row[0] == k], key=lambda k: k[2], reverse=True)], label="k="+str(k))

        plt.title("ERT for dim " + str(d))
        plt.xlabel("-log(prec)")
        plt.ylabel("ERT")
        plt.legend(loc='best')
else:
    splt=1
    for prec in sorted(list(set(np.array(funcs[func])[:,2]))):
        plt.subplot(2,2,splt)
        splt += 1
        data = np.array(sorted([row for row in funcs[func] if row[2] == prec]))

        X = sorted(list(set(data[:,0])))
        Y = sorted(list(set(data[:,1])))
        Z = np.ndarray(shape = (len(Y), len(X)), dtype=float)

        dataZ = {(row[0],row[1]): row[3] for row in data}
        for i in range(len(X)):
            for j in range(len(Y)):
                Z[j][i] = dataZ[(X[i],Y[j])]

        V = [int(min(itervalues(dataZ)))]
        preci = False
        for i in sorted(itervalues(dataZ)):
            if preci:
                V.append(int((i+preci)/2))
            V.append(int(i)+1)
            preci = i

        CS = plt.contour(X, Y, Z, V, colors=[tuple(np.random.rand(3)) for fakeit in range(len(V))])
        plt.clabel(CS, inline=1, fontsize=10)
        plt.colorbar(CS, shrink=0.8, extend='both')
        plt.title("ERT for precision " + str(prec))
        plt.xlabel("k")
        plt.ylabel("d")
#    plt.yscale('log')

plt.show(block=True)
