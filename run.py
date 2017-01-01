from inspect import getargspec
from mpi4py import MPI
import suma
import cma
import stdrun

#runparameters = {'ks': [2, 5, 10],
#                 'dims' : [50, 100, 500],
#                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
#                            cma.fcts.rosen],
#                 'tests' : range(15)}

runparameters = {'ks': [2, 5, 10],
                 'dims' : [50, 100, 200, 500],
                 'funcs' : [cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
                            cma.fcts.rosen],
                 'tests' : range(15)}

def doWork(k, dim, func, test):
    print("running test %d for k=%d d=%d f=%s" % (test, k, dim, func.__name__))
#   ests = suma.Suma([0.5] * dim, 0.5,
#                    {'SUMA_k': k,
#                     'verb_disp': 0,
#                     'verb_log': 100,
#                     'verb_filenameprefix': "suma_oracle_k=%d_d=%d_f=%s_%d_" %
#                                            (k, dim, func.__name__, test)})
#   if 'rot' in getargspec(func)[0]:
#       ests.optimize(func, args=(1,))
#   else:
#       ests.optimize(func)

stdrun.launch(runparameters, doWork)
