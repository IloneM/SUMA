from inspect import getargspec
import suma
import cma
import stdrun

runparameters = {'ks': [1, 2, 5, 10, 20, 50],
                 'dims' : [50, 100, 200, 500, 1000, 5000, 10000],
                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
                            cma.fcts.rosen],
                }
#                 'tests' : range(10)}

#def doWork(k, dim, func, test):
def doWork(k, dim, func):
    print("running test for k=%d d=%d f=%s" % (k, dim, func.__name__))
    ests = suma.Suma([0.5] * dim, 0.5,
                     {'SUMA_k': k,
                      'verb_disp': 0,
                      'verb_filenameprefix': "suma_oracle_k=%d_d=%d_f=%s_" %
                                             (k, dim, func.__name__)})
    if 'rot' in getargspec(func)[0]:
        ests.optimize(func, args=(1,))
    else:
        ests.optimize(func)

stdrun.launch(runparameters, doWork)
