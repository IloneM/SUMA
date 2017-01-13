from inspect import getfullargspec
#from mpi4py import MPI
import suma
import cma
#import stdrun
import erttools

runparameters = {'ks': [5],
                 'dims' : [50],
                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow],
                 'tests' : range(15)}


#runparameters = {'ks': [2, 5, 10],
#                 'dims' : [50, 100, 500],
#                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
#                            cma.fcts.rosen],
#                 'tests' : range(15)}

#path_prefix_input  = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/SUMA-14-12-data/"
#path_prefix_output = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/ert/"

#    runparameters = {'ks': [2, 5, 10],
#                    'dims' : [50, 100, 200, 500],
#                    'funcs' : [cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
#                            cma.fcts.rosen],
#                    'tests' : range(15)}

path = '/home/ilone/Documents/Studies/Recherche/sumapy/results/12.1.17/'
#precs = [1e-4, 1e-6, 1e-8]
precs = list(range(6,13,2))


def doWork(k, dim, func, test, logger):
    if 'rot' in getfullargspec(func)[0]:
        args = (1,)
    else:
        args = ()

    if func.__name__ == 'rosen':
        x0 = [0] * dim
    else:
        x0 = [1] * dim

    print("running test %d for k=%d d=%d f=%s" % (test, k, dim, func.__name__))
    ests = suma.Suma(x0, 0.5, {'SUMA_k': k, 'verb_disp': 0})
    ests.optimize(func, args=args, logger=logger)

#for it in range(5):
    #opts = {'stdbehaviour': lambda: not comm.Get_rank() == 1, 'specialbehaviour': specialbehaviour}

opts = {'manager': {'path': path, 'specific_runparameters': {'funcs': [f.__name__ for f in runparameters['funcs']]}, 'splitkeys': {'funcs'}}}

erttools.ERTMPIWrapper(runparameters, doWork, precs, opts).launch()

#stdrun.launch(runparameters, doWork, opts)

