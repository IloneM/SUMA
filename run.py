#from inspect import getfullargspec
#from mpi4py import MPI
import sumatricheur as st
#import cma
import stdrun
#import erttools

runparameters = {'Ds' : [50, 100, 200],
		 'ds' : [5, 10, 20, 30, 40],
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

path = '/home/ilone/SUMA-chal/results/'
#precs = [1e-4, 1e-6, 1e-8]
#precs = list(range(6,13,2))

def genpath(d,D,test):
    return path + 'D=%d_d=%d_test=%d' % (D,d,test)

def doWork(D, d, test):
    print("running test %d for d=%d D=%d" % (test, d,D))
    es = st.Suma([1] * D, 0.5, {'SUMA_d': d, 'verb_log': 0, 'verb_disp': 0, 'SUMA_nrp': 25, 'verb_filenameprefix': genpath(d,D,test)})
    es.optimize(st.cma.fcts.ellirot)

#for it in range(5):
    #opts = {'stdbehaviour': lambda: not comm.Get_rank() == 1, 'specialbehaviour': specialbehaviour}

#opts = {'manager': {'path': path, 'specific_runparameters': {'funcs': [f.__name__ for f in runparameters['funcs']]}, 'splitkeys': {'funcs'}}}

#erttools.ERTMPIWrapper(runparameters, doWork, precs, opts).launch()

stdrun.launch(runparameters, doWork)#, opts)

