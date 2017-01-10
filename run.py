from inspect import getfullargspec
from mpi4py import MPI
import suma
import cma
import stdrun
import erttools

#runparameters = {'ks': [2, 5, 10],
#                 'dims' : [50, 100, 500],
#                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
#                            cma.fcts.rosen],
#                 'tests' : range(15)}

comm = MPI.COMM_WORLD
suffix = "fit.dat"
path_prefix_input  = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/SUMA-14-12-data/"
path_prefix_output = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/ert/"


def doWork(k, dim, func, test, rawdata):
    if 'rot' in getfullargspec(func)[0]:
        args = (1,)
    else:
        args = ()

    if fun.__name__ == 'rosen':
        x0 = [0] * dim
    else:
        x0 = [1] * dim

    print("running test %d for k=%d d=%d f=%s" % (test, k, dim, func.__name__))
    ests = suma.Suma(x0, 0.5, {'SUMA_k': k, 'verb_disp': 0})
#                     'verb_log': 100,
#                     'verb_filenameprefix': "suma_oracle_k=%d_d=%d_f=%s_%d_" %
#                                            (k, dim, func.__name__, test)})
<<<<<<< HEAD
#   if 'rot' in getfullargspec(func)[0]:
#       ests.optimize(func, args=(1,))
#   else:
#       ests.optimize(func)
=======
    ests.optimize(func, args, logger=erttools.ERTLogger(precs))

    comm.send((ests.logger.data(), rawdata), dest=1, tag=2)

def specialbehaviour(rank):
    #TODO one file per func
    #resultfiles = {fct.__name__: "result_%s.dat" % (fct.__name__) for fct in runparameters['funcs']}
    resultfiles = {fct.__name__: open(path_prefix_output + "result_%s.dat" % (fct.__name__), 'w') for fct in runparameters['funcs']}

    totalData = 1
    for item in [len(list(iterable)) for iterable in itervalues(runparameters)]:
        totalData *= item

    for resultfile in itervalues(resultfiles):
#        with open(path_prefix_output + resultfile) as resfilestream:
 #           resfilestream.write("# k d prec ert\n")
        resultfile.write("# k d prec ert\n")

    for fakeit in range(totalData):
        data = comm.recv(tag=2, source=MPI.ANY_SOURCE)
        datastr = ''
        for i in range(len(precs)):
            datastr += str(data['k']) + ' ' + str(data['d']) + ' ' + str(precs[i]) + ' ' + str(data['erts'][i]) + '\n'
        resultfiles[data['func']].write(datastr)

    for resultfile in itervalues(resultfiles):
        resultfile.close()


for it in range(5):
    opts = {'stdbehaviour': lambda: not comm.Get_rank() == 1, 'specialbehaviour': specialbehaviour}
    runparameters = {'ks': [2, 5, 10],
                    'dims' : [50, 100, 200, 500],
                    'funcs' : [cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
                            cma.fcts.rosen],
                    'tests' : range(15)}

    precs = list(range(6,13,2))

    stdrun.launch(runparameters, doWork, opts)
>>>>>>> dcc82dd9143392f6c0e5333e10221034176e2724

