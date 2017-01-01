from mpi4py import MPI
import cma
import stdrun
import numpy as np
#from operator import mul
from six import itervalues

nbTest = 15
precs = [1e-4, 1e-6, 1e-8]
suffix = "fit.dat"
path_prefix_input  = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/SUMA-14-12-data/"
path_prefix_output = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/ert/"

runparameters = {'ks': [2, 5, 10],
                 'dims' : [50, 100, 500],
                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
                            cma.fcts.rosen],
                }

#runparameters = {'ks': [2, 5],
#                 'dims' : [500],
#                 'funcs' : [cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
#                            cma.fcts.rosen],
#                }
comm = MPI.COMM_WORLD

def doWork(k, dim, func):
#    #extract all datas in a matrix..
#    rawdata = np.array(cma._fileToMatrix(path_prefix + "suma_oracle_k=%d_d=%d_f=%s.dat" % (k, dim, func.__name__)))
#    #..and then format all the relevant in an array of tuples
#    rawdata = [tuple(it) for it in rawdata[:,(1,5)]]
#    RTs, RTus, nbS, nbUNS = (0,0,0,0)
#    for fctevals, best in rawdata:
#        if best > 
    #BE CAREFUL: using tuple([[0] * len(precs)]*4) does not work because then all items in tuple refers to same array
    RTs, RTus, nbS, nbUS = tuple([0] * len(precs) for fakeit in range(4)) #BE CAREFUL: see above
    #RTs, RTus, nbS, nbUNS = tuple([{precit: 0 for precit in range(len(precs))}] * 4)
    for test in range(nbTest):
        #extract all datas in a matrix..
        rawdata = np.array(cma._fileToMatrix(path_prefix_input + "suma_oracle_k=%d_d=%d_f=%s_%d_%s" % (k, dim, func.__name__, test, suffix)))
        #..and then format all the relevant in an array of tuples
        rawdata = [tuple(it) for it in rawdata[:,(1,5)]]
        precit = 0
        lastfcteval = 0
        for fcteval, best in rawdata:
            if best < precs[precit]:
                nbS[precit] += 1
                RTs[precit] += lastfcteval
                precit += 1
                if precit >= len(precs):
                    break
            lastfcteval = fcteval
        for remainingprecit in range(precit, len(precs)):
            nbUS[remainingprecit] += 1
            RTus[remainingprecit] += lastfcteval
    ERT = [0] * len(precs)
    for precit in range(len(precs)):
        assert nbS[precit] + nbUS[precit] == nbTest
        if nbS[precit] > 0:
            RTs[precit]  /= nbS[precit]
        if nbUS[precit] > 0:
            RTus[precit] /= nbUS[precit]
        ps = nbS[precit] / nbTest
        if ps > 0:
            ERT[precit] = RTs[precit] + (1 - ps) / ps * RTus[precit]
        else:
            ERT[precit] = float("inf")
    comm.send({'func': func.__name__, 'k': k, 'd': dim, 'erts': ERT}, dest=1, tag=2)

def specialbehaviour():
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
    #files_michele = [[open(path_prefix_output + 'michele/results_f=%s_prec=%.1f' %
    #                  (runparameters['funcs'][j].__name__, precs[i]), 'w')
    #                  for i in range(len(precs))]
    #                  for j in range(len(runparameters['funcs']))]
    #files_marc = 
    
#                selectedLine = ''
#                if len(line) > 0 and line[0] not in ('%', '#'):
#                    selectedLine = line
#            fw.write(selectedLine)
#            fr.close()
#    except IOError:
#        with open(path_prefix2 + "expe_not_done", 'a') as fwe:
#            fwe.write("k=%d dim=%d func=%s test=%d\n" % (k, dim, func.__name__, test))
#            fwe.close()
#fw.close()

opts = {'stdbehaviour': lambda: not comm.Get_rank() == 1, 'specialbehaviour': specialbehaviour}

stdrun.launch(runparameters, doWork, opts)
