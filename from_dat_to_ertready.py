import cma
#import numpy as np
#from runiterator import runGen
#from six import itervalues
import stdrun
#from collections import OrderedDict

nbTest = 15

#runparameters = OrderedDict([('dims', [50, 100, 500]),
#                ('funcs', [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
#                           cma.fcts.rosen]),
#                ('ks', [2, 5, 10]),
#                ('tests', range(nbTest)),
#               ])

runparameters = {'ks': [2, 5, 10],
                 'dims' : [50, 100, 500],
                 'funcs' : [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow,
                            cma.fcts.rosen],
                }

path_prefix = "/home/ilone/Documents/Studies/Recherche/sumapy/data-14-12/SUMA-14-12-data/"
path_prefix2 = "/home/ilone/Documents/Studies/Recherche/sumapy/data-14-12/"

#runparameterbounds = tuple(len(list(iterable)) for iterable in itervalues(runparameters))

suffix = "fit.dat"

#valuesForSpecificParams = {'fbever':[], 'fblast': [], ''}
## ##TODO link params with names of the the values represented

###lastvalues = np.array()

###maxTest = runparameterbounds[-1]

#f = None

#for (dim, func, k, test) in runGen(runparameterbounds):#mind that "test" MUST be the last parameter for stats purposes
#    dim = runparameters['dims'][dim]
#    func = runparameters['funcs'][func]
#    k = runparameters['ks'][k]
#    test = runparameters['tests'][test]
#    if test == 0:
#        f = open(path_prefix2 + "suma_oracle_k=%d_d=%d_f=%s_%s" % (k, dim, func.__name__, suffix), 'w')
#    for line in open(path_prefix + "suma_oracle_k=%d_d=%d_f=%s_%d_%s" % (k, dim, func.__name__, test, suffix), 'r').readlines():
#        selectedLine = ''
#        if len(line) > 0 and line[0] not in ('%', '#'):
#            selectedLine = line
#    f.write(selectedLine)
#    if test >= nbTest -1:
#        f.close()

#lres.append(list(map(float, line.split())))
###    lastvalues.append(
###     np.array(cma._fileToMatrix("suma_oracle_k=%d_d=%d_f=%s_%s" % (k, dim, func.__name__, suffix)))[-1,:]
###        )
###   if(test >= maxTest-1) {
###       #compute stats
###       
###       lastvalues = np.array()
###   }

def doWork(k, dim, func):
    fw = open(path_prefix2 + "suma_oracle_k=%d_d=%d_f=%s.dat" % (k, dim, func.__name__), 'w')
    for test in range(nbTest):
        try:
            with open(path_prefix + "suma_oracle_k=%d_d=%d_f=%s_%d_%s" % (k, dim, func.__name__, test, suffix), 'r') as fr:
                for line in fr.readlines():
                    selectedLine = ''
                    if len(line) > 0 and line[0] not in ('%', '#'):
                        selectedLine = line
                fw.write(selectedLine)
                fr.close()
        except IOError:
            with open(path_prefix2 + "expe_not_done", 'a') as fwe:
                fwe.write("k=%d dim=%d func=%s test=%d\n" % (k, dim, func.__name__, test))
                fwe.close()
    fw.close()

stdrun.launch(runparameters, doWork)
