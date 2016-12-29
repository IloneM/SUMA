import cma
import stdrun

nbTest = 15
suffix = "fit.dat"
path_prefix = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/SUMA-14-12-data/"
path_prefix2 = "/home/ilone/Documents/Studies/Recherche/sumapy/data-std/"

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
