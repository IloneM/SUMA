import suma
import cma
from mpi4py import MPI
from inspect import getargspec
from os.path import isfile
from future.utils import iteritems

ks = [1, 2, 5, 10, 20, 50]
dims = [50, 100, 200, 500, 1000, 5000, 10000]
funcs = [cma.fcts.cigar, cma.fcts.tablet, cma.fcts.elli, cma.fcts.diffpow, cma.fcts.rosen]
nb_test = 10

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

runsqueue = [(len(ks), len(dims), len(funcs), nb_test)] * (comm.Get_size()-1)#because first is master and does not work

if isfile('current_run.dat'):
    lrnf = eval(open('current_run.dat', 'r').read())#last run not finished
else:
    lrnf = (0, 0, 0, 0)

if rank == 0:
    for k in range(len(ks)):
        for d in range(len(dims)):
            for f in range(len(funcs)):
                for test in range(nb_test):
                    if((k, d, f, test) < lrnf):
                        continue#because changing vars directly has no effect on for loops; Note: this is a tricky and ugly workaround!!
                    runner = comm.recv(tag=0, source=MPI.ANY_SOURCE)
                    run_data = {'ki': k, 'di': d, 'fi': f, 'ti': test}
                    comm.send(run_data, dest=runner, tag=1)
                    print("running test %d/%d for k=%d d=%d f=%s" % (test+1, nb_test, ks[k], dims[d], funcs[f].__name__))

                    runsqueue[runner-1] = (k, d, f, test)
                    if lrnf < min(runsqueue):
                        lrnf = min(runsqueue)
                        open('current_run.dat', 'w').write(str(lrnf))

    for i in range(comm.Get_size()-1):
        runner = comm.recv(tag=0, source=MPI.ANY_SOURCE)
        comm.send(False, dest=runner, tag=1)

        runsqueue[runner-1] = (len(ks), len(dims), len(funcs), nb_test)
        if lrnf < min(runsqueue):
            lrnf = min(runsqueue)
            open('current_run.dat', 'w').write(str(lrnf))

#if rank == 1:
#TODO: compute stats

else:
    while True:
        comm.send(rank, dest=0, tag=0)
        data = comm.recv(source=0, tag=1)
        if not data:
            break
        ests = suma.Suma([0.5] * dims[data['di']], 0.5, {'verb_disp': 0, 'SUMA_k': ks[data['ki']], 'verb_filenameprefix': "suma_oracle_k=%d_d=%d_f=%s_run=%d_" % (ks[data['ki']], dims[data['di']], funcs[data['fi']].__name__, data['ti'])})
        if 'rot' in getargspec(funcs[data['fi']])[0]:
            ests.optimize(funcs[data['fi']], [1])
        else:
            ests.optimize(funcs[data['fi']])

