from mpi4py import MPI
from os.path import isfile
from six import iteritems,itervalues,iterkeys
from runiterator import runGen
from collections import OrderedDict
from os import rename
import hashlib

def launch(runparameters, workfunc=lambda *args, **kwargs: None):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    #if is std dict then sort it by key so that rp becomes identificable
    if not isinstance(runparameters, OrderedDict):
        runparameters = OrderedDict(sorted(iteritems(runparameters), key=lambda t: t[0]))
    runprogressfile = "currentrun_" + hashlib.md5(str(runparameters).encode('utf-8')).hexdigest()
    if isfile(runprogressfile):
        lrnf = eval(open(runprogressfile, 'r').read())#last run not finished
    else:
        lrnf = (0,) * len(runparameters)

    #if it's a master
    if rank == 0:
        #runparameterslst = list(list(itervalues(v)) if type(v) == dict else list(v) for v in itervalues(runparameters))
        #runparameterbounds = tuple(len(iterable) for iterable in runparameterslst)
        runparameterkeys = tuple(k for k in iterkeys(runparameters))
        runparameterbounds = tuple(len(list(iterable)) for iterable in itervalues(runparameters))

        runsqueue = [runparameterbounds] * (comm.Get_size()-1)#because first is master and does not work

        for it in runGen(runparameterbounds):
            if it < lrnf:
                continue#because changing vars directly has no effect on for loops; Note: this is a tricky and ugly workaround!!
            runner = comm.recv(tag=0, source=MPI.ANY_SOURCE)
            run_data = {runparameterkeys[i]: it[i] for i in range(len(runparameters))}
            comm.send(run_data, dest=runner, tag=1)

            runsqueue[runner-1] = it
            if lrnf < min(runsqueue):
                lrnf = min(runsqueue)
                open(runprogressfile, 'w').write(str(lrnf))

        #manage the process still running after the end of the running loop
        for i in range(comm.Get_size()-1):
            runner = comm.recv(tag=0, source=MPI.ANY_SOURCE)
            comm.send(False, dest=runner, tag=1)

            runsqueue[runner-1] = runparameterbounds
            if lrnf < min(runsqueue):
                lrnf = min(runsqueue)
                open(runprogressfile, 'w').write(str(lrnf))

        rpfileback = runprogressfile
        if isfile(runprogressfile + '.bak'):
            nbSaves = 2
            while isfile(runprogressfile + '_' + str(nbSaves) + '.bak'):
                nbSaves += 1
            rpfileback += '_' + str(nbSaves)
            print("You have %d old backups. You may considerer deleting it!" % (nbSaves-1))
        rpfileback += '.bak'
        rename(runprogressfile, rpfileback)

    #if rank == 1:
    #TODO: compute stats

    #otherwise is a worker
    else:
        while True:
            comm.send(rank, dest=0, tag=0)
            data = comm.recv(source=0, tag=1)
            if not data:
                break
            data = {k[:-1]: runparameters[k][v] for k, v in iteritems(data)}
            workfunc(**data)
