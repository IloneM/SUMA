from mpi4py import MPI
from os import rename,remove
from os.path import isfile
from six import iteritems,itervalues,iterkeys
from runiterator import runGen
from collections import OrderedDict
from hashlib import md5
from inspect import getfullargspec

#TODO implement differenciation between worker and others with mpi communicators and broadcast (dict of workers)
#TODO implement consistent manager for tags
#TODO implement rather a class than a simple func?
#TODO make launch using non blocking MPI functions in way to be able to do anything else while not computing for the run to go on

def launch(runparameters, workfunc=lambda *args, **kwargs: None, opts={}):
    defopts = {'stdbehaviour': True, 'specialbehaviour': workfunc, 'specialbehaviour_data': {}}#, 'workfunc_needrawdata': False}
    defopts.update(opts)
    opts = defopts

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    #if is std dict then sort it by key so that rp becomes identificable
    if not isinstance(runparameters, OrderedDict):
        runparameters = OrderedDict(sorted(iteritems(runparameters), key=lambda t: t[0]))
    runparameterkeys = tuple(k for k in iterkeys(runparameters))

    #if it's a master
    if rank == 0:
        runparameterbounds = tuple(len(list(iterable)) for iterable in itervalues(runparameters))

        if 'runprogressfile' in opts:
            runprogressfile = opts['runprogressfile']
        else:
            invariant = str(runparameterkeys) + str(runparameterbounds)
            runprogressfile = "currentrun_" + md5(invariant.encode()).hexdigest()
        if isfile(runprogressfile):
            lrnf = eval(open(runprogressfile, 'r').read())#last run not finished
        else:
            lrnf = (0,) * len(runparameters)

        runsqueue = {}

        for it in runGen(runparameterbounds, lrnf):
            runner = comm.recv(tag=0, source=MPI.ANY_SOURCE)
            comm.send(it, dest=runner, tag=1)
#           run_data = {runparameterkeys[i]: it[i] for i in range(len(runparameters))}
#           comm.send(run_data, dest=runner, tag=1)

            runsqueue[runner] = it
            if lrnf < min(itervalues(runsqueue)):
                lrnf = min(itervalues(runsqueue))
                open(runprogressfile, 'w').write(str(lrnf))

        #manage the process still running after the end of the running loop
        while True:
            runner = comm.recv(tag=0, source=MPI.ANY_SOURCE)
            comm.send(False, dest=runner, tag=1)

            runsqueue.pop(runner, None)
            if not bool(runsqueue):#i.e. if runsqueue empty
                open(runprogressfile, 'w').write(str(runparameterbounds))#in case algo stops right now.. (paranoic)
                break
            if lrnf < min(itervalues(runsqueue)):
                lrnf = min(itervalues(runsqueue))
                open(runprogressfile, 'w').write(str(lrnf))

        if 'rpfbackup' in opts and bool(opts['rpfbackup']):
            if isinstance(opts['rpfbackup'], str):
                rpfileback = opts['rpfbackup']
            else:
                rpfileback = runprogressfile
                if isfile(runprogressfile + '.bak'):
                    nbSaves = 2
                    while isfile(runprogressfile + '_' + str(nbSaves) + '.bak'):
                        nbSaves += 1
                    rpfileback += '_' + str(nbSaves)
                    print("You have %d old backups. You may considerer deleting it!" % (nbSaves-1))
                rpfileback += '.bak'
            rename(runprogressfile, rpfileback)
        else:
            remove(runprogressfile)

    #otherwise is a worker
    else:
        #TODO: Maybe recconsider spliting on start but rather dynamicly as e.g. mpi brain decide the role of each one. NOTE: must then also modify rank 0 behaviour (design for instance a message passed by the worker to signal if he "gets some rest or continue working")
        if callable(opts['stdbehaviour']):
            behavestd = opts['stdbehaviour']()
        else:
            behavestd = opts['stdbehaviour']

        if behavestd:
            while True:
                comm.send(rank, dest=0, tag=0)
                rawdata = comm.recv(source=0, tag=1)
#           run_data = {runparameterkeys[i]: it[i] for i in range(len(runparameters))}
                if not rawdata:
                    break
                #data = {k[:-1]: runparameters[k][v] for k, v in iteritems(data)}
                data = {runparameterkeys[i][:-1]: runparameters[runparameterkeys[i]][rawdata[i]] for i in range(len(rawdata))}
                if 'rawdata' in getfullargspec(workfunc)[0]:
                    workfunc(**data, rawdata=rawdata)
                else:
                    workfunc(**data)
        else:
            #you may consider putting special behaviours in this func
            if opts['specialbehaviour_data']:
                opts['specialbehaviour'](**opts['specialbehaviour_data'])
            else:
                opts['specialbehaviour']()


