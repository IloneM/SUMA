from six import iteritems,iterkeys,itervalues
from inspect import getfullargspec
from runiterator import runGen
from mpi4py import MPI
from collections import OrderedDict,Iterable
import cma
import copy
import stdrun

#TODO rewrite it with special processing for test and self management of the nbtest with own processing via dedicated MPI comm/group
#WHY? Because test management as part of the runparameters is tricky (see below..) and inconvenient as the data is stored as computation of all results throughout nbtest runs, then all the recovery stuffs designed in stdrun are useless/dangerous because not thinked for that specific usage
#TODO integrate also better management of precs

class ERTLogger(cma.BaseDataLogger):
    def __init__(self, precs, inopts={'intrusivestop': True, '-log_precs': True}):
        self.precs = copy.copy(precs);
        if inopts['-log_precs']:
            self.workprecs = [10 ** (-prec) for prec in precs]
        else:
            self.workprecs = copy.copy(precs)
        self.lastfcteval = 0
        self.currentprecit = 0
        self.resultsforert = {}
        self.nbprecs = len(precs)
        self.opts = copy.copy(inopts)
        self.modulo = 1#the calculations are so thin that there is no need for logging only part of iter but creating var prevents some 'errors'

    def add(self, es, modulo=None):#modulo param in way to avoid errors as specified above
        newf = es.best.f
        for precit in range(self.currentprecit, self.nbprecs):
            if newf < self.workprecs[precit]:
                self.resultsforert[self.precs[precit]] = self.lastfcteval
                self.currentprecit += 1

        self.resultsforert['last'] = self.lastfcteval = es.countevals
        if self.opts['intrusivestop'] and self.currentprecit == self.nbprecs:
            es.opts['termination_callback'] = lambda x: True

    def data(self):
        return self.resultsforert

class ERTMPIManager:
    def __init__(self, runparameters, precs, mpicomm, inopts=dict()):
        """
        runparameters must be an OrderedDict so that keys are names of parameters as they appear in file.s and each key must point to a list of elements used to identify the run (as func name or whatever)
        """
        defopts = {'path': '.', 'splitkeys': {}, 'testkey': 'test', 'filename': 'results', 'fileext': 'dat', 'keysplural': True}
        defopts.update(inopts)
        self.opts = opts = defopts

        self.comm = mpicomm
        
        assert isinstance(runparameters, OrderedDict)
        if opts['keysplural']:
            runparameters = OrderedDict((k[:-1], runparameters[k]) for k in iterkeys(runparameters))
            opts['splitkeys'] = {k[:-1] for k in self.opts['splitkeys']}#also modify self.opts beacause same ref
        else:#avoid runparemeters to be modified globally which is not wanted
            runparameters = copy.copy(runparameters)

        keys = list(iterkeys(runparameters))

        self.testkeyid = None
        self.splitkeyids = []
        ##maybe used if we choose to write a func for computing header before writing data; see below
        ##self.infilekeyids = []
        infilekeys = []

        for keyid in range(len(keys)):
            if keys[keyid] == opts['testkey']:
                self.testkeyid = keyid
            elif keys[keyid] in opts['splitkeys']:
                self.splitkeyids.append(keyid - 0 if self.testkeyid is None else 1)
            else:
                infilekeys.append(keys[keyid])
                ##self.iterkeyids.append(keyid - 0 if self.testkeyid is None else 1)

        if isinstance(runparameters[opts['testkey']], Iterable):
            self.nbtest = len(runparameters.pop(opts['testkey']))
        else:
            self.nbtest = runparameters.pop(opts['testkey'])
        
        #Not useful
        #self.splitkeys = opts['splitkeys']
        self.iterparams = [list(iterable) for iterable in itervalues(runparameters)]
#        self.splitparams = [list( iterable) for k,iterable in iteritems(runparameters) if k in self.splitkeys]

        self.precs = precs

#        infilekeys = [k for k in keys if k not in self.splitkeys and not k == opts['testkey']]
        ##self.__initfiles([keys[infilekeyid] for infilekeyid in self.infilekeyids])
        self.__initfiles(infilekeys)

    @staticmethod
    def computeERT(precs, data):
        """
        precs are the precisions for ERTs
        data is a 1-d array of dict which are designed as {'precision1': fcteval, 'precision2': fcteval,..,'last':lasfcteval}
        where precisionN must be in precs but it is not neccessary to last precs to be in the dict
    """
        RTs, RTus, nbS, nbUS = tuple({prec: 0 for prec in precs} for fakeit in range(4))
        nbtest = len(data)

        for exp in data:
            lastfcteval = exp.pop('last')
            remainingprecs = set(precs)
            for prec,fcteval in iteritems(exp):
                nbS[prec] += 1
                RTs[prec] += fcteval
                remainingprecs.discard(prec)
            for remainingprec in remainingprecs:
                nbUS[remainingprec] += 1
                RTus[remainingprec] += lastfcteval

        ERT = {prec: 0 for prec in precs}
        for prec in precs:
            assert nbS[prec] + nbUS[prec] == nbtest
            if nbS[prec] > 0:
                RTs[prec]  /= nbS[prec]
            if nbUS[prec] > 0:
                RTus[prec] /= nbUS[prec]
            ps = nbS[prec] / nbtest
            if ps > 0:
                ERT[prec] = RTs[prec] + (1 - ps) / ps * RTus[prec]
            else:
                ERT[prec] = float("inf")

        return ERT

    def __initfiles(self, infilekeys):
        header = '#'
        for key in infilekeys:
            header += ' ' + str(key)
        header += " prec ert\n"

        if self.splitkeyids:
            for fileparams in runGen(tuple(len(self.iterparams[splitkeyid]) for splitkeyid in self.splitkeyids)):
                with open(self.__getfilename(fileparams), 'w') as resfilestream:
                    resfilestream.write(header)
        else:#empty i.e. one single file
            with open(self.__getfilename(), 'w') as resfilestream:
                resfilestream.write(header)

    def __getfilename(self, fileparams=tuple()):
        suffix = ''
        for paramid in range(len(fileparams)):
            suffix += '_' + str(self.iterparams[self.splitkeyids[paramid]][fileparams[paramid]])
        return self.opts['path'] + self.opts['filename'] + suffix + '.' + self.opts['fileext']


    def __storeERT(self, runid, ert):
        """
        runparameters must be an OrderedDict so that keys are names of parameters as they appear in file.s and each key must point to a list of elements used to identify the run (as func name or whatever)
        """
        splitid = 0
        fileparams = []
        ##Maybe write a func for creating data header? Then create infilekeyids as commented above in __init__
        header = ''

        if self.splitkeyids:
            lensplitkeyids = len(self.splitkeyids)
            for paramid in range(len(runid)):
                if lensplitkeyids > splitid and paramid == self.splitkeyids[splitid]:
                    fileparams.append(runid[paramid])
                    splitid += 1
                else:
                    header += str(self.iterparams[paramid][runid[paramid]]) + ' '
        else:
            for paramid in range(len(runid)):
                header += str(self.iterparams[paramid][runid[paramid]]) + ' '
        datastr = ''
        for prec,ert in iteritems(ert):
            datastr += header + ' ' + str(prec) + ' ' + str(ert) + "\n"
        #datastr += str(ert) + "\n"
        with open(self.__getfilename(fileparams), 'a') as resfilestream:
            resfilestream.write(datastr)

    def __call__(self):
        #TODO: use seprated comm for this work
        storeddata = dict()

        #find better solution than computing all the amount of data expected (above all with this ugly writing)
        for fakeit in runGen([len(iterparam) for iterparam in self.iterparams] + [self.nbtest]):
            data = list(self.comm.recv(tag=self.opts['MPI']['tag'], source=self.opts['MPI']['workers']))

            runid = list(data[1])
            del runid[self.testkeyid]
            runid = tuple(runid)

            if runid not in storeddata:
                #set cannot be used because dict (type of data[0]) aren't hashable
                storeddata[runid] = list()
            #storeddata[runid].add(data[0])
            storeddata[runid].append(data[0])
            if len(storeddata[runid]) >= self.nbtest:
                self.__storeERT(runid, self.computeERT(self.precs, storeddata.pop(runid)))

class ERTMPIWorker:
    def __init__(self, precs, workfunc, mpicomm, inopts):
        defopts = dict()
        defopts.update(inopts)
        self.opts = opts = defopts

        self.precs = precs
        self.workfunc = workfunc
        self.comm = mpicomm

    def __senddata(self, ertdata, rawrunit):
        self.comm.send((ertdata, rawrunit), dest=self.opts['MPI']['manager'], tag=self.opts['MPI']['tag'])

    #consider storring logger and write a reset func
    def __call__(self, rawdata, **kwargs):
        if 'logger' in getfullargspec(self.workfunc)[0]:
            kwargs['logger'] = logger = ERTLogger(self.precs)
            self.workfunc(**kwargs)
        else:
            logger = self.workfunc(**kwargs)
        self.__senddata(logger.data(), rawdata)
        #TODO: Maybe add suport for cases where value returned by func is the dict as in logger.data

class ERTMPIWrapper:
    def __init__(self, runparameters, workfunc, precs, inopts):
        #TODO: create specific MPI comm or group
        defopts = {'MPI': {'tag': 2, 'manager': 1, 'workers': MPI.ANY_SOURCE}}
        defopts.update(inopts)
        self.opts = opts = defopts

        #if is std dict then sort it by key so that rp becomes identificable
        if not isinstance(runparameters, OrderedDict):
            self.runparameters = OrderedDict(sorted(iteritems(runparameters), key=lambda t: t[0]))
        else:
            self.runparameters = runparameters
        self.precs = precs
        self.workfunc = workfunc

        #TODO: manage MPI specific stuffs here
        comm = MPI.COMM_WORLD
        self.rank = comm.Get_rank()

#        self.__manager = None
        if isinstance(opts['MPI']['manager'], Iterable):
            self.__manager = self.rank in list(opts['MPI']['manager'])
        else:
            self.__manager = self.rank == opts['MPI']['manager']
        self.rank = None
#        self.workinstance = None

        #TODO add case master i.e. process sending tasks to everyone where all this stuff is not needed
        #TODO find user friendler names than 'master', 'worker', 'manager'
        if self.__manager:
            #TODO: Put this inside manager constructor?
            manager_runparameters = copy.copy(self.runparameters)
            if 'manager' in self.opts:
                opts = self.opts['manager']
                #if the parameters for storing data are different from those used for work (e.g. function name instead of the actual function) they must be specified in this option. NOTE: the update of params is only available on lists i.e. you cannot update {'f': [foo, bar]} with {'f': [foo.__name__]} but with {'f': [foo.__name__, bar.__name__]} e.g. 
                manager_runparameters.update(opts.pop('specific_runparameters', dict()))
            else:
                opts = dict()
            opts['MPI'] = self.opts['MPI']

            self.workinstance = ERTMPIManager(manager_runparameters, self.precs, comm, opts)

        else:
            #TODO: Put this inside worker constructor?
            if 'worker' in self.opts:
                opts = self.opts['worker']
            else:
                opts = dict()
            opts['MPI'] = self.opts['MPI']

            self.workinstance = ERTMPIWorker(self.precs, self.workfunc, comm, opts)

    #TODO: add support for direct args feeding here
    def launch(self):
        #TODO take a look at optimization of useless operations (if they are..) executed each time a new rawdata is given
        stdrun.launch(self.runparameters, self.workinstance, {'stdbehaviour': not self.amimanager()})


    #TODO: define rather whoami?
    def amimanager(self):
        #if self.__manager is None:
        #    self.__manager = self.rank in list(self.opts['MPI']['manager'])
        return self.__manager

    #def __call__(self, **kwargs):
    #    #TODO: manage MPI specific stuffs here
    #    comm = MPI.COMM_WORLD
    #    self.rank = comm.Get_rank()

    #    if self.workinstance is None:
    #        if self.__amimanager():
    #            #TODO: Put this inside worker constructor?
    #            if 'worker' in self.opts:
    #                opts = self.opts['worker']
    #            else:
    #                opts = dict()
    #            opts['MPI'] = self.opts['MPI']

    #            self.workinstance = ERTMPIManager(self.precs, self.workfunc, comm, opts)

    #        else:
    #            #TODO: Put this inside manager constructor?
    #            if 'manager' in self.opts:
    #                opts = self.opts['manager']
    #            else:
    #                opts = dict()
    #            opts['MPI'] = self.opts['MPI']

    #            self.workinstance = ERTMPIManager(self.runparameters, self.precs, comm, opts)

    #idea: send str which is the key to look for in the dict sent as "rawdata" and then store data in splited files (splited according to a list of keys used for splitting e.g.) and then store the others as in compute_ert and use a specific path and put all this in the opt dict used in stdrun main funct and modify stdrun in way to send this stored datas to the stdbehaviour
