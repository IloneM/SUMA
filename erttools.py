from six import iteritems,iterkeys,itervalues
from inspect import getfullargspec
from runiterator import runGen
from mpi4py import MPI
import cma
import copy
import collections

class ERTLogger(cma.BaseDataLogger):
    def __init__(self, precs, inopts={'intrusivestop': True, '-log_precs': True}):
        self.precs = copy.copy(precs);
        if inopts['-log_precs']:
            print(precs)
            self.workprecs = [10 ** (-prec) for prec in precs]
        else:
            self.workprecs = copy.copy(precs)
        self.lastfcteval = 0
        self.currentprecit = 0
        self.resultsforert = {}
        self.nbprecs = len(precs)
        self.opts = copy.copy(inopts)

    def add(self, es):
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
        defopts = {'path': '.', 'splitkeys': {}, 'testkey': 'test', 'filename': 'results', 'fileext': 'dat', 'keysplural': True, 'data_tag': 2}
        defopts.update(inopts)
        self.opts = opts = defopts

        self.comm = mpicomm
        
        if opts['keysplural']:
            runparameters = {k[:-1]: runparameters[k] for k in iterkeys(runparameters)}

        keys = iterkeys(runparameters)

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

        if isinstance(runparameters[opts['testkey']], collections.Iterable):
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
        header += " ert\n"

        for fileparams in runGen(tuple(len(self.iterparams[splitkeyid]) for splitkeyid in self.splitkeyids)):
            with open(self.__getfilename(fileparams), 'w') as resfilestream:
                resfilestream.write(header)

    def __getfilename(self, fileparams):
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
        datastr = ''
        for paramid in range(len(runid)):
            if paramid == self.splitkeyids[splitid]:
                fileparams.append(runid[paramid])
            else:
                datastr += str(self.iterparams[paramid][runid[paramid]]) + ' '
        datastr += str(ert) + "\n"
        with open(self.__getfilename(fileparams), 'w') as resfilestream:
            resfilestream.write(datastr)

    def __call__(self):
        #TODO: use seprated comm for this work
        storeddata = dict()

        #find better solution than computing all the amount of data expected
        for fakeit in runGen(tuple(len(iterparam) for iterparam in self.iterparams)):
            data = list(self.comm.recv(tag=self.opts['MPI']['tag'], source=self.opts['MPI']['workers']))

            runid = data[1]
            del runid[self.testkeyid]
            runid = tuple(runid)

            if runid not in storeddata:
                storeddata[runid] = set()
            storeddata[runid].add(data[0])
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
        #TODO: create specific MPI comm
        defopts = {'MPI': {'tag': 2, 'manager': 1, 'workers': MPI.ANY_SOURCE}}
        defopts.update(inopts)
        self.opts = opts = defopts

        self.__manager = None

    #TODO: define rather whoami?
    def __amimanager():
        if self.__manager is None:
            self.__manager = self.rank in list(self.opts['MPI']['manager'])
        return self.__manager

    def __call__(**kwargs):
        if self.__amimanager():


    #idea: send str which is the key to look for in the dict sent as "rawdata" and then store data in splited files (splited according to a list of keys used for splitting e.g.) and then store the others as in compute_ert and use a specific path and put all this in the opt dict used in stdrun main funct and modify stdrun in way to send this stored datas to the stdbehaviour
