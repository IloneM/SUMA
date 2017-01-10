#from mpi4py import MPI
#import suma
import cma
import copy
from six import iteritems

##class ERTManagerForWorkers:
#class ERTWrapper:
class ERTLogger(cma.BaseDataLogger):
##    def __init__(self, precs, inopts={'intrusivestop': True, '-log_precs': True}):
    def __init__(self, precs, inopts={'intrusivestop': True, '-log_precs': True}):
##        self.cr = currentRun#must be dict as {'k': 0.3, test: 3, ..}
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
##        self.runparams = runparams

#    def add(self, es=None, more_data=[], modulo=None):
    def add(self, es):
#    def storedata(self, es):
        newf = es.best.f
        #newf = es.fit.fit[0]
        for precit in range(self.currentprecit, self.nbprecs):
            if newf < self.workprecs[precit]:
                self.resultsforert[self.precs[precit]] = self.lastfcteval
                self.currentprecit += 1

        self.resultsforert['last'] = self.lastfcteval = es.countevals
        if self.opts['intrusivestop'] and self.currentprecit == self.nbprecs:
#        if self.currentprecit == self.nbprecs:
            es.opts['termination_callback'] = lambda x: True

#   def __call__(self, es, objective_fct, args=(), call_back=[], iterations=None, min_iterations=1):
#       #es.opts['termination_callback'] = lambda x: True
#       if isinstance(call_back, list):
#           call_back.append(self.storedata)
#       else:
#           call_back = [call_back, self.storedata]
#       es.optimize(objective_fct, args=args, call_back=call_back,
#                   iterations=iterations,min_iterations=min_iterations)

#    def done(self, es):
#        return self.currentprecit == self.nbprecs

#    def getresults():
    def data(self):
        return self.resultsforert
##    def sendresults(self, comm, dest=1, tag=10):
##        comm.send((self.runparams, self.resultsforert), dest=dest, tag=tag)
#        return self.resultsforert

#class ERTLoggerMPI(ERTLogger):
#    def senddata(comm, dest, tag, more_data=None):
#        comm.send((self.resultsforert,more_data) if more_data else self.resultsforert, dest, tag)

def computeERT(precs, data):
    """
    precs are the precisions for ERTs
    data is a 1-d array of dict which are designed as {'precision1': fcteval, 'precision2': fcteval,..,'last':lasfcteval}
    where precisionN must be in precs but it is not neccessary to last precs to be in the dict
"""
    RTs, RTus, nbS, nbUS = tuple({prec: 0 for prec in precs} for fakeit in range(4))
    nbtest = len(data)

    for exp in data:
#        precit = 0
        lastfcteval = exp.pop('last')
        remainingprecs = set(precs)
        for prec,fcteval in iteritems(exp):
            nbS[prec] += 1
            RTs[prec] += fcteval
            remainingprecs.discard(prec)
#            precit += 1
#            if precit >= nbprecs:
#                break
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
#    comm.send({'func': func.__name__, 'k': k, 'd': dim, 'erts': ERT}, dest=1, tag=2)

#def sendERTMPI()

#def storeERTMPI(rank):


#class ERTMPIWorker():
#    def __init__():

#class ERTManagerFor:
