import cma
import copy
from six import iteritems

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

