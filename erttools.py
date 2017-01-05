from mpi4py import MPI
#import suma
#import cma
import copy

class ERTManagerForWorkers:
    def __init__(self, runparams, precs, inopts={'intrusivestop': True}):
        self.cr = currentRun#must be dict as {'k': 0.3, test: 3, ..}
        self.precs = copy.copy(precs);
        self.workprecs = [10 ** (-prec) for prec in precs]
        self.lastfcteval = 0
        self.currentprecit = 0
        self.resultsforert = {}
        self.nbprecs = len(precs)
        self.opts = copy.copy(opts)

    def storedata(self, es):
        newf = es.fit.fit[0]
        for precit in range(self.currentprecit, self.nbprecs):
            if newf < self.workprecs[precit]:
                self.resultsforert[self.precs[precit]] = self.lastfcteval
                self.currentprecit += 1

        self.lastfcteval = es.countevals
        if self.opts['intrusivestop'] and self.currentprecit == self.nbprecs:
            es.opts['termination_callback'] = lambda x: True

    def __call__(self, es):
        self.storedata(es)

    def done(self, es):
        return self.currentprecit == self.nbprecs

    def sendresults(self):
        return self.resultsforert
        #comm.send({'func': func.__name__, 'k': k, 'd': dim, 'erts': ERT}, dest=1, tag=2)

