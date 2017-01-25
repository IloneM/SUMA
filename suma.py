import numpy as np
import owncma
import copy

# mandatory opts :
# (SUMA_d,) SUMA_spc (find if possible to choose by process), SUMA_a0 (see TODO below), SUMA_calphaeps, SUMA_k

class SumaSmall():
#TODO check whether x0 is copy of obj or not
    #def __init__(self, x0, inopts, suma_inopts):
    def __init__(self, x0, inopts):
#TODO write logger for this
#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
#TODO check if param exist
        #TODO: is it relevant?
        #self.opts = self.inopts = {**inopts, **suma_inopts}
        self.opts = self.inopts = inopts
##        self.opts = {**self.opts, **suma_opts}

        self.d = inopts['SUMA_d']
        self.rp = self.grp(self.d, inopts['SUMA_D'], inopts['SUMA_spc'])
        self.trp = np.transpose(self.rp)
        self.mean = self.rp @ x0
        #TODO: find better way manage triangulation of C; see below
        self.C = np.eye(self.d, dtype=float) * 0.5
        #TODO: must be found an fixed if possible
        self.alpha = inopts['SUMA_a0']
        self.calpha = 1. + inopts['SUMA_calphaeps']
        
#TODO does not work; must be before calling super so check if exists in inopts an if not take default instead plus Small
#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
    @staticmethod
    def grp(d, D, spc, eps=1e-6):
        res = np.zeros(shape=(d,D))
        i = 0
        while i < d:
            rand_vector = new_basis_candidate = np.zeros(shape=(D), dtype=float)
            for j in range(D):
                if np.random.uniform() > spc:
                    rand_vector[j] = new_basis_candidate[j] = 2*np.random.randint(1)-1
            j=0;
            for j in range(i):
               new_basis_candidate -= res[j] * np.dot(rand_vector,res[j])
               #we don't need to divide by result.row(j).norm() cause all vectors are normalized
            if np.linalg.norm(new_basis_candidate) > eps:
                res[i] = new_basis_candidate / np.linalg.norm(new_basis_candidate)
                i += 1
        return res;
    
    def regenrp(self):
        self.rp = self.grp(self.inopts['SUMA_d'], self.inopts['SUMA_D'], self.inopts['SUMA_spc'])
        self.trp = np.transpose(self.rp)
        
    #def logproba(self, X, C=None, mean=None):
    def logproba(self, X, projected=False, centered=False):
        #if C is None:
        #    C = self.C
        #if mean is None:
        #    mean = self.mean
        #Px = self.rp @ X - mean
        if not projected:
            X = self.rp @ X
        if not centered:
            X -= self.mean

        if([1 for coor in X if not coor == 0]):
            #renormalize in way to avoid centering effect (No because Michele said norm is part of the feature we are looking for)
            #TODO find better way to manage that C is inf triangular
            Cx = (self.C + self.C.T) @ X
            return -X.dot(Cx)
        else:
            return 0.0

    def foretellbest(self, candidates, projected=False, centered=False):
        if len(candidates) < 1:
            raise cma._Error('candidates vector cannot be empty');
        else:
            best_seen = candidates[0]
            score_bs = self.logproba(best_seen, projected, centered)
            for i in range(1,len(candidates)):
                scorei = self.logproba(candidates[i], projected, centered)
                if  scorei > score_bs:
                    best_seen = candidates[i]
                    score_bs = scorei
            return best_seen

    def fitsmalldist(self, solutions, projected=False, centered=False):
        fit = 0.
        for s in solutions:
            fit += self.logproba(s, projected, centered)
        return fit

    #receives big.pop_sorted[:self.sp.mu]
    def tell(self, solutions, mean):
        self.mean = self.rp @ mean
        #csols = np.array(sol - self.mean for sol in (solutions @ self.trp))
        csols = solutions @ self.trp - self.mean
#        sols = sols[ordersols[::-1]]
        self.U = np.zeros(shape=(self.d, self.d), dtype=float)
        for i in range(1,self.d):
            for j in range(1,i):
                self.U[i,j] = i * j * (csols[i] @ csols[j])
            self.U[i, i] = i**2 * (csols[i] @ csols[i]) / 2

        #NOTE: be careful: self.alpha must be in last pos so that C adapted with alpha is current alpha is in last pos in Ccnadidates so that the last update of self.C gets the correct value for this iteration
        alphacandidates = (min(self.alpha * self.calpha, 1), self.alpha / self.calpha, self.alpha)
        Ccandidates = tuple((1-alpha) * self.C + alpha * self.U for alpha in alphacandidates)
        #TODO: make in const i.e multiply by 3
        #res = [0.] * len(alphacandidates)
        res = []

        for C in Ccandidates:
            self.C = C 
            res.append(self.fitsmalldist(csols, projected=True, centered=True))

        #alpha plus is best
        if res[0] > res[2] and res[2] > res[1]:
            self.alpha = alphacandidates[0]
        #alpha minus is best
        elif res[1] > res[2] and res[2] > res[0]:
            self.alpha = alphacandidates[1]
        #else alpha remains unchanged

class Suma(owncma.CMAEvolutionStrategy):
#TODO check whether x0 is copy of obj or not
    def __init__(self, x0, sigma0, inopts={}):
        suma_inopts = self.extract_suma_inopts(inopts);
#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
#TODO check if param exist
        #if 'SUMA_D' not in inopts:
        inopts['CMA_diagonal'] = True
        #inopts['CMA_const_trace'] = True
        #inopts['CMA_on'] = True

        super().__init__(x0, sigma0, inopts)

        suma_inopts['SUMA_D'] = len(x0)
        suma_inopts['SUMA_d'] = self.sp.mu
        #suma_inopts['SUMA_d'] = self.popsize
        suma_opts = copy.copy(suma_inopts)
        #else if inopts['SUMA_D'] != len(x0);
        #    raise _error('suma_d has to match  = ' + str(n)) 

        #self._smallstrat = SumaSmall(x0, copy.copy(inopts), suma_inopts);
        self._smallstrat = SumaSmall(x0, suma_opts)

        self.inopts = {**self.inopts, **suma_inopts};
        self.opts = {**self.opts, **suma_opts};

    def extract_suma_inopts(self, inopts):
        suma_keys = [k for k in inopts if k.startswith('SUMA')];
        suma_inopts = {};
        for key in suma_keys:
            suma_inopts[key] = inopts.pop(key, None);
        return suma_inopts;

#    def ask(self, number=None, xmean=None, sigma_fac=1, gradf=None, args=()):
    def ask_geno(self, number=None, xmean=None, sigma_fac=1):
        if number is None or number < 1: 
            number = self.sp.popsize;
        pop = np.ndarray(shape=(number, self.N));

        for i in range(number):
#TODO see EDIT above
#TODO check if param exist
            candidates = super().ask_geno(self.opts['SUMA_k']);
            pop[i] = self._smallstrat.foretellbest(candidates);
        return pop;

    def tell(self, solutions, function_values, check_points=None, copy=False):
        self.C /= np.mean(self.C)
        self.dC = self.C
        super().tell(solutions, function_values, check_points, copy);
        self._smallstrat.tell(self.pop_sorted[:self.sp.mu], self.mean);
        #self._smallstrat.tell(self.pop_sorted, self.mean);

