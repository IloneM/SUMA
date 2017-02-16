import numpy as np
import cma
import copy
from math import log,exp

class SumaSmall():
#TODO check whether x0 is copy of obj or not
    def __init__(self, x0, sigma0, inopts, suma_inopts, bigstrat):
        #inopts['CMA_diagonal'] = False
        #inopts['AdaptSigma'] = False
        #inopts['CMA_eigenmethod'] = cma.Misc.eig
        
        #inopts['verb_filenameprefix'] = inopts['verb_filenameprefix'] + 'Small';

#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
#TODO check if param exist
        suma_opts = suma_inopts;
        self.nrp = max(suma_opts['SUMA_nrp'],2)
        self.mean = []
        self.Csi = [np.zeros(shape=(suma_opts['SUMA_d'],suma_opts['SUMA_D'])) for fakit in range(self.nrp)]

        cma.rotate(x0)
        self._random_projection = [cma.rotate.dicMatrices[str(suma_opts['SUMA_D'])][-suma_opts['SUMA_d']:,:],
                                   cma.rotate.dicMatrices[str(suma_opts['SUMA_D'])][:suma_opts['SUMA_d'],:]]
        self.mean.append(self._random_projection[0] @ x0)
        self.mean.append(self._random_projection[1] @ x0)
        for i in range(2,self.nrp):
            self._random_projection.append(self.grp(suma_opts['SUMA_d'], suma_opts['SUMA_D']))
            self.mean.append(self._random_projection[i] @ x0)

        self.bs = bigstrat
        #self.popsize = bigstrat.popsize
        #self.mu = bigstrat.mu

        
        #super().__init__(self._random_projection @ x0, sigma0, inopts);

#TODO does not work; must be before calling super so check if exists in inopts an if not take default intead plus Small
#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
        #self.opts['verb_filenameprefix'] = self.opts['verb_filenameprefix'] + 'Small';

        #self.inopts = {**self.inopts, **suma_inopts};
        self.inopts = inopts
        #self.opts = {**self.opts, **suma_opts};

        #test 
        #self._first_identity = True;
        #self.gp.isidentity = False

    @staticmethod
    def grp(d, D, seed=None):
        rstate = np.random.get_state()
        np.random.seed(seed) if seed else np.random.seed()
        B = np.random.randn(d, D)
        for i in range(d):
            for j in range(0, i):
                B[i] -= np.dot(B[i], B[j]) * B[j]
            B[i] /= sum(B[i]**2)**0.5

        np.random.set_state(rstate)
        return B
        ##res = np.zeros(shape=(d,D));
        ##i = 0;
        ##while i < d:
        ##    rand_vector = new_basis_candidate = np.zeros(shape=(D), dtype=float);
        ##    for j in range(D):
        ##        if np.random.uniform() > spc:
        ##            rand_vector[j] = new_basis_candidate[j] = 2*np.random.randint(1)-1;
        ##    j=0;
        ##    for j in range(i):
        ##       new_basis_candidate -= res[j] * np.dot(rand_vector,res[j]);
        ##       #we don't need to divide by result.row(j).norm() cause all vectors are normalized
        ##    if np.linalg.norm(new_basis_candidate) > eps:
        ##        res[i] = new_basis_candidate / np.linalg.norm(new_basis_candidate);
        ##        i += 1;
        ##return res;
    
    def relativeProba(self, X, rpi):
        #renormalize in way to avoid centering effect
#        Px = self._random_projection @ X - self.mean;
        Px = self._random_projection[rpi] @ X - self.mean[rpi];
        if(np.linalg.norm(Px) > 0):
            Csx = self.Csi[rpi] @ (Px / np.linalg.norm(Px));
            return -Csx.dot(Csx);
        else:
            return 0.0;
        #Csx = self.Csi @ (self._random_projection @ X - self.mean);
        ##the probability evoluates as -||Csx||^2 in a guassian distrib because of exp
        #return -Csx.dot(Csx);

    def perf(self, sortedcand, fit, relevance):
        #if not self.perfdoneonce:
        #    slef.fmin = min(fit)
        #    self.fmax = max(fit)
        #    self.dminmax = self.fmax-self.fmin
#       #     self.popsize = len(fit)
        popsize = len(fit)

        rel = relevance(fit, popsize)

        res = []
        for rpi in range(self.nrp):
            if self.setneword:
                newordt = sorted([(i, self.relativeProba(sortedcand[i], rpi)) for i in range(popsize)], key=lambda t: t[1], reverse=True)
                self.neword.append([t[0] for t in newordt])
            neword = self.neword[rpi]
            #fitord = fit[neword]
            
            perfval = 0
            cumuldistr = 0
            cumulf = 0
            for i in range(popsize):
                cumuldistr += rel[neword[i]]/log(i+2, 2)
                cumulf += rel[i]/log(i+2, 2)
                #cumuldistr += (2 ** (fmin/fitord[i]) -1)/log(i+2, 2)
                #cumulf += (2 ** (fmin/fit[i]) -1)/log(i+2, 2)
                if cumulf > 0:
                    perfval += cumuldistr / cumulf
                #perfval += cumuldistr / cumulf
            res.append(perfval / popsize)
        if self.setneword:
            self.setneword = False
        return res

    def begin_iteration(self):
#        self._flgtelldone = False;
#        self.evaluations_per_f_value = 1
#        self._tmp_candidates = np.transpose(self._random_projection @ np.transpose(big_candidates));
        #if all([self.D[i] == 1.0 for i in range(self.N)]) and not self._first_identity:
        #    raise cma._Error('identity again');
        #else:
        #    self._first_identity = False;
        #    self._savedC = self.C;

        self.perfdoneonce = False
        self.neword = []
        self.setneword = True
        for i in range(self.nrp):
            self.Csi[i] = self._random_projection[i] @ self.bs.Csi @ self._random_projection[i].T
            self.mean[i] = self._random_projection[i] @ self.bs.mean
        #self.Csi = np.dot(self.B, (self.B / self.D).T)  # square root of inverse of C according to cma.py update_exponential

#    def end_iteration(self):
#            self._flgtelldone = False

class Suma(cma.CMAEvolutionStrategy):
#TODO check whether x0 is copy of obj or not
    def __init__(self, x0, sigma0, inopts={}):
        suma_inopts = self.extract_suma_inopts(inopts);
#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
#TODO check if param exist
        #if 'SUMA_D' not in inopts:
        suma_inopts['SUMA_D'] = len(x0);
        #else if inopts['SUMA_D'] != len(x0);
        #    raise _error('suma_d has to match  = ' + str(n)) 
        suma_opts = suma_inopts;

#        self._smallstrat_inopts = copy.copy(inopts);
        #EDIT we do modify in SumaSmall __init__()
        #we do not need to copy inopts because it is not modified by CMAEvolutionStrategy.__init__()
#        self._smallstrat = SumaSmall(x0, sigma0, copy.copy(inopts), suma_inopts);
#        inopts['CMA_diagonal'] = True;

        smallinopts = copy.copy(inopts)

        super().__init__(x0, sigma0, inopts);

        self._smallstrat = SumaSmall(x0, 1.0, smallinopts, suma_inopts, self);

#        self.Pqualbound = (self.N / suma_opts['SUMA_d']) ** 0.5
        self.mu = self.sp.mu
        
        self._f = None
        if 'SUMA_f' in suma_opts:
            self._f = suma_opts['SUMA_f']

        self.inopts = {**self.inopts, **suma_inopts};
        self.opts = {**self.opts, **suma_opts};

    def extract_suma_inopts(self, inopts):
        suma_keys = [k for k in inopts if k.startswith('SUMA')];
        suma_inopts = {};
        for key in suma_keys:
            suma_inopts[key] = inopts.pop(key, None);
        return suma_inopts;

#    def optimize(self, objective_fct, iterations=None, min_iterations=1, args=(), verb_disp=None, logger=None, call_back=None):
#        self._smallstrat.logger.add(self._smallstrat);
#        super().optimize(objective_fct, iterations, min_iterations, args, verb_disp, logger, call_back);


#    def ask(self, number=None, xmean=None, sigma_fac=1, gradf=None, args=()):
    def ask_geno(self, number=None, xmean=None, sigma_fac=1):
        if number is None or number < 1: 
            number = self.popsize;
        
        #pop = np.matrix();
        pop = np.ndarray(shape=(number, self.N));
        self.Csi = np.dot((self.B / self.D), self.B.T) / self.sigma_vec / self.sigma
        self._smallstrat.begin_iteration()
        self.neword = None

###        for i in range(number):
####TODO see EDIT above
####TODO check if param exist
###            candidates = super().ask_geno(self.opts['SUMA_k']);
###            pop[i] = self._smallstrat.bestCandidate(candidates);

        #return pop;
        return super().ask_geno(number,xmean,sigma_fac)

    def relproba(self, X, isproba=True):
        #renormalize in way to avoid centering effect
#        Px = self._random_projection @ X - self.mean;
        Px = X - self.mean
        if(np.linalg.norm(Px) > 0):
            #Csx = self.Csi @ (Px / np.linalg.norm(Px));
            Csx = self.Csi @ Px
            if isproba:
            	return exp(-Csx.dot(Csx));
            else:
            	return -Csx.dot(Csx);
        else:
            #return 1.;
            if isproba:
            	return 1.
            else:
            	return 0.
        #Csx = self.Csi @ (self._random_projection @ X - self.mean);
        ##the probability evoluates as -||Csx||^2 in a guassian distrib because of exp
        #return -Csx.dot(Csx);

    def relevancearith(self, fit, popsize):
        return [(2 ** ((self.fmax-fi)/self.dminmax) -1) for fi in fit]

    def relevancegeom(self, fit, popsize):
        return [(2 ** (self.fmin/fi) -1) for fi in fit]

    @staticmethod
    def relevancefromproba(fit, popsize):
        return [(2 ** fi -1) for fi in fit]

    @staticmethod
    def relevanceordersoft(fit, popsize):
        return [(1/(i+1)) for i in range(popsize)]

    @staticmethod
    def relevanceorderhard(fit, popsize):
        return [2 ** (-i) for i in range(popsize)]

    @staticmethod
    def relevanceorderexpsoft(fit, popsize):
        return [(1/(i+1)) for i in range(popsize)]

    @staticmethod
    def relevanceorderexphard(fit, popsize):
        return [2 ** reli -1 for reli in Suma.relevanceorderhard(None, popsize)]

    @staticmethod
    def relevanceorderexpsoft(fit, popsize):
        return [2 ** reli -1 for reli in Suma.relevanceordersoft(None, popsize)]

    @staticmethod
    def relevancemuonly(fit, popsize):
        return [(1 if i < popsize // 2 else 0) for i in range(popsize)]

    def perf(self, sortedcand, fit, relevance):
        #fmin = min(fit)
        #fmax = max(fit)
        #dminmax = fmax-fmin
        popsize = len(fit)

        if self.neword is None:
            self.newordt = sorted([(i, self.relproba(sortedcand[i])) for i in range(popsize)], key=lambda t: t[1], reverse=True)
            self.neword = [t[0] for t in self.newordt]
        neword = self.neword
        #fitord = fit[neword]
        
        rel = relevance(fit, popsize)
        perfval = 0
        cumuldistr = 0
        cumulf = 0

        for i in range(popsize):
            #cumuldistr += (2 ** (fmin/fitord[i]) -1)/log(i+2, 2)
            #cumulf += (2 ** (fmin/fit[i]) -1)/log(i+2, 2)
            cumuldistr += rel[neword[i]]/log(i+2, 2)
            cumulf += rel[i]/log(i+2, 2)
            if cumulf > 0:
                perfval += cumuldistr / cumulf
        return perfval / popsize

    def tell(self, solutions, function_values, check_points=None, copy=False):
#        Psolutions = np.transpose(self._smallstrat._random_projection @ np.transpose(solutions));
       # print(solutions)
       # for i in range(len(solutions)):
       #     print(Psolutions[i])
        #    print(self._smallstrat._random_projection.dot(solutions[i]))
        #self._smallstrat.tell(Psolutions, function_values, check_points, copy);
        #self._smallstrat.logger.add(self._smallstrat);
        #self._smallstrat.tell(self._smallstrat.ask(), function_values, check_points, copy);
        super().tell(solutions, function_values, check_points, copy);

        #mp = sum(self.pop_sorted[:self.mu])  / self.mu
        #mm = sum(self.pop_sorted[-self.mu:]) / self.mu
        #vmp = mp - mm
        #Pvmp = np.dot(self._smallstrat._random_projection, vmp)
        
        ##for f in [self.relevanceorderexphard, self.relevanceorderhard, self.relevanceordersoft, self.relevanceorderexpsoft, self.relevancemuonly]:
        ##    with open(self.opts['verb_filenameprefix'] + '_small_' + f.__name__ + '.dat', 'a') as resfilestream:
        ##        resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(self.pop_sorted, self.fit.fit, f)])+ '\n')
        ##    with open(self.opts['verb_filenameprefix'] + '_big_' + f.__name__ + '.dat', 'a') as resfilestream:
        ##        resfilestream.write(str(self.perf(self.pop_sorted, self.fit.fit, f))+ '\n')

        ##self._smallstrat.setneword = True
        ##self._smallstrat.neword = []

##        self.newordt = sorted([(i, self.relproba(self.pop_sorted[i])) for i in range(self.popsize)], key=lambda t: t[1], reverse=True)
##        self.neword = [t[0] for t in self.newordt]
##        newpop = self.pop_sorted[self.neword]
##        for f in [self.relevanceorderexphard, self.relevanceorderhard, self.relevanceordersoft, self.relevanceorderexpsoft, self.relevancemuonly, self.relevancefromproba]:
##            with open(self.opts['verb_filenameprefix'] + '_small_acc_big_' + f.__name__ + '.dat', 'a') as resfilestream:
##                resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(newpop, [t[1] for t in self.newordt], f)])+ '\n')

        #self.neword = None
        #self._smallstrat.setneword = True
        #self._smallstrat.neword = []
        
        #for f in [self.relevanceorderexphard, self.relevanceorderhard, self.relevanceordersoft, self.relevanceorderexpsoft]:
        #    with open(self.opts['verb_filenameprefix'] + '_small_onlymu_' + f.__name__ + '.dat', 'a') as resfilestream:
        #        resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(self.pop_sorted[:self.mu], self.fit.fit[:self.mu], f)])+ '\n')
        #    with open(self.opts['verb_filenameprefix'] + '_big_onlymu_' + f.__name__ + '.dat', 'a') as resfilestream:
        #        resfilestream.write(str(self.perf(self.pop_sorted[:self.mu], self.fit.fit[:self.mu], f))+ '\n')

        #self._smallstrat.setneword = True
        #self._smallstrat.neword = []
        #
        #for f in [self.relevancefromproba]:
        #    with open(self.opts['verb_filenameprefix'] + '_small_acc_big_onlymu_' + f.__name__ + '.dat', 'a') as resfilestream:
        #        resfilestream.write(str(self.perf(self.pop_sorted[self.neword][:self.mu], [t[1] for t in self.newordt][:self.mu], f))+ '\n')
        poptestN = 100
        poptest = np.random.randn(poptestN, self.N)
        for i in range(poptestN):
#            poptest[i] *= self.sigma
            poptest[i] += self.mean

        self.ordft = sorted([(i, self._f(poptest[i])) for i in range(poptestN)], key=lambda t: t[1])
        self.ordf = [t[0] for t in self.ordft]
        self.newordt = sorted([(i, self.relproba(poptest[i], False)) for i in range(poptestN)], key=lambda t: t[1], reverse=True)
#        self.newordt = sorted([(i, self.relproba(self.pop_sorted[i])) for i in range(self.mu)], key=lambda t: t[1])
        self.neword = [t[0] for t in self.newordt]

        self._smallstrat.setneword = True
        self._smallstrat.neword = []

        #newpop = self.pop_sorted[self.neword]
        newpop = poptest[self.neword]
        newfit = [t[1] for t in self.newordt]
        #if not any(newfit):
            #print(self.newordt)
            #print([np.linalg.norm(p-self.mean) for p in poptest])
            #print(self.sigma)
            #print([self.relproba(poptest[i], False) for i in range(poptestN)])
        self.fmin = newfit[-1]
        self.fmax = newfit[0]
        self.dminmax = self.fmax-self.fmin

        for f in (self.relevanceorderhard, self.relevanceordersoft, self.relevanceorderexphard, self.relevanceorderexpsoft, self.relevancemuonly):
            with open(self.opts['verb_filenameprefix'] + '_small_acc_big_' + f.__name__ + '.dat', 'a') as resfilestream:
                resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(newpop, newfit, f)])+ '\n')

        self._smallstrat.setneword = True
        self._smallstrat.neword = []

        #newpop = self.pop_sorted[self.neword]
        newpop = poptest[self.ordf]
        newfit = [t[1] for t in self.ordft]
        self.fmin = newfit[0]
        self.fmax = newfit[-1]
        self.dminmax = self.fmax-self.fmin

        #for f in [self.relevancearith, self.relevancegeom, self.relevanceorderexphard, self.relevanceorderexpsoft, self.relevancemuonly]:
        for f in (self.relevanceorderhard, self.relevanceordersoft, self.relevanceorderexphard, self.relevanceorderexpsoft, self.relevancemuonly):
            with open(self.opts['verb_filenameprefix'] + '_small_acc_f_' + f.__name__ + '.dat', 'a') as resfilestream:
                resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(newpop, newfit, f)])+ '\n')
            with open(self.opts['verb_filenameprefix'] + '_big_acc_f_' + f.__name__ + '.dat', 'a') as resfilestream:
                resfilestream.write(str(self.perf(newpop, newfit, f))+ '\n')

#        f = self.relevancefromproba
#
#        with open(self.opts['verb_filenameprefix'] + '_small_acc_big_onlymu.dat', 'a') as resfilestream:
#            resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(self.pop_sorted[self.neword], [t[1] for t in self.newordt], f)])+ '\n')
        
#        with open('smallfitaccbig.dat', 'a') as resfilestream:
#            resfilestream.write(' '.join([str(s) for s in self._smallstrat.perf(self.pop_sorted, self.fit.fit)])+ '\n')

    def optimize(self, iterations=None, min_iterations=1,
                 args=(), verb_disp=None, logger=None, call_back=None):
        super().optimize(self._f, iterations, min_iterations, args, verb_disp, logger, call_back)

    def optimize(self, objective_fct, iterations=None, min_iterations=1,
                 args=(), verb_disp=None, logger=None, call_back=None):
        if self._f is None:
            self._f = objective_fct
        super().optimize(self._f, iterations, min_iterations, args, verb_disp, logger, call_back)

