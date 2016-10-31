import numpy as np
import cma
import copy

class SumaSmall(cma.CMAEvolutionStrategy):
#TODO check whether x0 is copy of obj or not
    def __init__(self, x0, sigma0, inopts, suma_inopts):
        inopts['CMA_diagonal'] = False
        inopts['AdaptSigma'] = False
        #inopts['CMA_eigenmethod'] = cma.Misc.eig
        
        inopts['verb_filenameprefix'] = inopts['verb_filenameprefix'] + 'Small';

#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
#TODO check if param exist
        suma_opts = suma_inopts;

        self._random_projection = self.grp(suma_opts['SUMA_d'], suma_opts['SUMA_D'], suma_opts['SUMA_spc']);
        
        super().__init__(self._random_projection @ x0, sigma0, inopts);

#TODO does not work; must be before calling super so check if exists in inopts an if not take default intead plus Small
#TODO EDIT better: inherit CMAOptions in way to set SUMA variables
        self.opts['verb_filenameprefix'] = self.opts['verb_filenameprefix'] + 'Small';

        self.inopts = {**self.inopts, **suma_inopts};
        self.opts = {**self.opts, **suma_opts};

        #test 
        self._first_identity = True;
        #self.gp.isidentity = False

    @staticmethod
    def grp(d, D, spc, eps=1e-6):
        res = np.zeros(shape=(d,D));
        i = 0;
        while i < d:
            rand_vector = new_basis_candidate = np.zeros(shape=(D), dtype=float);
            for j in range(D):
                if np.random.uniform() > spc:
                    rand_vector[j] = new_basis_candidate[j] = 2*np.random.randint(1)-1;
            j=0;
            for j in range(i):
               new_basis_candidate -= res[j] * np.dot(rand_vector,res[j]);
               #we don't need to divide by result.row(j).norm() cause all vectors are normalized
            if np.linalg.norm(new_basis_candidate) > eps:
                res[i] = new_basis_candidate / np.linalg.norm(new_basis_candidate);
                i += 1;
        return res;
    
    def reGrp(self):
        self._random_projection = self.grp(self.inopts['SUMA_d'], self.inopts['SUMA_D'], self.inopts['SUMA_spc']);
        
    def relativeProba(self, X):
        #renormalize in way to avoid centering effect
        Px = self._random_projection @ X - self.mean;
        if(np.linalg.norm(Px) > 0):
            Csx = self.Csi @ (Px / np.linalg.norm(Px));
            return -Csx.dot(Csx);
        else:
            return 0.0;
        #Csx = self.Csi @ (self._random_projection @ X - self.mean);
        ##the probability evoluates as -||Csx||^2 in a guassian distrib because of exp
        #return -Csx.dot(Csx);

    def bestCandidate(self, candidates):
        if len(candidates) < 1:
            raise cma._Error('candidates vector cannot be empty');
        else:
            best_seen = candidates[0];
            score_bs = self.relativeProba(best_seen);
            for i in range(1,len(candidates)):
                if self.relativeProba(candidates[i]) > score_bs:
                    best_seen = candidates[i];
                    score_bs = self.relativeProba(candidates[i]);
            return best_seen;

    def store_candidates(self, big_candidates):
        self._tmp_candidates = np.transpose(self._random_projection @ np.transpose(big_candidates));

    def begin_iteration(self):
#        self._flgtelldone = False;
#        self.evaluations_per_f_value = 1
#        self._tmp_candidates = np.transpose(self._random_projection @ np.transpose(big_candidates));
        #if all([self.D[i] == 1.0 for i in range(self.N)]) and not self._first_identity:
        #    raise cma._Error('identity again');
        #else:
        #    self._first_identity = False;
        #    self._savedC = self.C;

        self._tmp_candidates = None;
        self.Csi = np.dot(self.B, (self.B / self.D).T)  # square root of inverse of C according to cma.py update_exponential

    # ____________________________________________________________
    # ____________________________________________________________
    def ask_geno(self, number=None, xmean=None, sigma_fac=1):
        #TODO: document it
        """get new candidate solutions in genotyp, sampled from a
        multi-variate normal distribution.

        Arguments are
            `number`
                number of returned solutions, by default the
                population size `popsize` (AKA lambda).
            `xmean`
                distribution mean
            `sigma_fac`
                multiplier for internal sample width (standard
                deviation)

        `ask_geno` returns a list of N-dimensional candidate solutions
        in genotyp representation and is called by `ask`.

        Details: updates the sample distribution and might change
        the geno-pheno transformation during this update.

        :See: `ask`, `ask_and_eval`

        """
        # update distribution, might change self.mean
        if self.sp.CMA_on and (
                (self.opts['updatecovwait'] is None and
                 self.countiter >=
                     self.itereigenupdated + 1. / (self.sp.c1 + self.sp.cmu) / self.N / 10
                 ) or
                (self.opts['updatecovwait'] is not None and
                 self.countiter > self.itereigenupdated + self.opts['updatecovwait']
                 ) or
                (self.sp.neg.cmuexp * (self.countiter - self.itereigenupdated) > 0.5
                )  # TODO (minor): not sure whether this is "the right" criterion
            ):
            self.updateBD()

        if self._flgtelldone:  # could be done in tell()!?
            self._flgtelldone = False

        self.evaluations_per_f_value = 1

        return self._tmp_candidates;

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
        self._smallstrat = SumaSmall(x0, 1.0, copy.copy(inopts), suma_inopts);

        inopts['CMA_diagonal'] = True;

        super().__init__(x0, sigma0, inopts);

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
        
        #pop = np.matrix();
        pop = np.ndarray(shape=(number, self.N));
        self._smallstrat.begin_iteration();

        for i in range(number):
#TODO see EDIT above
#TODO check if param exist
            candidates = super().ask_geno(self.opts['SUMA_k']);
            pop[i] = self._smallstrat.bestCandidate(candidates);

        self._smallstrat.store_candidates(pop);

        return pop;

    def tell(self, solutions, function_values, check_points=None, copy=False):
#        Psolutions = np.transpose(self._smallstrat._random_projection @ np.transpose(solutions));
       # print(solutions)
       # for i in range(len(solutions)):
       #     print(Psolutions[i])
        #    print(self._smallstrat._random_projection.dot(solutions[i]))
        #self._smallstrat.tell(Psolutions, function_values, check_points, copy);
        self._smallstrat.tell(self._smallstrat.ask(), function_values, check_points, copy);
        super().tell(solutions, function_values, check_points, copy);
        #put in optimize
        self._smallstrat.logger.add(self._smallstrat);

