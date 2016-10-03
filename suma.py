import cma
import copy

class SumaSmall(cma.CMAEvolutionStrategy):
#TODO check wheter x0 is copy of obj or not
    def __init__(self, x0, sigma0, inopts={}):
#        inopts['CMA_diagonal'] = False
#TODO check if param exist
        self._random_projection = grp(inopts['SUMA_d'], len(x0));
        super(Suma, self).__init__(self, x0, sigma0, inopts);

    def grp(d, D):
        for i in range(d):
            new_basis_candidate = np.ndarray((D,));
            for j in range(D):
                np.


class Suma(cma.CMAEvolutionStrategy):
#TODO check wheter x0 is copy of obj or not
    def __init__(self, x0, sigma0, inopts={}):
        self._smallstrat_inopts = copy.copy(inopts);
        self._smallstrat = SumaSmall(x0, sigma0, self._smalldim_inopts);

        inopts['CMA_diagonal'] = True;

        super(Suma, self).__init__(self, x0, sigma0, inopts);

#    def ask(self, number=None, xmean=None, sigma_fac=1, gradf=None, args=()):
    def ask_geno(self, number=None, xmean=None, sigma_fac=1):
        if number is None or number < 1: 
            number = self.sp.popsize;
        
        pop = np.matrix();

        for i in range(number):
#TODO check if param exist
            candidates = super(Suma, self).ask_geno(self, self.opts['SUMA_k']);
            pop[i] = self._smallstrat.bestCandidate(candidates);

        return pop;

