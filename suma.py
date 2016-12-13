import numpy as np
import cma

class Suma(cma.CMAEvolutionStrategy):
    def __init__(self, x0, sigma0, inopts={}):
        self.suma_opts = self.extract_suma_inopts(inopts)

        inopts['CMA_diagonal'] = True

        super(Suma, self).__init__(x0, sigma0, inopts)

        self._f = None

    def extract_suma_inopts(self, inopts):
        suma_keys = [k for k in inopts if k.startswith('SUMA')]
        suma_inopts = {}
        for key in suma_keys:
            suma_inopts[key] = inopts.pop(key, None)
        return suma_inopts

    def ask_geno(self, number=None, xmean=None, sigma_fac=1):
        if number is None or number < 1: 
            number = self.sp.popsize
        pop = np.ndarray(shape=(number, self.N))

        for i in range(number):
            candidates = super(Suma, self).ask_geno(self.suma_opts['SUMA_k'])

            if len(candidates) < 1:
                raise cma._Error('candidates vector cannot be empty')
            else:
                best_seen = candidates[0]
                score_bs = self._f(best_seen)
                for j in range(1, len(candidates)):
                    cs = self._f(candidates[j])
                    if cs < score_bs:
                        best_seen = candidates[j]
                        score_bs = cs

            pop[i] = best_seen

        return pop

    def tell(self, solutions, function_values, check_points=None, copy=False):
        s = np.mean(self.C)
        self.C /= s
        super(Suma, self).tell(solutions, function_values, check_points, copy)

    def optimize(self, objective_fct, iterations=None, min_iterations=1,
                 args=(), verb_disp=None, logger=None, call_back=None):
        self._f = objective_fct
        return super(Suma, self).optimize(objective_fct, iterations, min_iterations, args, verb_disp, logger, call_back)

