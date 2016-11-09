import cma,suma
from multiprocessing import Process

def launch_optim(flavour, opts, suma_inopts, inopts, number=10):
    proc = []

    verb_filenameprefix_base = inopts['verb_filenameprefix']

    for i in range(number):
        es = None
        inopts['verb_filenameprefix'] = verb_filenameprefix_base + '_' + str(i)

        if flavour == 'A':
            inopts['CMA_diagonal'] = True
            es = cma.CMAEvolutionStrategy([opts['x0']] * (opts['dim'] * suma_inopts['SUMA_k']), opts['sigma0'], inopts)
        elif flavour == 'B':
            es = suma.Suma([opts['x0']] * opts['dim'], opts['sigma0'], {**inopts, **suma_inopts})
        elif flavour == 'C':
            es = suma.SumaC([opts['x0']] * opts['dim'], opts['sigma0'], {**inopts, **suma_inopts})
        elif flavour == 'D':
            es = suma.SumaD([opts['x0']] * opts['dim'], opts['sigma0'], {**inopts, **suma_inopts})

        #args must be iterable
        p = Process(target=es.optimize, args=[opts['func']])
        p.start()
        proc.append(p)

    for p in proc:
        p.join()

launch_optim(flavour='A', opts={'x0': 1.0, 'sigma0': 0.5, 'dim': 50, 'func': cma.fcts.ellirot}, inopts={'verb_filenameprefix': 's0:0.5_d:200_func:ellirot'}, suma_inopts={'SUMA_k': 40, 'SUMA_d': 40, 'SUMA_spc': 0.6}, number=2)

