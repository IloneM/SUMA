class RunIterator():
    def __init__(self, bounds, start=None):#bounds
        #if start is not None:
        #   if not len(start) == len(bounds):
        #TODO       raise + verif if in bounds
        self._bounds = bounds
        self.current = [0] * len(bounds) if start is None else list(start)
#        self.ended = len(bounds) <= 0 or min(bounds) <= 0 or (start is not None and len(start) is not len(bounds))
        self.ended = len(bounds) <= 0 or min(bounds) <= 0 or \
                                         (start is not None and \
                                         (len(start) is not len(bounds) or \
                                          [i for i in range(len(bounds)) if start[i] >= bounds[i]]))
    def __iter__(self):
        return self
    def __next__(self):
        return self.next()
    def next(self): # Python 3: def __next__(self)
        if self.ended:
            raise StopIteration
        self.ended = True
        res = tuple(self.current)
        for i in range(len(self._bounds)-1, -1, -1):
#        for i in range(len(self._bounds)): dont use this because lexico order is not respected in that case
            if self.current[i] < self._bounds[i]-1:
                self.current[i] += 1
                self.ended = False
                break
            else:
                self.current[i] = 0
        return res

def runGen(bounds, start = None):
    #if start is not None:
    #   if not len(start) == len(bounds):
    #TODO       raise + verif if in bounds
    current = [0] * len(bounds) if start is None else list(start)
    ended = len(bounds) <= 0 or min(bounds) <= 0 or \
                                (start is not None and \
                                (len(start) is not len(bounds) or \
                                 [i for i in range(len(bounds)) if start[i] >= bounds[i]]))

    while not ended:
        ended = True
        res = tuple(current)
        for i in range(len(bounds)-1, -1, -1):
#        for i in range(len(bounds)): dont use this because lexico order is not respected in that case
            if current[i] < bounds[i]-1:
                current[i] += 1
                ended = False
                break
            else:
                current[i] = 0
        yield res
