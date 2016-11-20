import math

base_complements = {
    'A' : 'T',
    'C' : 'G',
    'R' : 'Y',
    'S' : 'S',
    'W' : 'W',
    'K' : 'M',
    'B' : 'V',
    'D' : 'H',
    'N' : 'N'
}
for k,v in list(base_complements.items()):
    base_complements[v] = k
    base_complements[k.lower()] = v.lower()
    base_complements[v.lower()] = k.lower()

MAGNITUDE = dict(
    G=1E9,
    M=1E6,
    K=1E3
)

log2 = math.log(2)

def complement(seq):
    return "".join(base_complements[base] for base in seq)

def reverse_complement(seq):
    return "".join(base_complements[base] for base in reversed(seq))

def sequence_complexity(seq):
    """A simple measure of sequence complexity"""
    seq = seq.upper()
    seqlen = float(len(seq))
    term = 0
    for base in ('A','C','G','T'):
        c = seq.count(base)
        if c > 0:
            d = c / seqlen
            term += d * math.log(d) / log2
    return -term

def qual2int(q, base=33):
    """Convert a quality charater to a phred-scale int.
    """
    return ord(q) - base

def quals2ints(quals, base=33):
    """Convert an iterable of quality characters to phred-scale ints.
    """
    return (ord(q) - base for q in quals)

def enumerate_range(collection, start, end):
    'Generates an indexed series:  (0,coll[0]), (1,coll[1]) ...'
    i = start
    it = iter(collection)
    while i < end:
        yield (i, next(it))
        i += 1

def mean(data):
    return sum(data) / len(data)

def median(data):
    """
    Median function borrowed from python statistics module, and sped up by
    in-place sorting of the array.
    """
    n = len(data)
    if n == 0:
        raise Exception("no median for empty data")
    
    data.sort()
    
    i = n // 2
    if n % 2 == 1:
        return data[i]
    else:
        return (data[i - 1] + data[i]) / 2

class RandomMatchProbability(object):
    """
    Class for computing random match probability for DNA sequences
    based on binomial expectation. Maintains a cache of factorials
    to speed computation.
    """
    def __init__(self, init_size=150):
        self.cache = {}
        self.factorials = [1] * init_size
        self.max_n = 1
        self.cur_array_size = init_size
    
    def __call__(self, matches, size, p1=0.25, p2=0.75):
        # First see if we have the result in the cache
        key = (matches, size)
        p = self.cache.get(key, None)
        if p:
            return p
        
        # When there are no mismatches, the probability is
        # just that of observing a specific sequence of the
        # given length by chance.
        if matches == size:
            p = p1 ** matches
        
        else:
            nfac = self.factorial(size)
            p = 0.0
            for i in range(matches, size+1):
                j = size - i
                p += (
                    (p2 ** j) *
                    (p1 ** i) *
                    nfac /
                    self.factorial(i) /
                    self.factorial(j)
                )

        self.cache[key] = p
        return p
    
    def factorial(self, n):
        if n > self.max_n:
            self._fill_upto(n)
        return self.factorials[n]

    def _fill_upto(self, n):
        if n >= self.cur_array_size:
            extension_size = n - self.cur_array_size + 1
            self.factorials += [1] * extension_size
        i = self.max_n
        next_i = i + 1
        while i < n:
            self.factorials[next_i] = next_i * self.factorials[i]
            i = next_i
            next_i += 1
        self.max_n = i
