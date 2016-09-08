import math

complement = {
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
for k,v in list(complement.items()):
    complement[v] = k
    complement[k.lower()] = v.lower()
    complement[v.lower()] = k.lower()

def reverse_complement(seq):
    return "".join(complement[base] for base in reversed(seq))

log2 = math.log(2)
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

MAGNITUDE = dict(
    G=(1E9, "000000000"),
    M=(1E6, "000000"),
    K=(1E3, "000")
)

def magnitude_formatter(magnitude):
    suffix = ""
    if magnitude is None:
        div = 1.0
    else:
        div = float(MAGNITUDE[magnitude.upper()][0])
        suffix = magnitude
    return lambda val: "{:.1f} {}".format(val / div, suffix)

def int_or_str(x):
    if x is None or isinstance(x, int):
        return x
    elif isinstance(x, str):
        x = x.upper()
        for a, mag in MAGNITUDE.items():
            x = x.replace(a, mag[1])
        return int(x)
    else:
        raise Exception("Unsupported type {}".format(x))

def enumerate_range(collection, start, end):
        'Generates an indexed series:  (0,coll[0]), (1,coll[1]) ...'
        i = start
        it = iter(collection)
        while i < end:
            yield (i, next(it))
            i += 1

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
        self.factorials = [1] * init_size
        self.max_n = 1
        self.cur_array_size = init_size
    
    def __call__(self, matches, size):
        # When there are no mismatches, the probability is
        # just that of observing a specific sequence of the
        # given length by chance.
        if matches == size:
            return 0.25 ** matches
        
        nfac = self.factorial(size)
        p = 0.0
        i = matches

        while i <= size:
            j = size - i
            p += (
                (0.75 ** j) *
                (0.25 ** i) *
                nfac /
                self.factorial(i) /
                self.factorial(j)
            )
            i += 1

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
