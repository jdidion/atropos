from collections import OrderedDict
import logging
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

def check_equal_merger(dest, src):
    assert dest == src
    return dest

def nested_dict_merger(d_dest1, d_src1):
    assert isinstance(d_src1, dict)
    for k1, d_src2 in d_src1.items():
        if k1 in d_dest1:
            d_dest2 = d_dest1[k1]
            for k2, v_src in d_src2.items():
                if k2 in d_dest2:
                    d_dest2[k2] += v_src
                else:
                    d_dest2[k2] = v_src
        else:
            d_dest1[k1] = d_src2
    return d_dest1

MERGERS = dict(
    check_equal=check_equal_merger,
    nested_dict=nested_dict_merger
)

class MergingDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mergers = {}
    
    def set_with_merger(self, key, value, merger):
        self[key] = value
        self.mergers[key] = merger
    
    def merge(self, src):
        for k, v_src in src.items():
            if k in self and self[k] is not None:
                if v_src is None:
                    continue
                v_self = self[k]
                if k in self.mergers:
                    merger = MERGERS[self.mergers[k]]
                    self[k] = merger(v_self, v_src)
                # default behavior: lists have two integers, which are summed;
                # dicts have integer values, which are summed; strings must be
                # identical; otherwise must be numeric and are summed
                elif isinstance(v_self, dict):
                    assert isinstance(v_src, dict)
                    for kk,vv in v_src.items():
                        if kk in v_self:
                            v_self[kk] += vv
                        else:
                            v_self[kk] = vv
                elif isinstance(v_self, list):
                    assert isinstance(v_src, list)
                    self[k] = [v_src[0]+v_self[0], v_src[1]+v_self[1]]
                elif isinstance(v_self, str):
                    assert v_self == v_src
                else:
                    self[k] = v_src + v_self
            else:
                self[k] = v_src
                if isinstance(src, MergingDict) and k in src.mergers:
                    self.mergers[k] = src.mergers[k]
    
    def handle_nested_dict(self, key, value):
        d = {}
        for k,v in value.items():
            d[k] = dict(v)
        self.set_with_merger(key, d, "nested_dict")

class CountingDict(dict):
    def __getitem__(self, name):
        return self.get(name, 0)
    
    def sorted_items(self, by='key'):
        """Sort items"""
        return tuple(sorted(
            self.items(), key=lambda item: item[0 if by == 'key' else 1]))

class NestedDict(dict):
    def __getitem__(self, name):
        if name not in self:
            self[name] = CountingDict()
        return self.get(name)
    
    def flatten(self, shape="long"):
        keys1 = sorted(self.keys())
        if shape == "long":
            return [
                (key1, key2, value)
                for key1 in keys1
                for key2, value in self[key1].items()
            ]
        else:
            keys2 = set()
            for child in self.values():
                keys2.update(child.keys())
            keys2 = tuple(sorted(keys2))
            return (keys2, [
                (key1,) + tuple(self[key1].get(key2, 0) for key2 in keys2)
                for key1 in keys1
            ])

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

def truncate_string(s, n=100):
    """Shorten string s to at most n characters, appending "..." if necessary."""
    if s is None:
        return None
    if len(s) > n:
        s = s[:n-3] + '...'
    return s

def run_interruptible_with_result(func, *args, **kwargs):
    # Return code (0=normal, anything else is an error)
    rc = 0
    result = None
    try:
        result = func(*args, **kwargs)
    except KeyboardInterrupt as e:
        logging.getLogger().error("Interrupted")
        rc = 130
    except IOError as e:
        if e.errno == errno.EPIPE:
            rc = 1
        else:
            raise
    except (AtroposError, EOFError) as e:
        logging.getLogger().error("Atropos error", exc_info=True)
        rc = 1
    except Exception as e:
        logging.getLogger().error("Unknown error", exc_info=True)
        rc = 1
    return (rc, result)

def run_interruptible(func, *args, **kwargs):
    rc, result = run_interruptible_with_result(func, *args, **kwargs)
    return rc
