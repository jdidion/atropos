from collections import OrderedDict, Iterable
from datetime import datetime
import logging
import math
from numbers import Number
import time
from atropos import AtroposError

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

class Timestamp(object):
    def __init__(self):
        self.dt = datetime.now()
        self.clock = time.clock()
    
    def timestamp(self):
        return self.dt.timestamp()
    
    def isoformat(self):
        return self.dt.isoformat()
    
    def __sub__(self, other):
        return dict(
            wallclock=self.timestamp() - other.timestamp(),
            cpu=self.clock - other.clock)

class Timing(object):
    def __init__(self):
        self.start_time = None
        self.cur_time = None
    
    def __enter__(self):
        self.start_time = Timestamp()
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        self.update()
    
    def update(self):
        self.cur_time = Timestamp()
    
    def summarize(self):
        if not self.cur_time:
            self.update()
        assert self.start_time is not None
        summary = dict(start=self.start_time.isoformat())
        summary.update(self.cur_time - self.start_time)
        return summary

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

class Mergeable(object):
    def merge(self, other):
        raise NotImplementedError()

class Const(Mergeable):
    def __init__(self, value):
        self.value = value
    
    def merge(self, other):
        if isinstance(other, Const):
            other = other.value
        if type(other) != type(self.value):
            raise ValueError("{} cannot be merged with {}".format(self, other))
        if self.value != other:
            raise ValueError("{} != {}".format(self, other))
        return self.value
    
    def __repr__(self):
        return self.value

class MergingDict(OrderedDict):
    def merge(self, src):
        merge_dicts(self, src)

def merge_dicts(dest, src):
    """Merge corresponding items in `src` into `dest`. Values in `src` missing
    in `dest` are simply added to `dest`. Values that appear in both `src` and
    `dest` are merged using `merge_values`.
    
    Raises:
        ValueError if a value is not one of the accepted types.
    """
    for k, v_src in src.items():
        if dest.get(k, None) is None:
            dest[k] = v_src
        elif v_src is not None:
            dest[k] = merge_values(dest[k], v_src)

def merge_values(v_dest, v_src):
    """Merge two values based on their types, as follows:
    
    - Mergeable: merging is done by the dest object's merge function.
    - dict: merge_dicts is called recursively.
    - str: Treated as a Const (i.e. must be identical).
    - Number: values are summed.
    - Iterable: First src and dest values are converted to tuples; they must be
    the same length. Then, corresponding values are handled as above. The value
    is returned as a list.
    """
    if isinstance(v_dest, Mergeable):
        v_dest = v_dest.merge(v_src)
    elif isinstance(v_dest, dict):
        assert isinstance(v_src, dict)
        merge_dicts(v_dest, v_src)
    elif isinstance(v_dest, str):
        assert v_dest == v_src
    elif isinstance(v_dest, Number):
        v_dest += v_src
    elif isinstance(v_dest, Iterable):
        i_dest = tuple(v_dest)
        i_src = tuple(v_src)
        assert len(i_dest) == len(i_src)
        v_dest = [merge_values(d, s) for d, s in zip(i_dest, i_src)]
    else:
        raise ValueError("Unmergable value {}".format(v_dest))
    return v_dest

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
        raise ValueError("no median for empty data")
    
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

def run_interruptible(func, *args, **kwargs):
    # Return code (0=normal, anything else is an error)
    rc = 0
    try:
        func(*args, **kwargs)
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
    return rc
