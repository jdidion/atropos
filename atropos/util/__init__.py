"""Widely useful utility methods.
"""
from collections import OrderedDict, Iterable
from datetime import datetime
import errno
import logging
import math
from numbers import Number
import time
from atropos import AtroposError

def build_iso_nucleotide_table():
    """Generate a dict mapping ISO nucleotide characters to their complements,
    in both upper and lower case.
    """
    nuc = {
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
    for base, comp in tuple(nuc.items()):
        nuc[comp] = base
        nuc[base.lower()] = comp.lower()
        nuc[comp.lower()] = base.lower()
    return nuc

BASE_COMPLEMENTS = build_iso_nucleotide_table()

MAGNITUDE = dict(
    G=1E9,
    M=1E6,
    K=1E3
)

LOG2 = math.log(2)

class RandomMatchProbability(object):
    """Class for computing random match probability for DNA sequences based on
    binomial expectation. Maintains a cache of factorials to speed computation.
    
    Args:
        init_size: Initial cache size.
    """
    def __init__(self, init_size=150):
        self.cache = {}
        self.factorials = [1] * init_size
        self.max_n = 1
        self.cur_array_size = init_size
    
    def __call__(self, matches, size, match_prob=0.25, mismatch_prob=0.75):
        """Computes the random-match probability for a given sequence size and
        number of matches.
        
        Args:
            match_prob: Probability of two random bases matching.
            mismatch_prob: Probability of two random bases not matcing.
        
        Returns:
            The probability.
        """
        # First see if we have the result in the cache
        key = (matches, size)
        prob = self.cache.get(key, None)
        if prob:
            return prob
        
        # When there are no mismatches, the probability is
        # just that of observing a specific sequence of the
        # given length by chance.
        if matches == size:
            prob = match_prob ** matches
        
        else:
            nfac = self.factorial(size)
            prob = 0.0
            for i in range(matches, size+1):
                j = size - i
                prob += (
                    (mismatch_prob ** j) *
                    (match_prob ** i) *
                    nfac /
                    self.factorial(i) /
                    self.factorial(j)
                )
        
        self.cache[key] = prob
        return prob
    
    def factorial(self, num):
        """Returns `num`!.
        """
        if num > self.max_n:
            self._fill_upto(num)
        return self.factorials[num]

    def _fill_upto(self, num):
        if num >= self.cur_array_size:
            extension_size = num - self.cur_array_size + 1
            self.factorials += [1] * extension_size
        idx = self.max_n
        next_i = idx + 1
        while idx < num:
            self.factorials[next_i] = next_i * self.factorials[idx]
            idx = next_i
            next_i += 1
        self.max_n = idx

class Timestamp(object):
    """Records datetime and clock time at object creation.
    """
    def __init__(self):
        self.dtime = datetime.now()
        self.clock = time.clock()
    
    def timestamp(self):
        """Returns the unix timestamp.
        """
        return self.dtime.timestamp()
    
    def isoformat(self):
        """Returns the datetime in ISO format.
        """
        return self.dtime.isoformat()
    
    def __sub__(self, other, minval=0.01):
        """Subtract another timestamp from this one.
        
        Args:
            other: The other timestamp.
            minval: The minimum difference.
        
        Returns:
            A dict of {wallclock=<datetime_diff>, cpu=<clock_diff>}.
        """
        return dict(
            wallclock=max(minval, self.timestamp() - other.timestamp()),
            cpu=max(minval, self.clock - other.clock))

class Timing(object):
    """Context manager that maintains timing information using
    :class:`Timestamp`s. Maintains a start time on __enter__, and can be updated
    with the current time by the user. Does a final update on __exit__.
    """
    def __init__(self):
        self.start_time = None
        self.cur_time = None
    
    def __enter__(self):
        self.start_time = Timestamp()
        return self
    
    def __exit__(self, exception_type, exception_value, traceback):
        self.update()
    
    def update(self):
        """Set :attr:`self.cur_time` to the current time.
        """
        self.cur_time = Timestamp()
    
    def summarize(self):
        """Returns a summary dict
        {start=<start_time>, wallclock=<datetime_diff>, cpu=<clock_diff>}.
        """
        if not self.cur_time:
            self.update()
        assert self.start_time is not None
        summary = dict(start=self.start_time.isoformat())
        summary.update(self.cur_time - self.start_time)
        return summary

class CountingDict(dict):
    """A dictionary that always returns 0 on get of a missing key.
    """
    def __getitem__(self, name):
        return self.get(name, 0)
    
    def sorted_items(self, sort_by=0):
        """Returns a tuple of sorted items.
        
        Args:
            sort_by: Whether to sort by key (0) or value (1).
        
        Returns:
            A tuple of sorted items.
        """
        return tuple(sorted(self.items(), key=lambda item: item[sort_by]))

class NestedDict(dict):
    """A dict that initalizes :class:`CountingDict`s for missing keys.
    """
    def __getitem__(self, name):
        if name not in self:
            self[name] = CountingDict()
        return self.get(name)
    
    def summarize(self, shape="long"):
        """Returns a flattened version of the nested dict.
        
        Args:
            shape: The flattened shape: 'long' or 'wide'.
        
        Returns:
            When `shape=='long'`, a list of (key1, key2, value) tuples.
            When `shape=='wide'`, a dict of
                {columns:keys2, rows: {key1, values}}, where `keys2` is the set
                of keys in the child dicts.
        """
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
            return dict(
                columns=keys2,
                rows=ordered_dict(
                    (key1, tuple(self[key1].get(key2, 0) for key2 in keys2))
                    for key1 in keys1))

class Mergeable(object):
    """Marker class for an object that can merge itself with another.
    """
    def merge(self, other):
        """Merges `other` with `self` and returns the merged value.
        """
        raise NotImplementedError()

class Const(Mergeable):
    """A :class:`Mergeable` that is a constant value. Merging simply checks
    that two values are identical.
    
    Args:
        value: The value to treat as a constant.
    """
    def __init__(self, value):
        self.value = value
    
    def merge(self, other):
        """Check that `self==other`.
        
        Raises:
            ValueError
        """
        if self != other:
            raise ValueError("{} != {}".format(self, other))
        return self
    
    def __eq__(self, other):
        """Returns True if `self.value==other` (or `other.value` if `other` is
        a :class:`Const`).
        """
        if isinstance(other, Const):
            other = other.value
        return self.value == other
    
    def __repr__(self):
        return str(self.value)

class MergingDict(OrderedDict, Mergeable):
    """An :class:`collections.OrderedDict` that implements :class:`Mergeable`.
    """
    def merge(self, other):
        """Merge `other` with `self` using :method:`merge_dicts`.
        """
        merge_dicts(self, other)
        return self

def merge_dicts(dest, src):
    """Merge corresponding items in `src` into `dest`. Values in `src` missing
    in `dest` are simply added to `dest`. Values that appear in both `src` and
    `dest` are merged using `merge_values`.
    
    Args:
        dest: The dict to merge into.
        src: The dict to merge from.
    
    Raises:
        ValueError if a value is not one of the accepted types.
    """
    for key, v_src in src.items():
        if dest.get(key, None) is None:
            dest[key] = v_src
        elif v_src is not None:
            dest[key] = merge_values(dest[key], v_src)

def merge_values(v_dest, v_src):
    """Merge two values based on their types, as follows:
    
    - Mergeable: merging is done by the dest object's merge function.
    - dict: merge_dicts is called recursively.
    - Number: values are summed.
    - Iterable (non-string): First src and dest values are converted to tuples;
    they must be the same length. Then, corresponding values are handled as
    above. The value is returned as a list.
    - Otherwise: Treated as a Const (i.e. must be identical).
    
    Args:
        v_dest: The dest value.
        v_src: The src value.
    
    Returns:
        The merged value.
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
        if len(i_dest) == 0:
            v_dest = i_src
        elif len(i_src) > 0:
            v_dest = [merge_values(d, s) for d, s in zip(i_dest, i_src)]
    else:
        assert v_dest == v_src
    return v_dest

def ordered_dict(iterable):
    d = OrderedDict()
    for key, value in iterable:
        d[key] = value
    return d

def complement(seq):
    """Returns the complement of nucleotide sequence `seq`.
    """
    return "".join(BASE_COMPLEMENTS[base] for base in seq)

def reverse_complement(seq):
    """Returns the reverse complement of nucleotide sequence `seq`.
    """
    return "".join(BASE_COMPLEMENTS[base] for base in reversed(seq))

def sequence_complexity(seq):
    """Computes a simple measure of sequence complexity.
    
    Args:
        seq: The sequence to measure.
    
    Returns:
        Complexity, as a value [0,2], where 0 = a homopolymer and
        2 = completely random.
    """
    seq = seq.upper()
    seqlen = float(len(seq))
    term = 0
    for base in ('A','C','G','T'):
        count = seq.count(base)
        if count > 0:
            frac = count / seqlen
            term += frac * math.log(frac) / LOG2
    return -term

def qual2int(qual, base=33):
    """Convert a quality charater to a phred-scale int.
    
    Args:
        q: The quality value.
        base: The offset of the first quality value (Old Illumina = 64,
            new Illumina = 33).
    
    Returns:
        The integer quality.
    """
    return ord(qual) - base

def quals2ints(quals, base=33):
    """Convert an iterable of quality characters to phred-scale ints.
    
    Args:
        quals: The qualities.
        base: The offset of the first quality value (Old Illumina = 64,
            new Illumina = 33).
    
    Returns:
        A tuple of integer qualities.
    """
    return (ord(q) - base for q in quals)

def qual2prob(qchar):
    """Converts a quality char to a probability.
    """
    return 10 ** (-qual2int(qchar) / 10)

def enumerate_range(collection, start, end):
    """Generates an indexed series:  (0,coll[0]), (1,coll[1]) ...
    
    Only up to (start-end+1) items will be yielded.
    
    Args:
        collection: The collection to enumerate.
        start: The starting index.
        end: The ending index.
    
    Yields:
        (index, item) tuples.
    """
    idx = start
    itr = iter(collection)
    while idx < end:
        yield (idx, next(itr))
        idx += 1

def mean(data):
    """Computes the mean of a sequence of numeric values.
    
    Args:
        data: Sequence of numeric values.
    
    Returns:
        The mean (floating point).
    """
    return sum(data) / len(data)

def median(data):
    """Median function borrowed from python statistics module, and sped up by
    in-place sorting of the array.
    
    Args:
        data: Sequence of numeric values.
    
    Returns:
        The median (floating point).
    """
    datalen = len(data)
    if datalen == 0:
        raise ValueError("no median for empty data")
    
    data.sort()
    
    idx = datalen // 2
    if datalen % 2 == 1:
        return data[idx]
    else:
        return (data[idx - 1] + data[idx]) / 2

def truncate_string(string, max_len=100):
    """Shorten string s to at most n characters, appending "..." if necessary.
    """
    if string is None:
        return None
    if len(string) > max_len:
        string = string[:max_len-3] + '...'
    return string

def run_interruptible(func, *args, **kwargs):
    """Run a function, gracefully handling keyboard interrupts.
    
    Args:
        func: The function to execute.
        args, kwargs: Positional and keyword arguments to pass to `func`.
    
    Returns:
        A (unix-style) return code (0=normal, anything else is an error).
    
    Raises:
        
    """
    retcode = 0
    try:
        func(*args, **kwargs)
    except KeyboardInterrupt as err:
        logging.getLogger().error("Interrupted")
        retcode = 130
    except IOError as err:
        if err.errno == errno.EPIPE:
            retcode = 1
        else:
            raise
    except (AtroposError, EOFError) as err:
        logging.getLogger().error("Atropos error", exc_info=True)
        retcode = 1
    except Exception as err: # pylint: disable=broad-except
        logging.getLogger().error("Unknown error", exc_info=True)
        retcode = 1
    return retcode
