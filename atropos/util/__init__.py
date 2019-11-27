"""Widely useful utility methods.
"""
from collections import OrderedDict, Iterable, Sequence
from datetime import datetime
import errno
import functools
import logging
import math
from numbers import Number
import time
from atropos import AtroposError

# TODO: the nucleotide table should be implemented as an alphabet.

class NotInAlphabetError(Exception):
    def __init__(self, character):
        super().__init__()
        self.character = character

class Alphabet():
    def __init__(self, valid_characters, default_character):
        if not isinstance(valid_characters, set):
            valid_characters = set(valid_characters)
        if not default_character in valid_characters:
            valid_characters.add(default_character)
        self.valid_characters = valid_characters
        self.default_character = default_character

    def __contains__(self, character):
        return character in self.valid_characters

    def validate(self, character):
        """Raises NotInAlphabetError if the character is not in the alphabet.
        """
        if not character in self:
            raise NotInAlphabetError(character)

    def validate_string(self, string):
        """Raises NotInAlphabetError if any character in 'string' is not in the
        alphabet.
        """
        for character in string:
            self.validate(character)

    def resolve(self, character):
        """Returns 'character' if it's in the alphabet, otherwise the
        alphabet's default character.
        """
        if character in self.valid_characters:
            return character
        else:
            return self.default_character

    def resolve_string(self, string):
        """Returns a new string with any non-alphabet characters replaced
        with the default character.
        """
        return "".join(self.resolve(c) for c in string)

ALPHABETS = dict(
    dna=Alphabet('ACGT', 'N'),
    iso=None,
    colorspace=Alphabet('0123', None)
)

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

IUPAC_BASES = frozenset(('X',) + tuple(BASE_COMPLEMENTS.keys()))
"""Valid IUPAC bases, plus 'X'"""

GC_BASES = frozenset('CGRYSKMBDHVN')
"""IUPAC bases that include C or G."""

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
        key = (matches, size, match_prob)
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
                # use integer division in the case that the numbers are too
                # large for floating point division
                try:
                    div = nfac / self.factorial(i) / self.factorial(j)
                except OverflowError:
                    div = nfac // self.factorial(i) // self.factorial(j)
                prob += (mismatch_prob ** j) * (match_prob ** i) * div

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

class Mergeable(object):
    """Base class for objects that can merge themselves with another.
    """
    def merge(self, other):
        """Merges `other` with `self` and returns the merged value.
        """
        raise NotImplementedError()

class Summarizable(object):
    """Base class for objects that can summarize themselves.
    """
    def summarize(self):
        """Returns a summary dict.
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

class Timestamp(object):
    """Records datetime and clock time at object creation.
    """
    def __init__(self):
        self.dtime = datetime.now()
        self.process_time = time.process_time()

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
            cpu=max(minval, self.process_time - other.process_time))

class Timing(Summarizable):
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

class CountingDict(dict, Mergeable, Summarizable):
    """A dictionary that always returns 0 on get of a missing key.

    Args:
        sort_by: Whether summary is sorted by key (0) or value (1).
    """
    def __init__(self, keys=None, sort_by=0, summary_type='dict'):
        super().__init__()
        self.sort_by = sort_by
        self.summary_type = summary_type
        if keys:
            for key in keys:
                self.increment(key)

    def __getitem__(self, name):
        return self.get(name, 0)

    def increment(self, key, inc=1):
        """Increment the count of `key` by `inc`.
        """
        self[key] += inc

    def merge(self, other):
        if not isinstance(other, CountingDict):
            raise ValueError(
                "Cannot merge object of type {}".format(type(other)))
        for key, value in other.items():
            self[key] += value
        return self

    def get_sorted_items(self):
        """Returns an iterable of (key, value) sorted according to this
        CountingDict's `sort_by` param.
        """
        return sorted(self.items(), key=lambda item: item[self.sort_by])

    def summarize(self):
        """Returns an OrderedDict of sorted items.
        """
        summary_func = ordered_dict if self.summary_type == 'dict' else tuple
        return summary_func(self.get_sorted_items())

class Histogram(CountingDict):
    """Counting dict that returns a summary dict that contains summary stats.
    """
    def summarize(self):
        hist = super().summarize()
        return dict(
            hist=hist,
            summary=self.get_summary_stats())

    def get_summary_stats(self):
        """Returns dict with mean, median, and modes of histogram.
        """
        values = tuple(self.keys())
        counts = tuple(self.values())
        mu0 = weighted_mean(values, counts)
        return dict(
            mean=mu0,
            stdev=weighted_stdev(values, counts, mu0),
            median=weighted_median(values, counts),
            modes=weighted_modes(values, counts))

class NestedDict(dict, Mergeable, Summarizable):
    """A dict that initalizes :class:`CountingDict`s for missing keys.

    Args:
        shape: The flattened shape: 'long' or 'wide'.
    """
    def __init__(self, shape="wide"):
        super().__init__()
        self.shape = shape

    def __getitem__(self, name):
        if name not in self:
            self[name] = CountingDict()
        return self.get(name)

    def merge(self, other):
        if not isinstance(other, NestedDict):
            raise ValueError(
                "Cannot merge object of type {}".format(type(other)))
        for key, value in other.items():
            if key in self:
                self[key].merge(value)
            else:
                self[key] = value
        return self

    def summarize(self):
        """Returns a flattened version of the nested dict.

        Returns:
            When `shape=='long'`, a list of (key1, key2, value) tuples.
            When `shape=='wide'`, a dict of
                {columns:keys2, rows: {key1, values}}, where `keys2` is the set
                of keys in the child dicts.
        """
        keys1 = sorted(self.keys())
        if self.shape == "long":
            return tuple(
                (key1, key2, value)
                for key1 in keys1
                for key2, value in self[key1].items())
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
    """Create an OrderedDict from an iterable of (key, value) tuples.
    """
    ordict = OrderedDict()
    for key, value in iterable:
        ordict[key] = value
    return ordict

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

def mean(values):
    """Computes the mean of a sequence of numeric values.

    Args:
        values: Sequence of numeric values.

    Returns:
        The mean (floating point).
    """
    if len(values) == 0:
        raise ValueError("Cannot determine the mode of an empty sequence")
    return sum(values) / len(values)

def weighted_mean(values, counts):
    """Computes the mean of a sequence of numeric values weighted by counts.

    Args:
        values: Sequence of numeric values.

    Returns:
        The weighted mean (floating point).
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the mena of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")
    return sum(v * c for v, c in zip(values, counts)) / sum(counts)

def stdev(values, mu0=None):
    """Returns standard deviation of values having the specified mean.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the stdev of an empty sequence")
    if datalen == 1:
        return 0
    if mu0 is None:
        mu0 = mean(values)
    return math.sqrt(sum((val - mu0) ** 2 for val in values) / len(values))

def weighted_stdev(values, counts, mu0=None):
    """Returns standard deviation of values having the specified mean weighted
    by counts.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the stdev of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")
    if datalen == 1:
        return 0
    if mu0 is None:
        mu0 = weighted_mean(values, counts)
    return math.sqrt(
        sum(
            ((val - mu0) ** 2) * count
            for val, count in zip(values, counts)) /
        sum(counts))

def median(values):
    """Median function borrowed from python statistics module, and sped up by
    in-place sorting of the array.

    Args:
        data: Sequence of numeric values.

    Returns:
        The median (floating point).
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the median of an empty sequence")

    values.sort()

    idx = datalen // 2
    if datalen % 2 == 1:
        return values[idx]
    else:
        return (values[idx - 1] + values[idx]) / 2

def weighted_median(values, counts):
    """Compute the median of `values` weighted by `counts`.

    Args:
        values: Sequence of unique values.
        counts: Sequence of counts, where each count is the number of times the
            value at the corresponding position appears in the sample.

    Returns:
        The weighted median.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the median of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")
    counts_cumsum = functools.reduce(
        lambda c, x: c + [c[-1] + x], counts, [0])[1:]
    total = counts_cumsum[-1]
    if total == 0:
        return None
    mid1 = mid2 = (total // 2) + 1
    if total % 2 == 0:
        mid1 -= 1
    val1 = val2 = None
    for i, val in enumerate(counts_cumsum):
        if val1 is None and mid1 <= val:
            val1 = values[i]
        if mid2 <= val:
            val2 = values[i]
            break
    return float(val1 + val2) / 2

def modes(values):
    """Returns a sorted sequence of the modal (i.e. most frequent) values.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the mode of an empty sequence")
    elif datalen == 1:
        return values
    return _find_modes(CountingDict(values).items())

def weighted_modes(values, counts):
    """Returns a sorted sequence of the modal (i.e. most frequent) values
    weighted by counts.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the mode of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")
    if datalen == 1:
        return values
    return _find_modes(zip(values, counts))

def _find_modes(value_count_iter):
    sorted_counts = sorted(value_count_iter, key=lambda x: x[1], reverse=True)
    modal_values = [sorted_counts[0][0]]
    mode_count = sorted_counts[0][1]
    for value, count in sorted_counts[1:]:
        if count == mode_count:
            modal_values.append(value)
        else:
            break
    modal_values.sort()
    return modal_values

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
    except KeyboardInterrupt:
        logging.getLogger().error("Interrupted")
        retcode = 130
    except IOError as err:
        if err.errno == errno.EPIPE:
            retcode = 1
        else:
            raise
    except (AtroposError, EOFError):
        logging.getLogger().error("Atropos error", exc_info=True)
        retcode = 1
    except Exception: # pylint: disable=broad-except
        logging.getLogger().error("Unknown error", exc_info=True)
        retcode = 1
    return retcode
