import functools
import math
from typing import Iterable, Optional, Sequence, Tuple

from atropos.utils.collections import CountingDict


class RandomMatchProbability:
    """
    Class for computing random match probability for DNA sequences based on binomial
    expectation. Maintains a cache of factorials to speed computation.

    Args:
        init_size: Initial cache size.
    """

    def __init__(self, init_size: int = 150):
        self.cache = {}
        self.factorials = [1] * init_size
        self.max_n = 1
        self.cur_array_size = init_size

    def __call__(
        self,
        matches: int,
        size: int,
        match_prob: float = 0.25,
        mismatch_prob: float = 0.75,
    ) -> float:
        """
        Computes the random-match probability for a given sequence size and number of
        matches.

        Args:
            matches: Number of matches (numerator)
            size: Size of the sequence (denominator)
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
            for i in range(matches, size + 1):
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

    def factorial(self, num: int) -> int:
        """
        Returns `num`!.
        """
        if num > self.max_n:
            self._fill_upto(num)
        return self.factorials[num]

    def _fill_upto(self, num: int) -> None:
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


class Histogram(CountingDict):
    """
    CountingDict that returns a summary dict that contains summary stats.
    """

    def summarize(self) -> dict:
        hist = super().summarize()
        return dict(hist=hist, summary=self.get_summary_stats())

    def get_summary_stats(self) -> dict:
        """Returns dict with mean, median, and modes of histogram.
        """
        values = tuple(self.keys())
        counts = tuple(self.values())
        mu0 = weighted_mean(values, counts)
        return dict(
            mean=mu0,
            stdev=weighted_stdev(values, counts, mu0),
            median=weighted_median(values, counts),
            modes=weighted_modes(values, counts),
        )


def mean(values: Sequence[float]) -> float:
    """Computes the mean of a sequence of numeric values.

    Args:
        values: Sequence of numeric values.

    Returns:
        The mean (floating point).
    """
    if len(values) == 0:
        raise ValueError("Cannot determine the mode of an empty sequence")
    return sum(values) / len(values)


def weighted_mean(values: Sequence[float], counts: Sequence[int]) -> float:
    """
    Computes the mean of a sequence of numeric values weighted by counts.

    Args:
        values: Sequence of numeric values.
        counts: Sequence of counts.

    Returns:
        The weighted mean (floating point).
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the mena of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")
    return sum(v * c for v, c in zip(values, counts)) / sum(counts)


def stdev(values: Sequence[float], mu0: Optional[float] = None) -> float:
    """
    Computes the standard deviation of values having the specified mean.

    Args:
        values: The values from which to compute the standard deviation.
        mu0: The population mean. Computed from `values` if None.

    Returns:
        The standard deviation (floating point). Note: if `values` is of length 1,
        typically a ValueError would be raised, but this function returns 0.

    Raises:
        ValueError: If `values` is empty.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the stdev of an empty sequence")
    if datalen == 1:
        return 0
    if mu0 is None:
        mu0 = mean(values)
    return math.sqrt(sum((val - mu0) ** 2 for val in values) / (len(values) - 1))


def weighted_stdev(
    values: Sequence[float], counts: Sequence[int], mu0: Optional[float] = None
) -> float:
    """
    Computes the standard deviation of values having the specified mean weighted
    by counts.

    Args:
        values: The values from which to compute the standard deviation.
        counts: Sequence of counts.
        mu0: The population mean. Computed from `values` if None.

    Returns:
        The standard deviation (floating point). Note: if `values` is of length 1,
        typically a ValueError would be raised, but this function returns 0.

    Raises:
        ValueError: If `values` is empty or if `values` and `counts` are different
        lengths.
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
        sum(((val - mu0) ** 2) * count for val, count in zip(values, counts))
        / (sum(counts) - 1)
    )


def median(values: Sequence[float]) -> float:
    """
    Median function borrowed from python statistics module, and sped up by in-place
    sorting of the array.

    Args:
        values: Sequence of numeric values.

    Returns:
        The median (floating point).

    Raises:
        ValueError if `values` is empty.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the median of an empty sequence")

    values = sorted(values)
    idx = datalen // 2

    if datalen % 2 == 1:
        return values[idx]
    else:
        return (values[idx - 1] + values[idx]) / 2


def weighted_median(values: Sequence[float], counts: Sequence[int]) -> Optional[float]:
    """
    Computes the median of `values` weighted by `counts`.

    Args:
        values: Sequence of unique values.
        counts: Sequence of counts, where each count is the number of times the
            value at the corresponding position appears in the sample.

    Returns:
        The weighted median, or None if the total of `counts` is 0.

    Raises:
        ValueError if `values` is empty or a different length than `counts`.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the median of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")

    counts_cumsum = functools.reduce(lambda c, x: c + [c[-1] + x], counts, [0])[1:]
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


def modes(values: Sequence[float]) -> Sequence[float]:
    """
    Computes a sorted sequence of the modal (i.e. most frequent) values.

    Args:
        values: The values for which to find the mode(s).

    Returns:
        A sequence of the modal value(s).

    Raises:
        ValueError if `values` is empty.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the mode of an empty sequence")
    elif datalen == 1:
        return values
    return _find_modes(CountingDict(values).items())


def weighted_modes(values: Sequence[float], counts: Sequence[int]) -> Sequence[float]:
    """
    Computes a sorted sequence of the modal (i.e. most frequent) values weighted by
    counts.

    Args:
        values: Sequence of unique values.
        counts: Sequence of counts, where each count is the number of times the
            value at the corresponding position appears in the sample.

    Returns:
        A sequence of the modal value(s).

    Raises:
        ValueError if `values` is empty or a different length than `counts`.
    """
    datalen = len(values)
    if datalen == 0:
        raise ValueError("Cannot determine the mode of an empty sequence")
    if datalen != len(counts):
        raise ValueError("'values' and 'counts' must be the same length")
    if datalen == 1:
        return values
    return _find_modes(zip(values, counts))


def _find_modes(value_count_iter: Iterable[Tuple[float, int]]) -> Sequence[float]:
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
