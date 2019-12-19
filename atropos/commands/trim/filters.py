"""
Classes for writing and filtering of processed reads.

A Filter is a callable that has the read as its only argument. If it is called,
it returns True if the read should be filtered (discarded), and False if not.

To be used, a filter needs to be wrapped in one of the Wrapper classes. To determine
what happens to a read, a list of Wrappers with different filters is created and
each Wrapper is called in turn until one returns True. The main program will
determine whether and where to write the read(s) based on whether it was rejected
by a filter, and which one.
"""
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from typing import Callable, Optional, Type

from atropos.io.sequence import Sequence
from atropos.utils import classproperty
from atropos.utils.collections import Summarizable


# Constants used when returning from a Filter’s __call__ method to improve
# readability (it is unintuitive that "return True" means "discard the read").
DISCARD = True
KEEP = False


class Filter(metaclass=ABCMeta):
    """
    Marker type for Filters.
    """

    @classproperty
    def name(cls) -> str:
        return cls.__name__

    @abstractmethod
    def __call__(self, read: Sequence) -> bool:
        pass


class MergedReadFilter(Filter):
    """
    Returns True if the read was merged.
    """

    def __call__(self, read: Sequence) -> bool:
        return read.merged


class TooShortReadFilter(Filter):
    """
    Returns True if the read sequence is shorter than `minimum_length`.
    """

    @classproperty
    def name(cls) -> str:
        return "too_short"

    def __init__(self, minimum_length: int):
        self.minimum_length = minimum_length

    def __call__(self, read: Sequence) -> bool:
        return len(read) < self.minimum_length


class TooLongReadFilter(Filter):
    """
    Returns True if the read sequence is longer than `maximum_length`.
    """

    @classproperty
    def name(cls) -> str:
        return "too_long"

    def __init__(self, maximum_length: int):
        self.maximum_length = maximum_length

    def __call__(self, read: Sequence) -> bool:
        return len(read) > self.maximum_length


class NContentFilter(Filter):
    """
    Discards a reads that has a number of 'N's over a given threshold. It handles
    both raw counts of Ns as well as proportions. Note, for raw counts, it is a
    greater than comparison, so a cutoff of '1' will keep reads with a single N in it.
    """

    @classproperty
    def name(cls) -> str:
        return "too_many_n"

    def __init__(self, count: int):
        """
        Args:
            count: If it is below 1.0, it will be considered a proportion, and above
            and equal to 1 will be considered as discarding reads with a number of
            N's greater than this cutoff.
        """
        if count < 0:
            raise ValueError(f"Count {count} is < 0")

        self.is_proportion = count < 1.0
        self.cutoff = count

    def __call__(self, read: Sequence) -> bool:
        """
        Returns True when the read should be discarded.
        """
        n_count = read.sequence.lower().count("n")
        if not self.is_proportion:
            return n_count > self.cutoff
        elif len(read) == 0:
            return KEEP
        else:
            return n_count / len(read) > self.cutoff


class UntrimmedFilter(Filter):
    """
    Returns True if read is untrimmed.
    """

    def __call__(self, read: Sequence) -> bool:
        return read.match is None


class TrimmedFilter(Filter):
    """
    Returns True if read is trimmed.
    """

    def __call__(self, read: Sequence) -> bool:
        return read.match is not None


class NoFilter(Filter):
    """
    Always returns False.
    """

    @classproperty
    def name(cls) -> str:
        return "NoFilter"

    def __call__(self, read: Sequence) -> bool:
        return KEEP


class FilterWrapper(Summarizable, metaclass=ABCMeta):
    """
    Base wrapper around a filter.
    """

    def __init__(self, f: Filter):
        self._filtered = 0
        self._wrapped = f

    def __call__(self, read1, read2=None):
        """
        Calls the filter function.

        Returns:
            DISCARD if the filter function returns True, else KEEP
        """
        if self._filter(read1, read2):
            self._filtered += 1
            return DISCARD

        return KEEP

    @abstractmethod
    def _filter(self, read1: Sequence, read2: Optional[Sequence] = None) -> bool:
        """
        Calls the filter function.
        """

    @property
    def name(self) -> str:
        """
        The filter name.
        """
        return self._wrapped.name

    def summarize(self) -> dict:
        """
        Returns a summary dict.
        """
        return dict(records_filtered=self._filtered)


class SingleWrapper(FilterWrapper):
    """This is for single-end reads and for paired-end reads, using the 'legacy'
    filtering mode (backwards compatibility). That is, if the first read matches
    the filtering criteria, the pair is discarded. The second read is not
    inspected.
    """

    def _filter(self, read1: Sequence, read2: Optional[Sequence] = None) -> bool:
        return self._wrapped(read1)


class PairedWrapper(FilterWrapper):
    """
    This is for paired-end reads, using the 'new-style' filtering where both reads
    are inspected. That is, the entire pair is discarded if at least 1 or 2 of the
    reads match the filtering criteria.
    """

    def __init__(self, f, min_affected: int = 1):
        """
        Args:
            min_affected: values 1 and 2 are allowed.
                1 means: the pair is discarded if any read matches
                2 means: the pair is discarded if both reads match
        """
        if min_affected not in (1, 2):
            raise ValueError("min_affected must be 1 or 2")

        super().__init__(f)
        self.min_affected = min_affected

    def _filter(self, read1: Sequence, read2: Optional[Sequence] = None) -> bool:
        failures = 0

        if self._wrapped(read1):
            failures += 1

        if (self.min_affected - failures == 1) and (
            read2 is None or self._wrapped(read2)
        ):
            failures += 1

        return failures >= self.min_affected


class FilterFactory:
    """
    Factory that creates filters and wraps them in the appropriate (single-end or
    paired-end) wrapper.
    """

    def __init__(self, paired: bool, min_affected: int):
        self.paired = paired
        self.min_affected = min_affected

    def __call__(self, filter_type: Callable[..., Filter], *args, **kwargs):
        """
        Creates a wrapped filter of the specified type.
        """
        fltr = filter_type(*args, **kwargs)
        if self.paired == "both":
            return PairedWrapper(fltr, self.min_affected)
        else:
            return SingleWrapper(fltr)


class Filters(Summarizable):
    """
    Manages multiple filters.

    Args:
        filter_factory: The :class:`FilterFactory` to use for creating new filters.
    """

    def __init__(self, filter_factory: FilterFactory):
        self.filters = OrderedDict()
        self.filter_factory = filter_factory

    def add_filter(self, filter_type: Type[Filter], *args, **kwargs):
        """
        Add a filter of the specified type.
        """
        self.filters[filter_type] = self.filter_factory(filter_type, *args, **kwargs)

    def filter(self, read1: Sequence, read2: Optional[Sequence] = None) -> Type[Filter]:
        """
        Executes filters in the order they were added until one returns True.

        Args:
            read1: The first read to filter.
            read2: The second read to filter.

        Returns:
            The type of the first filter that returned True, or :class:NoFilter
            if none of the filters returned True.
        """
        for filter_type, fltr in self.filters.items():
            if fltr(read1, read2):
                return filter_type
        else:
            return NoFilter

    def __contains__(self, filter_type: Type[Filter]) -> bool:
        return filter_type in self.filters

    def __getitem__(self, filter_type: Type[Filter]) -> Filter:
        return self.filters[filter_type]

    def summarize(self) -> dict:
        """
        Returns a summary dict.
        """
        return dict((f.name, f.summarize()) for f in self.filters.values())
