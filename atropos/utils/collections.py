from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from datetime import datetime
from numbers import Number
import time
from typing import (
    Any,
    Generic,
    Iterable,
    Optional,
    Tuple,
    Type,
    TypeVar,
    Union,
)

from atropos.errors import AtroposError


class Mergeable(metaclass=ABCMeta):
    """
    Base class for objects that can merge themselves with another.
    """

    @abstractmethod
    def merge(self, other):
        """
        Merges `other` with `self` and returns the merged value.
        """
        pass


class Summarizable(metaclass=ABCMeta):
    """
    Base class for objects that can summarize themselves.
    """

    @abstractmethod
    def summarize(self) -> dict:
        """
        Returns a summary dict.
        """


ConstType = TypeVar("ConstType")


class Const(Generic[ConstType], Mergeable):
    """
    A :class:`Mergeable` that is a constant value. Merging simply checks that two
    values are identical.
    """

    def __init__(self, value: ConstType):
        """
        Args:
            value: The value to treat as a constant.
        """
        self.value = value

    def merge(self, other: "Const[ConstType]") -> "Const[ConstType]":
        """
        Checks that `self==other`.

        Raises:
            ValueError
        """
        if self != other:
            raise ValueError(f"{self} != {other}")

        return self

    def __eq__(self, other: Union[ConstType, "Const[ConstType]"]) -> bool:
        """
        Returns True if `self.value==other` (or `other.value` if `other` is a
        :class:`Const`).
        """
        if isinstance(other, Const):
            other = other.value
        return self.value == other

    def __repr__(self) -> str:
        return str(self.value)


class Timestamp:
    """
    Records datetime and clock time at object creation.
    """

    def __init__(self):
        self.dtime = datetime.now()
        self.process_time = time.process_time()

    @property
    def timestamp(self) -> float:
        """
        The unix timestamp.
        """
        return self.dtime.timestamp()

    @property
    def isoformat(self) -> str:
        """
        The datetime in ISO format.
        """
        return self.dtime.isoformat()

    def __sub__(self, other: "Timestamp", minval: float = 0.01) -> dict:
        """
        Subtracts another timestamp from this one.

        Args:
            other: The other timestamp.
            minval: The minimum difference.

        Returns:
            A dict of {wallclock=<datetime_diff>, cpu=<clock_diff>}.
        """
        return dict(
            wallclock=max(minval, self.timestamp - other.timestamp),
            cpu=max(minval, self.process_time - other.process_time),
        )


class Timing(Summarizable):
    """
    Context manager that maintains timing information using :class:`Timestamp`s.
    Maintains a start time on __enter__, and can be updated with the current time by
    the user. Does a final update on __exit__.
    """

    def __init__(self):
        self.start_time = None
        self.cur_time = None

    def __enter__(self):
        self.start_time = Timestamp()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.update()

    def update(self) -> None:
        """
        Sets :attr:`self.cur_time` to the current time.
        """
        self.cur_time = Timestamp()

    def summarize(self) -> dict:
        """
        Returns a summary dict
        {start=<start_time>, wallclock=<datetime_diff>, cpu=<clock_diff>}.
        """
        if not self.cur_time:
            self.update()

        if self.start_time is None:
            raise AtroposError(
                "Timing instance must be started before it can be summarized"
            )

        summary = dict(start=self.start_time.isoformat)
        summary.update(self.cur_time - self.start_time)

        return summary


SummaryType = TypeVar("SummaryType")


class BaseMergeableDict(dict, Mergeable, Summarizable, metaclass=ABCMeta):
    pass


class CountingDict(BaseMergeableDict):
    """
    A dictionary that always returns 0 on get of a missing key.
    """

    def __init__(
        self,
        keys: Optional[Iterable] = None,
        sort_by: int = 0,
        summary_type: Type[SummaryType] = dict,
    ):
        """
        Args:
            keys: The initial keys to count.
            sort_by: Whether summary is sorted by key (0) or value (1).
            summary_type:
        """
        super().__init__()
        self.sort_by = sort_by
        self.summary_type = summary_type
        if keys:
            for key in keys:
                self.increment(key)

    def __getitem__(self, key) -> int:
        return self.get(key, 0)

    def increment(self, key, inc: int = 1) -> None:
        """
        Increments the count of `key` by `inc`.
        """
        self[key] += inc

    def merge(self, other: "CountingDict") -> "CountingDict":
        if not isinstance(other, CountingDict):
            raise ValueError(f"Cannot merge object of type {type(other)}")
        for key, value in other.items():
            self[key] += value
        return self

    def get_sorted_items(self) -> Iterable[Tuple[Any, int]]:
        """
        Returns an iterable of (key, value) sorted according to this CountingDict's
        `sort_by` param.
        """
        return sorted(self.items(), key=lambda item: item[self.sort_by])

    def summarize(self) -> SummaryType:
        """
        Returns an OrderedDict of sorted items.
        """
        summary_func = ordered_dict if self.summary_type == dict else tuple
        return summary_func(self.get_sorted_items())


class NestedDict(BaseMergeableDict):
    """
    A dict that initalizes :class:`CountingDict`s for missing keys.
    """

    def __init__(self, summary_type: Type[SummaryType] = dict):
        """
        Args:
            summary_type: The flattened shape: list or dict.
        """
        super().__init__()
        self.summary_type = summary_type

    def __getitem__(self, key) -> CountingDict:
        if key not in self:
            self[key] = CountingDict()

        return self.get(key)

    def merge(self, other: "NestedDict") -> "NestedDict":
        if not isinstance(other, NestedDict):
            raise ValueError(f"Cannot merge object of type {type(other)}")
        for key, value in other.items():
            if key in self:
                self[key].merge(value)
            else:
                self[key] = value
        return self

    def summarize(self) -> SummaryType:
        """
        Returns a flattened version of the nested dict.

        Returns:
            When `self.summary_type==list`, a list of (key1, key2, value) tuples.
            When `self.summary_type==dict`, a dict of
                {columns:keys2, rows: {key1, values}}, where `keys2` is the set
                of keys in the child dicts.
        """
        keys1 = sorted(self.keys())
        if self.summary_type == list:
            return tuple(
                (key1, key2, value)
                for key1 in keys1
                for key2, value in self[key1].items()
            )
        else:
            keys2 = set()
            for child in self.values():
                keys2.update(child.keys())
            keys2 = tuple(sorted(keys2))
            return dict(
                columns=keys2,
                rows=ordered_dict(
                    (key1, tuple(self[key1].get(key2, 0) for key2 in keys2))
                    for key1 in keys1
                ),
            )


class MergingDict(OrderedDict, Mergeable):
    """
    An :class:`collections.OrderedDict` that implements :class:`Mergeable`.
    """

    def merge(self, other):
        """
        Merges `other` with `self` using :method:`merge_dicts`.
        """
        merge_dicts(self, other)
        return self


def merge_dicts(dest: dict, src: dict):
    """
    Merges corresponding items in `src` into `dest`. Values in `src` missing in `dest`
    are simply added to `dest`. Values that appear in both `src` and `dest` are
    merged using `merge_values`.

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
            dest[key] = merge_values(key, dest[key], v_src)


def merge_values(key, v_dest, v_src):
    """
    Merges two values based on their types, as follows:

    - Mergeable: merging is done by the dest object's merge function.
    - dict: merge_dicts is called recursively.
    - Number: values are summed.
    - Iterable (non-string): First src and dest values are converted to tuples;
    they must be the same length. Then, corresponding values are handled as
    above. The value is returned as a list.
    - Otherwise: Treated as a Const (i.e. must be identical).

    Args:
        key: The key being merged.
        v_dest: The dest value.
        v_src: The src value.

    Returns:
        The merged value.
    """

    def assert_of_type(_type=None):
        if _type is None:
            _type = v_dest.__class__
        if not all(isinstance(o, _type) for o in (v_src, v_dest)):
            raise TypeError(
                f"When merging {key}, source and/or dest were not of type {_type}"
            )
        return _type

    def assert_equal():
        _type = assert_of_type()
        if v_dest != v_src:
            raise ValueError(
                f"When merging, objects of type {_type} must be equal; "
                f"{v_src} != {v_dest}"
            )

    if isinstance(v_dest, Mergeable):
        v_dest = v_dest.merge(v_src)
    elif isinstance(v_dest, dict):
        assert_of_type(dict)
        merge_dicts(v_dest, v_src)
    elif isinstance(v_dest, str) or isinstance(v_dest, bool):
        assert_equal()
    elif isinstance(v_dest, Number):
        assert_of_type(Number)
        v_dest += v_src
    elif isinstance(v_dest, Iterable):
        assert_of_type(Iterable)
        i_dest = tuple(v_dest)
        i_src = tuple(v_src)
        if len(i_dest) == 0:
            v_dest = i_src
        elif len(i_src) > 0:
            v_dest = [merge_values(key, d, s) for d, s in zip(i_dest, i_src)]
    else:
        assert_equal()
    return v_dest


def ordered_dict(iterable: Iterable) -> OrderedDict:
    """
    Creates an OrderedDict from an iterable of (key, value) tuples.

    Todo: this should no longer be necessary given that dicts are implicitly ordered
     in python 3.6+
    """
    ordict = OrderedDict()
    for key, value in iterable:
        ordict[key] = value
    return ordict
