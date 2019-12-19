from abc import ABCMeta, abstractmethod
from enum import Enum
from functools import lru_cache
import re
from typing import List, Optional, Pattern, Union

from atropos.io.sequence import Sequence
from atropos.utils import classproperty
from atropos.utils.collections import (
    BaseMergeableDict,
    CountingDict,
    Mergeable,
    NestedDict,
    Summarizable,
    ordered_dict,
)
from atropos.utils.statistics import Histogram
from atropos.utils.ngs import qual2int


DEFAULT_TILE_KEY_REGEXP = r"^(?:[^\:]+\:){4}([^\:]+)"
"""Regexp for the default Illumina read name format."""


class MetricsMode(Enum):
    PRE = "pre"
    POST = "post"


class PositionDicts(Mergeable, Summarizable, metaclass=ABCMeta):
    """
    A sequence of dicts, one for each position in a sequence.
    """

    @classproperty
    def _create_dict(cls) -> BaseMergeableDict:
        pass

    def __init__(self, is_qualities: bool = False, quality_base: int = 33):
        """
        Args:
            is_qualities: Whether values are base qualities.
            quality_base: Base for quality values.
        """
        self.dicts: List[BaseMergeableDict] = []
        self.is_qualities = is_qualities
        self.quality_base = quality_base

    def __getitem__(self, idx: int) -> dict:
        if idx >= len(self.dicts):
            self.extend(idx + 1)
        return self.dicts[idx]

    def extend(self, size: int) -> None:
        """
        Extends the number of bases to `size`.
        """
        diff = size - len(self.dicts)
        if diff > 0:
            for _ in range(diff):
                self.dicts.append(self._create_dict())

    def merge(self, other: "PositionDicts") -> None:
        if not isinstance(other, PositionDicts):
            raise ValueError(f"Cannot merge object of type {type(other)}")

        other_len = len(other.dicts)
        min_len = min(len(self.dicts), other_len)

        for i in range(min_len):
            self.dicts[i].merge(other.dicts[i])

        if other_len > min_len:
            self.dicts.extend(other.dicts[min_len:other_len])

    @abstractmethod
    def summarize(self) -> dict:
        pass


class BaseCountingDicts(PositionDicts):
    """
    A PositionDicts in which the items are CountingDicts that count the number of
    items associated with each nucleotide base.
    """

    @classproperty
    def _create_dict(cls) -> BaseMergeableDict:
        return CountingDict()

    def summarize(self) -> dict:
        """
        Flattens into a table with N rows (where N is the size of the sequence) and the
        columns are counts by nucleotide.

        Returns:
            A tuple of (columns, [rows]), where each row is
            (position, (base_counts...))
        """
        keys = set()

        for dict_item in self.dicts:
            keys.update(dict_item.keys())

        if self.is_qualities:
            keys = tuple(sorted(keys))
            columns = tuple(qual2int(k, self.quality_base) for k in keys)
        else:
            acgt = tuple("ACGT")
            n_val = tuple("N")
            columns = keys = acgt + tuple(keys - set(acgt + n_val)) + n_val

        return dict(
            columns=columns,
            rows=ordered_dict(
                (idx, tuple(dict_item.get(key, 0) for key in keys))
                for idx, dict_item in enumerate(self.dicts, 1)
            ),
        )


class BaseNestedDicts(PositionDicts):
    """
    A PositionDicts in which items are NestedDicts.
    """

    @classproperty
    def _dict_class(cls) -> BaseMergeableDict:
        return NestedDict()

    def summarize(self) -> dict:
        """
        Flattens into a table of N*K rows, where N is the sequence size and K is the
        union of keys in the nested dicts, and the columns are counts by nucleotide.
        """
        keys1 = set()
        keys2 = set()

        for dict1 in self.dicts:
            keys1.update(dict1.keys())
            for dict2 in dict1.values():
                keys2.update(dict2.keys())

        keys1 = tuple(sorted(keys1))
        keys2 = tuple(sorted(keys2))

        if self.is_qualities:
            columns = tuple(qual2int(k, self.quality_base) for k in keys2)
        else:
            columns = keys2

        return dict(
            columns=columns,
            columns2=keys1,
            rows=ordered_dict(
                (
                    idx,
                    ordered_dict(
                        (key1, tuple(dict_item[key1].get(key2, 0) for key2 in keys2))
                        for key1 in keys1
                    ),
                )
                for idx, dict_item in enumerate(self.dicts, 1)
            ),
        )


class ReadMetricsCollector(Summarizable):
    """
    Accumulates metrics on sequencing reads.

    Todo: positional k-mer profiles
    """

    def __init__(
        self,
        qualities=None,
        quality_base: int = 33,
        tiles: Union[bool, str, Pattern] = None,
    ):
        """
        Args:
            qualities: Whether to collect base quality statistics.
            tiles: Whether to collect tile-level statistics. If True, the default
                regular expression is used, otherwise must be a regular expression
                string or compiled re for extracting the tile ID from the read name.
                Only applies to Illumina sequences.
        """
        # max read length
        self.max_read_len = 0
        # read count
        self.count = 0
        # read length distribution
        self.sequence_lengths = Histogram()
        # per-sequence GC percentage
        self.sequence_gc = Histogram()
        # per-position base composition
        self.bases = BaseCountingDicts()
        # whether to collect base quality metrics
        self.qualities = qualities
        self.quality_base = quality_base
        self.tile_key_regexp = None
        self.sequence_qualities = None
        self.base_qualities = None
        self.tile_base_qualities = None

        if qualities:
            tile_key_regexp = DEFAULT_TILE_KEY_REGEXP if tiles is True else tiles
            if isinstance(tile_key_regexp, str):
                tile_key_regexp = re.compile(tile_key_regexp)
            self.tile_key_regexp = tile_key_regexp
            self._init_qualities()

    def _init_qualities(self) -> None:
        # per-sequence mean qualities
        self.sequence_qualities = Histogram()
        # per-position quality composition
        self.base_qualities = BaseCountingDicts(
            is_qualities=True, quality_base=self.quality_base
        )
        if self.tile_key_regexp:
            self.tile_base_qualities = BaseNestedDicts(
                is_qualities=True, quality_base=self.quality_base
            )
            self.tile_sequence_qualities = NestedDict()

    @property
    @lru_cache(maxsize=None)
    def gc_pct(self) -> float:
        return sum(base["C"] + base["G"] for base in self.bases) / self.total_bases

    @property
    @lru_cache(maxsize=None)
    def total_bases(self) -> int:
        return sum(
            length * count for base in self.bases for length, count in base.items()
        )

    @property
    def track_tiles(self) -> bool:
        """
        Whether tile statistics are being accumulated.
        """
        return self.qualities and self.tile_key_regexp is not None

    def collect_record(self, record: Sequence):
        """
        Collects metrics on a single sequence record.
        """
        if self.qualities is None and record.qualities:
            self.qualities = True
            self._init_qualities()

        seq = record.sequence
        seqlen = len(seq)

        self.count += 1
        self.sequence_lengths[seqlen] += 1

        if seqlen > 0:
            gc_pct = round((seq.count("C") + seq.count("G")) * 100 / seqlen)
            self.sequence_gc[gc_pct] += 1

            if seqlen > self.max_read_len:
                self._extend_bases(seqlen)

            tile = None
            if self.track_tiles:
                tile_match = self.tile_key_regexp.match(record.name)
                if tile_match:
                    tile = tile_match.group(1)
                else:
                    raise ValueError(
                        f"{self.tile_key_regexp} did not match {record.name}"
                    )

            if self.qualities:
                quals = record.qualities

                # mean read quality
                # NOTE: we use round here, as opposed to FastQC which uses
                # floor, resulting in slightly different quality profiles
                meanqual = round(
                    sum(ord(q) - self.quality_base for q in quals) / seqlen
                )
                self.sequence_qualities[meanqual] += 1

                if self.track_tiles:
                    self.tile_sequence_qualities[tile][meanqual] += 1

                for i, (base, qual) in enumerate(zip(seq, quals)):
                    self.add_base(i, base, qual, tile)
            else:
                for i, base in enumerate(seq):
                    self.add_base(i, base, tile=tile)

    def add_base(self, i, base, qual=None, tile=None):
        """
        Adds per-base information.

        Args:
            i: Position
            base: Nucleotide
            qual: Quality
            tile: Tile ID
        """
        self.bases[i][base] += 1
        if qual:
            self.base_qualities[i][qual] += 1
            if tile:
                self.tile_base_qualities[i][tile][qual] += 1

    def _extend_bases(self, new_size):
        self.bases.extend(new_size)
        if self.qualities:
            self.base_qualities.extend(new_size)
            if self.track_tiles:
                self.tile_base_qualities.extend(new_size)

    def summarize(self):
        """
        Returns a summary dict.
        """
        summary = dict(
            counts=self.count,
            lengths=self.sequence_lengths.summarize(),
            gc=self.sequence_gc.summarize(),
            bases=self.bases,
        )
        if self.sequence_qualities:
            summary["qualities"] = self.sequence_qualities
        if self.base_qualities:
            summary["base_qualities"] = self.base_qualities
        if self.track_tiles:
            summary["tile_base_qualities"] = self.tile_base_qualities
            summary["tile_sequence_qualities"] = self.tile_sequence_qualities
        return summary


class ReadMetrics(Summarizable, metaclass=ABCMeta):
    @abstractmethod
    def collect(self, read1: Sequence, read2: Optional[Sequence] = None):
        """
        Collects statistics on a pair of reads.
        """


class SingleEndReadMetrics(ReadMetrics):
    """
    ReadStatistics for single-end data.
    """

    def __init__(self, **kwargs):
        self.read1 = ReadMetricsCollector(**kwargs)

    def collect(self, read1: Sequence, read2: Optional[Sequence] = None):
        self.read1.collect_record(read1)

    def summarize(self):
        return dict(read1=super().summarize())


class PairedEndReadMetrics(ReadMetrics):
    """
    ReadStatistics for paired-end data.
    """

    def __init__(self, **kwargs):
        self.read1 = ReadMetricsCollector(**kwargs)
        self.read2 = ReadMetricsCollector(**kwargs)

    def collect(self, read1: Sequence, read2: Optional[Sequence] = None):
        """Collect statistics on a pair of reads.
        """
        self.read1.collect_record(read1)
        self.read2.collect_record(read2)

    def summarize(self):
        """Returns a summary dict.
        """
        return dict(read1=self.read1.summarize(), read2=self.read2.summarize())
