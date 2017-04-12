# coding: utf-8
"""Collect statistics to use in the QC report.
"""
import re
from atropos.util import (
    CountingDict, NestedDict, Histogram, Mergeable, Summarizable, ordered_dict, 
    qual2int)

DEFAULT_TILE_KEY_REGEXP = r"^(?:[^\:]+\:){4}([^\:]+)"
"""Regexp for the default Illumina read name format."""

class PositionDicts(Mergeable, Summarizable):
    """A sequence of dicts, one for each position in a sequence.
    
    Args:
        is_qualities: Whether values are base qualities.
        quality_base: Base for quality values.
    """
    def __init__(self, is_qualities=False, quality_base=33):
        self.dicts = []
        self.is_qualities = is_qualities
        self.quality_base = quality_base
    
    def __getitem__(self, idx):
        if idx >= len(self.dicts):
            self.extend(idx+1)
        return self.dicts[idx]
    
    def extend(self, size):
        """Extend the number of bases to `size`.
        """
        diff = size - len(self.dicts)
        if diff > 0:
            for _ in range(diff):
                self.dicts.append(self.dict_class())
    
    def merge(self, other):
        if not isinstance(other, BaseCountingDicts):
            raise ValueError(
                "Cannot merge object of type {}".format(type(other)))
        other_len = len(other.dicts)
        min_len = min(len(self.dicts), other_len)
        for i in range(min_len):
            self.dicts[i].merge(other.dicts[i])
        if other_len > min_len:
            self.dicts.extend(other.dicts[min_len:other_len])
    
    def summarize(self):
        raise NotImplementedError()

class BaseCountingDicts(PositionDicts):
    """A PositionDicts in which the items are CountingDicts that count the
    number of items associated with each nucleotide base.
    """
    dict_class = CountingDict
    
    def summarize(self):
        """Flatten into a table with N rows (where N is the size of the
        sequence) and the columns are counts by nucleotide.
        
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
            acgt = ('A','C','G','T')
            n_val = ('N',)
            columns = keys = acgt + tuple(keys - set(acgt + n_val)) + n_val
        return dict(
            columns=columns,
            rows=ordered_dict(
                (idx, tuple(dict_item.get(key, 0) for key in keys))
                for idx, dict_item in enumerate(self.dicts, 1)))

class BaseNestedDicts(PositionDicts):
    """A PositionDicts in which items are NestedDicts.
    """
    dict_class = NestedDict
    
    def summarize(self):
        """Flatten into a table of N*K rows, where N is the sequence size and
        K is the union of keys in the nested dicts, and the columns are counts
        by nucleotide.
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
                (idx, ordered_dict(
                    (key1, tuple(dict_item[key1].get(key2, 0) for key2 in keys2))
                    for key1 in keys1))
                for idx, dict_item in enumerate(self.dicts, 1)))

class ReadStatistics(object):
    """Accumulates statistics on sequencing reads.
    
    Args:
        qualities: Whether to collect base quality statistics.
        tiles: Whether to collect tile-level statistics. If True, the default
            regular expression is used, otherwise must be a regular expression
            string or compiled re for extracting the tile ID from the read name.
            Only applies to Illumina sequences.
    """
    def __init__(self, qualities=None, quality_base=33, tiles=None):
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
        
        # whether to collect base quality stats
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
        
        # cache of computed values
        self._cache = {}
    
    def _init_qualities(self):
        # per-sequence mean qualities
        self.sequence_qualities = Histogram()
        # per-position quality composition
        self.base_qualities = BaseCountingDicts(
            is_qualities=True, quality_base=self.quality_base)
        if self.tile_key_regexp:
            self.tile_base_qualities = BaseNestedDicts(
                is_qualities=True, quality_base=self.quality_base)
            self.tile_sequence_qualities = NestedDict()
    
    # These are attributes that are computed on the fly. If called by name
    # (without leading '_'), __getattr__ uses the method to compute the value
    # if it is not already cached; on subsequent calls, the cached value is
    # returned.
    
    def _gc_pct(self):
        return (
            sum(base['C'] + base['G'] for base in self.bases) /
            self.total_bases)
    
    def _total_bases(self):
        return sum(
            length * count
            for base in self.bases
            for length, count in base.items())
    
    def __getattr__(self, name):
        if name not in self._cache:
            func_name = '_' + name
            if not hasattr(self, func_name):
                raise ValueError("No function named {}".format(func_name))
            func = getattr(self, func_name)
            self._cache[name] = func()
        return self._cache[name]
    
    @property
    def track_tiles(self):
        """Whether tile statistics are being accumulated.
        """
        return self.qualities and self.tile_key_regexp is not None
    
    def collect_record(self, record):
        """Collect stats on a single sequence record.
        """
        if self.qualities is None and record.qualities:
            self.qualities = True
            self._init_qualities()
        
        seq = record.sequence
        seqlen = len(seq)
        
        self.count += 1
        self.sequence_lengths[seqlen] += 1
        
        if seqlen > 0:
            gc_pct = round((seq.count('C') + seq.count('G')) * 100 / seqlen)
            self.sequence_gc[gc_pct] += 1
            
            if seqlen > self.max_read_len:
                self._extend_bases(seqlen)
            
            quals = tile = None
            
            if self.qualities:
                quals = record.qualities
                # mean read quality
                # NOTE: we use round here, as opposed to FastQC which uses
                # floor, resulting in slightly different quality profiles
                meanqual = round(
                    sum(ord(q) - self.quality_base for q in quals) / seqlen)
                self.sequence_qualities[meanqual] += 1
                # tile ID
                if self.track_tiles:
                    tile_match = self.tile_key_regexp.match(record.name)
                    if tile_match:
                        tile = tile_match.group(1)
                        self.tile_sequence_qualities[tile][meanqual] += 1
                    else:
                        raise ValueError("{} did not match {}".format(
                            self.tile_key_regexp, record.name))
        
            # per-base nucleotide and quality composition
            for i, (base, qual) in enumerate(zip(seq, quals)):
                self.add_base(i, base, qual, tile)
        
        # TODO: positional k-mer profiles
    
    def collect(self, read1, read2=None):
        """Collect statistics on a pair of reads.
        """
        raise NotImplementedError()
    
    def add_base(self, i, base, qual=None, tile=None):
        """Add per-base information.
        
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
        """Returns a summary dict.
        """
        summary = dict(
            counts=self.count,
            lengths=self.sequence_lengths.summarize(),
            gc=self.sequence_gc.summarize(),
            bases=self.bases)
        if self.sequence_qualities:
            summary['qualities'] = self.sequence_qualities
        if self.base_qualities:
            summary['base_qualities'] = self.base_qualities
        if self.track_tiles:
            summary['tile_base_qualities'] = self.tile_base_qualities
            summary['tile_sequence_qualities'] = self.tile_sequence_qualities
        return summary

class SingleEndReadStatistics(ReadStatistics):
    """ReadStatistics for single-end data.
    """
    def collect(self, read1, read2=None):
        self.collect_record(read1)
    
    def summarize(self):
        return dict(read1=super().summarize())

class PairedEndReadStatistics(object):
    """ReadStatistics for paired-end data.
    """
    def __init__(self, **kwargs):
        self.read1 = ReadStatistics(**kwargs)
        self.read2 = ReadStatistics(**kwargs)
    
    def collect(self, read1, read2):
        """Collect statistics on a pair of reads.
        """
        self.read1.collect_record(read1)
        self.read2.collect_record(read2)
    
    def summarize(self):
        """Returns a summary dict.
        """
        return dict(
            read1=self.read1.summarize(),
            read2=self.read2.summarize())
