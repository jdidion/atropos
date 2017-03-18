# coding: utf-8
"""Collect statistics to use in the QC report.
"""
import re
from atropos.util import CountingDict, NestedDict, qual2int

DEFAULT_TILE_KEY_REGEXP = "@(((?:[^\:]+)\:){5})"
"""Regexp for the default Illumina read name format."""

class BaseDicts(object):
    def __init__(self):
        self.dicts = []
    
    def __getitem__(self, idx):
        return self.dicts[idx]
    
    def extend(self, size):
        n = size - len(self.dicts)
        if n > 0:
            for i in range(n):
                self.dicts.append(self.dict_class())

class BaseCountingDicts(BaseDicts):
    dict_class = CountingDict
    
    def flatten(self, datatype=None):
        """
        Args:
            datatype: 'b' for bases, 'q' for qualities, or None.
        """
        keys = set()
        for d in self.dicts:
            keys.update(d.keys())
        if datatype == "n":
            acgt = ('A','C','G','T')
            N = ('N',) if 'N' in keys else ()
            keys = acgt + tuple(keys - set(acgt + N)) + N
        else:
            keys = tuple(sorted(keys))
        header = tuple(qual2int(k) for k in keys) if datatype == "q" else keys
        return (header, [
            (i,) + tuple(d.get(k, 0) for k in keys)
            for i, d in enumerate(self.dicts, 1)
        ])

class BaseNestedDicts(BaseDicts):
    dict_class = NestedDict
    
    def flatten(self, datatype=None):
        keys1 = set()
        keys2 = set()
        for d1 in self.dicts:
            keys1.update(d1.keys())
            for d2 in d1.values():
                keys2.update(d2.keys())
        keys1 = tuple(sorted(keys1))
        keys2 = tuple(sorted(keys2))
        header = tuple(qual2int(k) for k in keys2) if datatype == "q" else keys2
        return (header, [
            (i, k1,) + tuple(d[k1].get(k2, 0) for k2 in keys2)
            for k1 in keys1
            for i, d in enumerate(self.dicts, 1)
        ])

class ReadStatistics(object):
    def __init__(self, qualities=None, tile_key_regexp=None):
        # max read length
        self.max_read_len = 0
        # read count
        self.count = 0
        # read length distribution
        self.sequence_lengths = CountingDict()
        # per-sequence GC percentage
        self.sequence_gc = CountingDict()
        # per-position base composition
        self.bases = BaseCountingDicts()
        
        # whether to collect base quality stats
        self.qualities = qualities
        self.tile_key_regexp = None
        self.sequence_qualities = None
        self.base_qualities = None
        self.tile_base_qualities = None
        
        if qualities:
            if tile_key_regexp is True:
                tile_key_regexp = DEFAULT_TILE_KEY_REGEXP
            if isinstance(tile_key_regexp, str):
                tile_key_regexp = re.compile(regexp)
            self.tile_key_regexp = tile_key_regexp
            self._init_qualities()
        
        # cache of computed values
        self._cache = {}
    
    def _init_qualities(self):
        # per-sequence mean qualities
        self.sequence_qualities = CountingDict()
        # per-position quality composition
        self.base_qualities = BaseCountingDicts()
        if self.tile_key_regexp:
            self.tile_base_qualities = BaseNestedDicts()
            self.tile_sequence_qualities = NestedDict()
    
    # These are attributes that are computed on the fly. If called by name
    # (without leading '_'), __getattr__ uses the method to compute the value
    # if it is not already cached; on subsequent calls, the cached value is
    # returned.
    
    def _gc_pct(self):
        return sum(base['C'] + base['G'] for base in self.bases) / self.total_bases
    
    def _total_bases(self):
        return sum(length * count for length, count in self.bases.items())
    
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
        return self.qualities and self.tile_key_regexp is not None
    
    def collect_record(self, record):
        if self.qualities is None and record.qualities:
            self.qualities = True
            self._init_qualities()
        
        seq = record.sequence
        seqlen = len(seq)
        gc = round((seq.count('C') + seq.count('G')) * 100 / seqlen)
        
        self.count += 1
        self.sequence_lengths[seqlen] += 1
        self.sequence_gc[gc] += 1
        
        if seqlen > self.max_read_len:
            self._extend_bases(seqlen)
        
        qual = tile = None
        
        if self.qualities:
            quals = record.qualities
            # mean read quality
            meanqual = round(sum(ord(q) for q in quals) / seqlen)
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
    
    def add_base(self, i, base, qual=None, tile=None):
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
        summary = dict(
            count=self.count,
            length=self.sequence_lengths.sorted_items(),
            gc=self.sequence_gc.sorted_items(),
            bases=self.bases.flatten(datatype="n"))
        if self.sequence_qualities:
            summary['qualities'] = self.sequence_qualities.sorted_items()
        if self.base_qualities:
            summary['base_qualities'] = self.base_qualities.flatten(datatype="q")
        if self.track_tiles:
            summary['tile_base_qualities'] = self.tile_base_qualities.flatten(datatype="q")
            summary['tile_sequence_qualities'] = self.tile_sequence_qualities.flatten(shape="wide")
        return summary

class SingleEndReadStatistics(ReadStatistics):
    def collect(self, read1, read2=None):
        self.collect_record(read1)
    
    def summarize(self):
        return dict(read1=super().summarize())

class PairedEndReadStatistics(object):
    def __init__(self, **kwargs):
        self.read1 = ReadStatistics(**kwargs)
        self.read2 = ReadStatistics(**kwargs)
    
    def collect(self, read1, read2):
        self.read1.collect_record(read1)
        self.read2.collect_record(read2)
    
    def summarize(self):
        return dict(
            read1=self.read1.summarize(),
            read2=self.read2.summarize())
