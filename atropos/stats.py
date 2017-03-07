# coding: utf-8
"""Collect statistics to use in the QC report.
"""

# TODO:
# * Refactor:
#   * Output of stats.py is a dict (serializable to/from JSON)
#   * Each stat is generated by a separate class. User can specify which stats
#     they want collected and/or plug in their own stat classes.
# * Enhancments:
#   * More generic way for filters and modifiers to report summary information

from collections import OrderedDict, defaultdict
import re
from atropos.adapters import BACK, FRONT, PREFIX, SUFFIX, ANYWHERE, LINKED
from atropos.modifiers import *
from atropos.filters import *
from atropos.util import qual2int

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
        super(MergingDict, self).__init__(*args, **kwargs)
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
    """Manages :class:`ReadStatCollector`s for pre- and post-trimming stats.
    
    Args:
        tile_key_regexp: Regular expression to parse read names and capture the
            read's 'tile' ID.
    """
    def __init__(self, mode, paired, **kwargs):
        self.mode = mode
        self.paired = paired
        self.pre = None
        self.post = None
        self.collector_args = kwargs
        
        if mode in ('pre', 'both'):
            self.pre = self._make_collectors()
        if mode in ('post', 'both'):
            self.post = {}
    
    def _make_collectors(self):
        return [
            ReadStatCollector(**self.collector_args)
            for i in range(2 if self.paired else 1)]
    
    def pre_trim(self, record):
        if self.pre is None:
            return
        if self.paired:
            self.pre[0].collect(record[0])
            self.pre[1].collect(record[1])
        else:
            self.pre[0].collect(record)
    
    def post_trim(self, dest, record):
        if self.post is None:
            return
        if dest not in self.post:
            self.post[dest] = self._make_collectors()
        post = self.post[dest]
        post[0].collect(record[0])
        if self.paired:
            post[1].collect(record[1])
    
    def finish(self):
        result = {}
        if self.pre is not None:
            result['pre'] = dict(
                ('read{}'.format(read), stats.finish())
                for read, stats in enumerate(self.pre, 1))
        if self.post is not None:
            result['post'] = {}
            for dest, collectors in self.post.items():
                result[post][dest] = dict(
                    ('read{}'.format(read), stats.finish())
                    for read, stats in enumerate(collectors, 1))
        return result

class ReadStatCollector(object):
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
        self.tile_key_regexp = tile_key_regexp
        self.qualities = qualities
        self.sequence_qualities = self.base_qualities = self.tile_base_qualities = None
        if qualities:
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
    
    def collect(self, record):
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
                    raise Exception("{} did not match {}".format(
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
    
    def finish(self):
        result = dict(
            count=self.count,
            length=self.sequence_lengths.sorted_items(),
            gc=self.sequence_gc.sorted_items(),
            bases=self.bases.flatten(datatype="n"))
        if self.sequence_qualities:
            result['qualities'] = self.sequence_qualities.sorted_items()
        if self.base_qualities:
            result['base_qualities'] = self.base_qualities.flatten(datatype="q")
        if self.track_tiles:
            result['tile_base_qualities'] = self.tile_base_qualities.flatten(datatype="q")
            result['tile_sequence_qualities'] = self.tile_sequence_qualities.flatten(shape="wide")
        return result

def collect_process_statistics(N, total_bp1, total_bp2, modifiers, filters, formatters):
    """
    Collect statistics from filters, and modifiers
    """
    stats = defaultdict(lambda: [0,0])
    
    written, written_bp = formatters.summary()
    assert written is not None
    
    stats.update(dict(
        N=N,
        total_bp1=total_bp1,
        total_bp2=total_bp2,
        total_bp=total_bp1 + total_bp2,
        written=written,
        written_bp=written_bp,
        total_written_bp=sum(written_bp),
    ))
    
    stats["too_short"] = None
    if TooShortReadFilter in filters:
        stats["too_short"] = filters[TooShortReadFilter].filtered
    
    stats["too_long"] = None
    if TooLongReadFilter in filters:
        stats["too_long"] = filters[TooLongReadFilter].filtered
    
    stats["too_many_n"] = None
    if NContentFilter in filters:
        stats["too_many_n"] = filters[NContentFilter].filtered
    
    # TODO: generalize this
    if modifiers.has_modifier(InsertAdapterCutter):
        insert_cutter = modifiers.get_modifiers(InsertAdapterCutter)[0]
        stats["with_adapters"] = insert_cutter.with_adapters
        stats["corrected"] = insert_cutter.corrected_pairs
        stats["corrected_bp"] = insert_cutter.corrected_bp
        stats["total_corrected_bp"] = sum(insert_cutter.corrected_bp)
    else:
        stats["with_adapters"] = [0, 0]
        if modifiers.has_modifier(AdapterCutter):
            adapter_cutters = modifiers.get_modifiers(AdapterCutter)[0]
            for read, modifier in enumerate(adapter_cutters):
                if modifier:
                    stats["with_adapters"][read] += modifier.with_adapters
            
    for modifier_class in modifiers.get_trimmer_classes():
        for modifier in modifiers.get_modifiers(modifier_class, 1):
            key = "{}_bp".format(type(modifier).__name__)
            stats[key][0] = modifier.trimmed_bases
                
        for modifier in modifiers.get_modifiers(modifier_class, 2):
            key = "{}_bp".format(type(modifier).__name__)
            stats[key][1] = modifier.trimmed_bases
        
        name = modifier_class.__name__
        stats[name] = sum(stats["{}_bp".format(name)])
    
    return dict(stats)

def collect_adapter_statistics(adapters):
    """
    TODO: push this method into Adapter?
    """
    stats = OrderedDict()
    for adapter in adapters:
        name = adapter.name
        total_front = sum(adapter.lengths_front.values())
        total_back = sum(adapter.lengths_back.values())
        
        stats[name] = MergingDict(
            name=name,
            total_front=total_front,
            total_back=total_back,
            total=total_front + total_back
        )
        
        where = adapter.where
        assert (where in (ANYWHERE, LINKED) or
            (where in (BACK, SUFFIX) and total_front == 0) or
            (where in (FRONT, PREFIX) and total_back == 0)
        )
        stats[name].set_with_merger("where", where, "check_equal")
        
        def handle_nested_dict(key, value):
            d = {}
            for k,v in value.items():
                d[k] = dict(v)
            stats[name].set_with_merger(key, d, "nested_dict")
        
        if where == LINKED:
            stats[name]["front_sequence"] = adapter.front_adapter.sequence
            stats[name]["back_sequence"] = adapter.back_adapter.sequence
            stats[name].set_with_merger("front_max_error_rate",
                adapter.front_adapter.max_error_rate, "check_equal")
            stats[name].set_with_merger("back_max_error_rate",
                adapter.back_adapter.max_error_rate, "check_equal")
            stats[name]["front_lengths_front"] = dict(adapter.front_adapter.lengths_front)
            stats[name]["front_lengths_back"] = dict(adapter.front_adapter.lengths_back)
            stats[name]["back_lengths_front"] = dict(adapter.back_adapter.lengths_front)
            stats[name]["back_lengths_back"] = dict(adapter.back_adapter.lengths_back)
            # have to clone these nested dicts and set them
            # up with a custom merge function
            handle_nested_dict("front_errors_front", adapter.front_adapter.errors_front)
            handle_nested_dict("front_errors_back", adapter.front_adapter.errors_back)
            handle_nested_dict("back_errors_front", adapter.back_adapter.errors_front)
            handle_nested_dict("back_errors_back", adapter.back_adapter.errors_back)
        else:
            stats[name]["sequence"] = adapter.sequence
            stats[name].set_with_merger("max_error_rate",
                adapter.max_error_rate, "check_equal")
            if where in (ANYWHERE, FRONT, PREFIX):
                stats[name]["lengths_front"] = dict(adapter.lengths_front)
                handle_nested_dict("errors_front", adapter.errors_front)
            if where in (ANYWHERE, BACK, SUFFIX):
                stats[name]["lengths_back"] = dict(adapter.lengths_back)
                handle_nested_dict("errors_back", adapter.errors_back)
            if where in (BACK, SUFFIX):
                stats[name]["adjacent_bases"] = dict(adapter.adjacent_bases)
    
    return stats

class Summary(object):
    def __init__(self, process_stats=MergingDict(),
                 adapter_stats=[OrderedDict(),OrderedDict()],
                 trimmer_classes=[]):
        self.process_stats = process_stats
        self.adapter_stats = adapter_stats
        self.trimmer_classes = trimmer_classes
    
    def add_process_stats(self, stats):
        self.process_stats.merge(stats)
    
    def add_read_stats(self, stats):
        raise NotImplementedError()
    
    def add_adapter_stats(self, stats):
        """
        stats = [ {name : adapter_stats}, {name : adapter_stats} ]
        """
        for i in (0,1):
            for name in stats[i].keys():
                if name in self.adapter_stats[i]:
                    self.adapter_stats[i][name].merge(stats[i][name])
                else:
                    self.adapter_stats[i][name] = stats[i][name]
    
    def finish(self):
        stats = self.process_stats.copy()
        
        stats["written_fraction"] = 0
        stats["too_short_fraction"] = 0
        stats["too_long_fraction"] = 0
        stats["too_many_n_fraction"] = 0
        stats["with_adapters_fraction"] = [0, 0]
        stats["total_written_bp_fraction"] = 0.0
        for modifier_class in self.trimmer_classes:
            name = modifier_class.__name__
            stats["{}_fraction".format(name)] = 0.0
            
        N = stats["N"]
        if N > 0:
            stats["written_fraction"] = stats["written"] / N if stats["written"] else 0
            stats["too_short_fraction"] = stats["too_short"] / N if stats["too_short"] else 0
            stats["too_long_fraction"] = stats["too_long"] / N if stats["too_long"] else 0
            stats["too_many_n_fraction"] = stats["too_many_n"] / N if stats["too_many_n"] else 0
            stats["with_adapters_fraction"] = [ (v / N) for v in stats["with_adapters"] ]
            if "corrected" in stats:
                stats["corrected_fraction"] = stats["corrected"] / N
                
        if stats["total_bp"] > 0:
            N = stats["total_bp"]
            stats["total_written_bp_fraction"] = (stats["total_written_bp"] / N) if stats["total_written_bp"] else 0
            if "corrected" in stats:
                stats["corrected_bp_fraction"] = [ (c / N) for c in stats["corrected_bp"] ]
                stats["total_corrected_bp_fraction"] = stats["total_corrected_bp"] / N
            for modifier_class in self.trimmer_classes:
                name = modifier_class.__name__
                if stats[name]:
                    stats["{}_fraction".format(name)] = (stats[name] / N)
        
        stats["adapters"] = [
            self.adapter_stats[0].values(),
            self.adapter_stats[1].values()
        ]
        
        return stats
