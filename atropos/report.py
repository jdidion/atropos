# coding: utf-8
"""
Routines for printing a report.
"""
import sys
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
import textwrap
from .adapters import BACK, FRONT, PREFIX, SUFFIX, ANYWHERE, LINKED
from .modifiers import *
from .filters import *

# TODO: Would be good to have a more generic way for filters
# and modifiers to report summary information

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
        total_written_bp=sum(written_bp)
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
    insert_cutter = modifiers.get_modifiers(InsertAdapterCutter)
    if len(insert_cutter) > 0:
        stats["with_adapters"] = insert_cutter[0].with_adapters
    else:
        stats["with_adapters"] = [0, 0]
        for read in (1, 2):
            for modifier in modifiers.get_modifiers(AdapterCutter, read):
                stats["with_adapters"][read-1] += modifier.with_adapters
            
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

def summarize_adapters(mods):
    summary = [{}, {}]
    for a1, a2 in mods:
        if a1 is not None:
            if len(summary[0]) == 0:
                summary[0] = collect_adapter_statistics(a1.adapters)
            else:
                raise Exception("Did not expect multiple read1 adapter modfiers")
        if a2 is not None:
            if len(summary[1]) == 0:
                summary[1] = collect_adapter_statistics(a2.adapters)
            else:
                raise Exception("Did not expect multiple read2 adapter modfiers")
    return summary

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
        
        if stats["total_bp"] > 0:
            stats["total_written_bp_fraction"] = (stats["total_written_bp"] / stats["total_bp"]) if stats["total_written_bp"] else 0
            for modifier_class in self.trimmer_classes:
                name = modifier_class.__name__
                if stats[name]:
                    stats["{}_fraction".format(name)] = (stats[name] / stats["total_bp"])
        
        stats["adapters"] = [
            self.adapter_stats[0].values(),
            self.adapter_stats[1].values()
        ]
        
        return stats

def print_report(paired, options, wallclock_time, cpu_time, stats, trimmer_classes):
    outfile = options.report_file
    close = False
    
    if outfile is not None:
        outfile = open(outfile, "w")
        close = True
    elif not options.quiet:
        outfile = sys.stderr if options.output is None else sys.stdout
    else:
        return
    
    stats["paired"] = paired
    stats["wallclock_time"] = max(wallclock_time, 0.01)
    stats["cpu_time"] = max(cpu_time, 0.01)
    stats['pairs_or_reads'] = "Pairs" if paired else "Reads"
    
    try:
        # TODO: rewrite generate_report so that this call
        # receives a string result and writes it to outfile
        # outfile.write(generate_report(paired, time, stats))
        generate_report(stats, trimmer_classes, outfile)

    finally:
        if close and outfile is not None:
            outfile.close()
    
        # def get_values(self):
        #	values = var(self)
        #	return values

def generate_report(stats, trimmer_classes, outfile):
    def _print(*args, **kwargs): print(*args, file=outfile, **kwargs)
    
    ADAPTER_TYPES = {
        BACK: "regular 3'",
        FRONT: "regular 5'",
        PREFIX: "anchored 5'",
        SUFFIX: "anchored 3'",
        ANYWHERE: "variable 5'/3'",
        LINKED: "linked",
    }
    
    def print_error_ranges(adapter_length, error_rate):
        _print("No. of allowed errors:")
        prev = 0
        for errors in range(1, int(error_rate * adapter_length) + 1):
            r = int(errors / error_rate)
            _print("{0}-{1} bp: {2};".format(prev, r - 1, errors - 1), end=' ')
            prev = r
        if prev == adapter_length:
            _print("{0} bp: {1}".format(adapter_length, int(error_rate * adapter_length)))
        else:
            _print("{0}-{1} bp: {2}".format(prev, adapter_length, int(error_rate * adapter_length)))
        _print()

    def print_histogram(d, adapter_length, n, error_rate, errors):
        """
        Print a histogram. Also, print the no. of reads expected to be
        trimmed by chance (assuming a uniform distribution of nucleotides in the reads).
        d -- a dictionary mapping lengths of trimmed sequences to their respective frequency
        adapter_length -- adapter length
        n -- total no. of reads.
        """
        h = []
        for length in sorted(d):
            # when length surpasses adapter_length, the
            # probability does not increase anymore
            estimated = n * 0.25 ** min(length, adapter_length)
            h.append( (length, d[length], estimated) )
        
        def errs_to_str(l, e):
            if e in errors[l]:
                return str(errors[l][e])
            return "0"
        
        _print("length", "count", "expect", "max.err", "error counts", sep="\t")
        for length, count, estimate in h:
            max_errors = max(errors[length].keys())
            errs = ' '.join(errs_to_str(length, e) for e in range(max_errors+1))
            _print(length, count, "{0:.1F}".format(estimate), int(error_rate*min(length, adapter_length)), errs, sep="\t")
        _print()

    def print_adjacent_bases(bases, sequence):
        """
        Print a summary of the bases preceding removed adapter sequences.
        Print a warning if one of the bases is overrepresented and there are
        at least 20 preceding bases available.

        Return whether a warning was printed.
        """
        total = sum(bases.values())
        if total == 0:
            return False
        _print('Bases preceding removed adapters:')
        warnbase = None
        for base in ['A', 'C', 'G', 'T', '']:
            b = base if base != '' else 'none/other'
            fraction = 1.0 * bases[base] / total
            _print('  {0}: {1:.1%}'.format(b, fraction))
            if fraction > 0.8 and base != '':
                warnbase = b
        if total >= 20 and warnbase is not None:
            _print('WARNING:')
            _print('	The adapter is preceded by "{0}" extremely often.'.format(warnbase))
            _print('	The provided adapter sequence may be incomplete.')
            _print('	To fix the problem, add "{0}" to the beginning of the adapter sequence.'.format(warnbase))
            _print()
            return True
        _print()
        return False
    
    #----------------
    
    """Print report to standard output."""
    if stats["N"] == 0:
        _print("No reads processed! Either your input file is empty or you used the wrong -f/--format parameter.")
        return
    _print("Wallclock time: {0:.2F} s ({1:.0F} us/read; {2:.2F} M reads/minute).".format(
        stats["wallclock_time"], 1E6 * stats["wallclock_time"] / stats["N"], stats["N"] / stats["wallclock_time"] * 60 / 1E6))
    _print("CPU time (main process): {0:.2F} s".format(stats["cpu_time"]))
    
    report = "\n=== Summary ===\n\n"
    if stats["paired"]:
        report += textwrap.dedent("""\
        Total read pairs processed:		 {N:13,d}
          Read 1 with adapter:			 {with_adapters[0]:13,d} ({with_adapters_fraction[0]:.1%})
          Read 2 with adapter:			 {with_adapters[1]:13,d} ({with_adapters_fraction[1]:.1%})
        """)
    else:
        report += textwrap.dedent("""\
        Total reads processed:			 {N:13,d}
        Reads with adapters:			 {with_adapters[0]:13,d} ({with_adapters_fraction[0]:.1%})
        """)
    if stats["too_short"] is not None:
        report += "{pairs_or_reads} that were too short:		   {too_short:13,d} ({too_short_fraction:.1%})\n"
    if stats["too_long"] is not None:
        report += "{pairs_or_reads} that were too long:			   {too_long:13,d} ({too_long_fraction:.1%})\n"
    if stats["too_many_n"] is not None:
        report += "{pairs_or_reads} with too many N:			   {too_many_n:13,d} ({too_many_n_fraction:.1%})\n"

    report += textwrap.dedent("""\
    {pairs_or_reads} written (passing filters):			{written:13,d} ({written_fraction:.1%})

    Total basepairs processed:               {total_bp:13,d} bp
    """)
    if stats["paired"]:
        report += "	 Read 1: {total_bp1:13,d} bp\n"
        report += "	 Read 2: {total_bp2:13,d} bp\n"
    
    for modifier_class in trimmer_classes:
        name = modifier_class.__name__
        if stats[name] > 0:
            spacing = ' ' * (40 - len(modifier_class.display_str))
            report += "{0}:{1}{{{2}:13,d}} bp ({{{2}_fraction:.1%}})\n".format(modifier_class.display_str, spacing, name)
            if stats["paired"]:
                report += "	 Read 1: {{{}_bp[0]:13,d}} bp\n".format(name)
                report += "	 Read 2: {{{}_bp[1]:13,d}} bp\n".format(name)
    
    report += "Total written (filtered):                {total_written_bp:13,d} bp ({total_written_bp_fraction:.1%})\n"
    if stats["paired"]:
        report += "	 Read 1: {written_bp[0]:13,d} bp\n"
        report += "	 Read 2: {written_bp[1]:13,d} bp\n"
    
    try:
        report = report.format(**stats)
    except ValueError:
        # Python 2.6 does not support the comma format specifier (PEP 378)
        report = report.replace(",d}", "d}").format(**stats)
    _print(report)

    warning = False
    for which_in_pair in (0, 1):
        for adapter in stats["adapters"][which_in_pair]:
            if stats["paired"]:
                extra = 'First read: ' if which_in_pair == 0 else 'Second read: '
            else:
                extra = ''

            _print("=" * 3, extra + "Adapter", adapter["name"], "=" * 3)
            _print()
            if adapter["where"] == LINKED:
                _print("Sequence: {0}...{1}; Type: linked; Length: {2}+{3}; Trimmed: {4} times; Half matches: {5}".
                    format(
                        adapter["front_sequence"],
                        adapter["back_sequence"],
                        len(adapter["front_sequence"]),
                        len(adapter["back_sequence"]),
                        adapter["total_front"], adapter["total_back"]
                    ))
            else:
                _print("Sequence: {0}; Type: {1}; Length: {2}; Trimmed: {3} times.".
                    format(adapter["sequence"], ADAPTER_TYPES[adapter["where"]],
                        len(adapter["sequence"]), adapter["total"]))
    
            if adapter["total"] == 0:
                _print()
                continue
            if adapter["where"] == ANYWHERE:
                _print(adapter["total_front"], "times, it overlapped the 5' end of a read")
                _print(adapter["total_back"], "times, it overlapped the 3' end or was within the read")
                _print()
                print_error_ranges(len(adapter["sequence"]), adapter["max_error_rate"])
                _print("Overview of removed sequences (5')")
                print_histogram(adapter["lengths_front"], len(adapter["sequence"]), stats["N"],
                    adapter["max_error_rate"], adapter["errors_front"])
                _print()
                _print("Overview of removed sequences (3' or within)")
                print_histogram(adapter["lengths_back"], len(adapter["sequence"]), stats["N"],
                    adapter["max_error_rate"], adapter["errors_back"])
            elif adapter["where"] == LINKED:
                _print()
                print_error_ranges(len(adapter["front_sequence"]), adapter["front_max_error_rate"])
                print_error_ranges(len(adapter["back_sequence"]), adapter["back_max_error_rate"])
                _print("Overview of removed sequences at 5' end")
                print_histogram(adapter["front_lengths_front"], len(adapter["front_sequence"]),
                    stats["N"], adapter["front_max_error_rate"], adapter["front_errors_front"])
                _print()
                _print("Overview of removed sequences at 3' end")
                print_histogram(adapter["back_lengths_back"], len(adapter["back_sequence"]),
                    stats["N"], adapter["back_max_error_rate"], adapter["back_errors_back"])
            elif adapter["where"] in (FRONT, PREFIX):
                _print()
                print_error_ranges(len(adapter["sequence"]), adapter["max_error_rate"])
                _print("Overview of removed sequences")
                print_histogram(adapter["lengths_front"], len(adapter["sequence"]),
                    stats["N"], adapter["max_error_rate"], adapter["errors_front"])
            else:
                assert adapter["where"] in (BACK, SUFFIX)
                _print()
                print_error_ranges(len(adapter["sequence"]), adapter["max_error_rate"])
                warning = warning or print_adjacent_bases(adapter["adjacent_bases"], adapter["sequence"])
                _print("Overview of removed sequences")
                print_histogram(adapter["lengths_back"], len(adapter["sequence"]),
                    stats["N"], adapter["max_error_rate"], adapter["errors_back"])
    
    if warning:
        _print('WARNING:')
        _print('	One or more of your adapter sequences may be incomplete.')
        _print('	Please see the detailed output above.')
