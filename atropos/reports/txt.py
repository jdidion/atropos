N# coding: utf-8
"""Routines for printing a report. This is the legacy code for generating
text-based reports. This will eventually be deprecated in favor of the jinja
and multiqc reports.
"""
import math
import sys
import textwrap

# TODO: Fix https://github.com/marcelm/cutadapt/issues/128,
# https://github.com/marcelm/cutadapt/issues/112

INDENT = '    '
paragraph = textwrap.TextWrapper()
indented = textwrap.TextWrapper(initial_indent=INDENT, subsequent_indent=INDENT)

class HeaderPrinter(object):
    def __init__(self, outfile, levels=(('=', '='), ('-', None))):
        self.outfile = outfile
        self.levels = levels
    
    def __call__(self, *title, level=1, **kwargs):
        if level > len(self.levels):
            raise ValueError("Invalid level: {}".format(level))
        underline, overline = self.levels[level-1]
        if overline is True:
            overline = underline
        width = len(title)
        if overline:
            print(overline * width, file=self.outfile, **kwargs)
        print(*title, file=self.outfile, **kwargs)
        if underline:
            print(underline * width, file=self.outfile, **kwargs)

class RowPrinter(object):
    def __init__(
            self, outfile, colwidths=10, justification=('<', '>'), indent='',
            **kwargs):
        self.outfile = outfile
        self.colwidths = (colwidths,) if isinstance(colwidths, int) else colwidths
        self.justification = (justification,) if isinstance(justification, str) else justification
        self.indent = (indent,) if isinstance(indent, str) else indent
        self.print_args = kwargs
    
    def __call__(
            *args, colwidths=None, extra_width=None, justification=None,
            extra_justification=None, indent=None, extra_indent=None,
            underline=None, pct=False):
        ncols = len(args)
        if ncols == 0:
            print(file=self.outfile)
            return
        
        def adjust(a, extra=None):
            l = len(a)
            if l == ncols:
                return a
            elif l < ncols:
                return l[:ncols]
            else:
                return a + ([extra or a[-1]] * (l - ncols))
        
        colwidths, justification, indent = (
            adjust(a, extra)
            for a, extra in zip((
                    colwidths or self.colwidths,
                    justification or self.justification,
                    indent or self.indent
                ),
                (extra_width, extra_justification, extra_indent)))
        
        fmt_str = []
        fmt_args = []
        for i, (value, width, just, ind) in enumerate(
                zip(args, colwidths, justification, indent)):
            if isinstance(value, str):
                typ = 's'
                if len(value) > width:
                    value = textwrap.wrap(value, width)
            elif isinstance(value, float):
                typ = '.1' + '%' if pct else 'f'
            else:
                typ = ',d'
            fmt_str.append(ind + '{' + str(i) + ':' + just + str(width) + typ + '}')
            fmt_args.append(value)
        
        fmt_str = ' '.join(fmt_str)
        print(fmt_str.format(*fmt_args), file=self.outfile, **self.print_args)
        
        if underline:
            print(
                ' '.join(underline * width for width in colwidths),
                file=self.outfile, **self.print_args)

def generate_report(summary, outfile):
    close = False
    if outfile == 'stdout':
        outfile = sys.stdout
    elif outfile == 'stderr':
        outfile = sys.stderr
    else:
        outfile = open(outfile, "w")
        close = True
    
    try:
        print_report(summary, outfile)
    finally:
        if close:
            outfile.close()

def print_report(summary, outfile):
    """
    Args:
        summary: Summary dict.
        outfile: Open output stream.
    """
    
    paired = summary["paired"]
    pairs_or_reads = "Pairs" if paired else "Reads"
    total_bp = sum(summary['total_bp_counts'])
    max_width = len(str(total_bp))
    
    _print_header = HeaderPrinter(outfile)
    _print = RowPrinter(outfile, (25, max_width))
    
    N = summary["total_record_count"]
    if N == 0:
        _print(
            "No reads processed! Either your input file is empty or you "
            "used the wrong -f/--format parameter.")
        return
    
    timing = summary['timing']
    modifiers, filters, formatters = (
        summary['trim'][key]
        for key in ('modifiers', 'filters', 'formatters'))
    adapter_cutter = None
    if 'AdapterCutter' in modifiers:
        adapter_cutter = modifiers['AdapterCutter']
    elif 'InsertAdaptercutter' in modifiers:
        adapter_cutter = modifiers['InsertAdapterCutter']
    corrected = None
    trimmers = []
    for name, mod in modifers.items():
        if 'records_corrected' in mod:
            corrected = mod
        elif 'bp_trimmed' in mod:
            trimmers.append((name, mod))
    
    _print("Wallclock time: {0:.2F} s ({1:.0F} us/read; {2:.2F} M reads/minute).".format(
        timing["wallclock"],
        1E6 * summary["wallclock_time"] / N,
        N / timing["wallclock"] * 60 / 1E6))
    _print("CPU time (main process): {0:.2F} s".format(timing["cpu"]))
    
    _print_header("Summary", level=1)
    _print('', 'records', 'fraction', underline='-')
    _print("Total {} processed".format('read pairs' if paired else 'reads'), N)
    if adapter_cutter:
        for read in range(2 if paired else 1):
            _print(
                "Read{} with adapter".format(' {}'.format(read+1) if paired else 's'),
                adapter_cutter['records_with_adapters'][read],
                adapter_cutter['fraction_records_with_adapters'][read],
                indent='  ', pct=True)
    
    def _print_filter(name, sep):
        if name in filters:
            _print(
                "{} {} {}:".format(pairs_or_reads, sep, name.replace('_', ' ')),
                filters[name]['records_filtered'],
                filters[name]['fraction_records_filtered'], pct=True)
    
    _print_filter('too_short', 'that were')
    _print_filter('too_long', 'that were')
    _print_filter('too_many_n', 'with')
    
    _print(
        '{} written (passing filters)'.format(pairs_or_reads),
        formatters['records_written'], formatters['fraction_records_written'],
        pct=True)
    
    if corrected:
        _print(
            'Pairs corrected', corrected['records_corrected'],
            corrected['fraction_records_corrected'], pct=True)
    
    _print("Total bp processed", total_bp)
    if paired:
        for read in range(2):
            _print("Read {}".format(read+1), summary['total_bp_counts'][read],
            indent='  ')
    
    def _print_bp(title, d, key):
        _print(
            title, d['total_'.format(key)],
            d['fraction_total_{}'.format(key)], pct=True)
        if paired:
            for read in range(2):
                _print(
                    "Read {}".format(read+1), d[key][read],
                    d['fraction_{}'.format(key)][read], pct=True)
    
    for name, mod in trimmers:
        _print_bp(name, mod, 'bp_trimmed')
    
    _print_bp("Total bp written (filtered)", formatters, 'bp_written')
    
    if corrected:
        _print_bp("Total bp corrected", formatters, 'bp_corrected')
    
    if adapter_cutter:
        adapters = adapter_cutter['adapters']
        print_adapter_report(adapters, outfile, paired, N, max_width)

def print_adapter_report(adapters, outfile, paired, N, max_width):
    max_seq_len = max(
        (len(a['front_sequence'] + a['back_sequence'] + 3)
            if a['where']['name'] == 'linked' else len(a['sequence']))
        for pair in adapters
        for a in pair.values())
    
    _print_header = HeaderPrinter(outfile)
    _print = RowPrinter(outfile, (25, max_width))
    _print_seq = RowPrinter(outfile, (max_seq_len, 8, 3, max_width))
    
    def print_error_ranges(adapter_length, error_rate):
        _print("No. of allowed errors:")
        prev = 0
        for errors in range(1, int(error_rate * adapter_length) + 1):
            r = int(errors / error_rate)
            _print("{0}-{1} bp: {2};".format(prev, r - 1, errors - 1), end=' ')
            prev = r
        if prev == adapter_length:
            _print("{0} bp: {1}".format(
                adapter_length, int(error_rate * adapter_length)))
        else:
            _print("{0}-{1} bp: {2}".format(
                prev, adapter_length, int(error_rate * adapter_length)))
        _print()
    
    def print_histogram(d, adapter_length, n, error_rate, errors):
        """Print a histogram. Also, print the no. of reads expected to be
        trimmed by chance (assuming a uniform distribution of nucleotides in
        the reads).
        
        Args:
            d: A dictionary mapping lengths of trimmed sequences to their
                respective frequency
            adapter_length: adapter length
            n: total no. of reads.
            error_rate: Max error rate.
            errors: Histogram of actual numbers of errors.
        """
        def errs_to_str(l, e):
            if e in errors[l]:
                return str(errors[l][e])
            return "0"
        
        hist = []
        for length in sorted(d):
            # when length surpasses adapter_length, the
            # probability does not increase anymore
            estimate = n * 0.25 ** min(length, adapter_length)
            count = d[length]
            max_errors = max(errors[length].keys())
            errs = ' '.join(errs_to_str(length, e) for e in range(max_errors+1))
            hist.append((length, count, estimated, errors, errs))
        
        colwidths = [sizeof(*x) for x in zip(*h)]
        
        _print(
            "length", "count", "expect", "max.err", "error counts",
            colwidths=colwidths, underline='-')
        
        for h in hist:
            _print(*h, colwidths=colwidths)
        _print()
    
    def print_adjacent_bases(bases, sequence):
        """Print a summary of the bases preceding removed adapter sequences.
        Print a warning if one of the bases is overrepresented and there are
        at least 20 preceding bases available.
        
        Return:
            True if a warning was printed.
        """
        total = sum(bases.values())
        if total == 0:
            return False
        _print('Bases preceding removed adapters:')
        warnbase = None
        for base in ['A', 'C', 'G', 'T', '']:
            b = base if base != '' else 'none/other'
            fraction = 1.0 * bases[base] / total
            _print(b, fraction)
            if fraction > 0.8 and base != '':
                warnbase = b
        if total >= 20 and warnbase is not None:
            _print('WARNING:')
            _print(indented.wrap(
                'The adapter is preceded by "{0}" extremely often. The provided '
                'adapter sequence may be incomplete. To fix the problem, add '
                '"{0}" to the beginning of the adapter sequence.'.format(warnbase)))
            _print()
            return True
        _print()
        return False
        
    warning = False
    for pair in range(2 if paired else 1):
        header = "Adapter {}"
        if paired:
            header = ('First read: ' if pair == 0 else 'Second read: ') + header
        
        for name, adapter in adapters[pair]:
            _print_header(header.format(adapter["name"]), level=1)
            
            where_name = adapter["where"]["name"]
            if where_name == "linked":
                front_len, back_len = [
                    len(adapter[s])
                    for s in ('front_sequence', 'back_sequence')]
                _print_seq(
                    "Sequence", "Type", "Length", "Trimmed (x)", "Half matches (x)",
                    underline='-')
                _print_seq(
                    "{0}...{1}".format(
                        adapter["front_sequence"],
                        adapter["back_sequence"]),
                    "linked",
                    "{2}+{3}".format(front_len, back_len),
                    adapter["total_front"],
                    adapter["total_back"]
                )
            else:
                seq_len = len(adapter["sequence"])
                _print_seq(
                    "Sequence", "Type", "Length", "Trimmed (x)", underline='-')
                _print_seq(
                    adapter["sequence"], adapter["where"]["desc"], seq_len,
                    adapter["total"])
            
            _print()
            
            if adapter["total"] == 0:
                return
            
            if where_name == "anywhere":
                _print(adapter["total_front"], "times, it overlapped the 5' end of a read")
                _print(adapter["total_back"], "times, it overlapped the 3' end or was within the read")
                _print()
                print_error_ranges(seq_len, adapter["max_error_rate"])
                _print("Overview of removed sequences (5')")
                hist_args = (seq_len, N, adapter["max_error_rate"], adapter["errors_front"])
                print_histogram(adapter["lengths_front"], *hist_args)
                _print()
                _print("Overview of removed sequences (3' or within)")
                print_histogram(adapter["lengths_back"], *hist_args)
            
            elif where_name == "linked":
                print_error_ranges(front_len, adapter["front_max_error_rate"])
                print_error_ranges(back_len, adapter["back_max_error_rate"])
                _print("Overview of removed sequences at 5' end")
                print_histogram(
                    adapter["front_lengths_front"], front_len, N,
                    adapter["front_max_error_rate"], adapter["front_errors_front"])
                _print()
                _print("Overview of removed sequences at 3' end")
                print_histogram(
                    adapter["back_lengths_back"], back_len, N,
                    adapter["back_max_error_rate"], adapter["back_errors_back"])
            
            elif where_name in ("front", "prefix"):
                print_error_ranges(seq_len, adapter["max_error_rate"])
                _print("Overview of removed sequences")
                print_histogram(
                    adapter["lengths_front"], seq_len, N,
                    adapter["max_error_rate"], adapter["errors_front"])
            
            elif where_name in ("back", "suffix"):
                print_error_ranges(seq_len, adapter["max_error_rate"])
                warning = warning or print_adjacent_bases(
                    adapter["adjacent_bases"], adapter["sequence"])
                _print("Overview of removed sequences")
                print_histogram(
                    adapter["lengths_back"], seq_len, N,
                    adapter["max_error_rate"], adapter["errors_back"])
    
    if warning:
        _print('WARNING:')
        _print(indented.wrap(
            'One or more of your adapter sequences may be incomplete. '
            'Please see the detailed output above.'))

def print_read_stats(options, stats):
    outfile = options.output
    close = False
    
    if outfile is not None:
        outfile = open(outfile, "w")
        close = True
    elif not options.quiet:
        outfile = sys.stderr if options.output is None else sys.stdout
    else:
        return None
    
    try:
        #try:
        #    generate_mako(stats, outfile, "stats.txt")
        #except:
        generate_read_stats(stats, outfile)
        return stats
    finally:
        if close:
            outfile.close()

def generate_read_stats(stats, outfile):
    def _print_stats(title, data):
        
        paired = 'read2' in data
        max_width = len(str(max(data['read1']['count'], data['read2']['count'])))
        # add space for commas and column separation
        max_width += (max_width // 3) + 1
        
        _print_header = HeaderPrinter(outfile)
        _print = RowPrinter(outfile, (25, max_width))
        
        def _print_histogram(title, hist1, hist2=None):
            _print_header(title, 2)
            if hist2:
                hist = ((h1[0], h1[1], h2[1]) for h1, h2 in zip(hist1, hist2))
            else:
                hist = h1
            for h in hist:
                _print(*h)
        
        def _print_base_histogram(title, col_names, hist):
            _print_header(title, 2)
            _print('Position', *col_names, extra_width=4)
            for row in hist:
                total_count = sum(row[1:])
                base_pcts = (
                    round(count * 100 / total_count, 1)
                    for count in row[1:])
                _print(row[0], *base_pcts, extra_width=4)
        
        def _print_tile_histogram(
                title, value_cols, hist, index_cols=('Tile',), index_widths=(5,)):
            _print_header(title, 2)
            ncol = len(value_cols)
            # As far as I know, tiles are max 4 digits
            max_width = max(
                4, len(str(math.ceil(data['read1']['count'] / ncol)))) + 1
            widths = index_widths + ((max_width,) * ncol)
            _print(*(index_cols + value_cols), extra_width=max_width)
            for row in hist:
                _print(*row, extra_width=max_width)
        
        _print_header(title, 1)
        _print('', 'Read1', 'Read2')
        
        # Sequence-level stats
        _print(
            "Read pairs:" if paired else "Reads:",
            data['read1']['count'],
            data['read2']['count'])
        _print()
        _print_histogram(
            "Sequence lengths:",
            data['read1']['length'],
            data['read2']['length'])
        _print()
        if 'qualities' in data['read1']:
            _print_histogram(
                "Sequence qualities:",
                data['read1']['qualities'],
                data['read2']['qualities'])
            _print()
        _print_histogram(
            "Sequence GC content (%):",
            data['read1']['gc'],
            data['read2']['gc'])
        _print()
        
        if 'tile_sequence_qualities' in data['read1']:
            _print_tile_histogram(
                "Per-tile sequence qualities (Read 1):",
                *data['read1']['tile_sequence_qualities'])
            _print()
            _print_tile_histogram(
                "Per-tile sequence qualities (Read 2):",
                *data['read2']['tile_sequence_qualities'])
            _print()
        
        # Base-level stats
        if 'base_qualities' in data['read1']:
            _print_base_histogram(
                "Base qualities (Read 1):",
                *data['read1']['base_qualities'])
            _print()
            _print_base_histogram(
                "Base qualities (Read 2):",
                *data['read2']['base_qualities'])
            _print()
        _print_base_histogram(
            "Base composition (Read 1):",
            *data['read1']['bases'])
        _print()
        _print_base_histogram(
            "Base composition (Read 2):",
            *data['read2']['bases'])
        _print()
        if 'tile_base_qualities' in data['read1']:
            _print_tile_histogram(
                "Per-tile base qualities (Read 1):",
                *data['read1']['tile_base_qualities'],
                index_cols=('Position', 'Tile'),
                index_widths=(25, 5))
            _print()
            _print_tile_histogram(
                "Per-tile base qualities (Read 2):",
                *data['read2']['tile_base_qualities'],
                index_cols=('Position', 'Tile'),
                index_widths=(25, 5))
            _print()
    
    if 'pre' in stats:
        _print_stats('Pre-trimming stats', stats['pre'])
    if 'post' in stats:
        _print_stats('Post-trimming stats', stats['post'])

def sizeof(*x):
    if isinstance(x[0], str):
        return max(len(s) for s in x)
    elif isinstance(x[0], int):
        return len(str(max(x)))
    else:
        return len(str(round(max(x), 1)))
