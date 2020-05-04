# coding: utf-8
"""Routines for printing a text report. This is the legacy code for generating
text-based reports. This will eventually be deprecated in favor of the Jinja2
and MultiQC reports.
"""
from abc import ABCMeta, abstractmethod
import math
import textwrap

from atropos.io import open_output
from atropos.util import truncate_string, weighted_median
from .reports import BaseReportGenerator

INDENT = '  '
PARAGRAPH = textwrap.TextWrapper()
INDENTED = textwrap.TextWrapper(initial_indent=INDENT, subsequent_indent=INDENT)

class Printer(object):
    """Manages printing to a file.

    Args:
        outfile: The output file.
        kwargs: Additional keyword arguments passed to the print function.
    """
    def __init__(self, outfile, indent=None, **kwargs):
        self.outfile = outfile
        self.indent = indent
        self.print_args = kwargs

    def __call__(self, *args, indent=None, **kwargs):
        if isinstance(indent, int):
            indent = self.indent * indent
        else:
            indent = indent or self.indent
        if indent:
            self._print(indent, end='')
        self._print(*args, **kwargs)

    def _print(self, *args, **kwargs):
        if self.print_args:
            print_args = self.print_args.copy()
            print_args.update(kwargs)
        else:
            print_args = kwargs
        print(*args, file=self.outfile, **print_args)

    def newline(self):
        """Print a newline.
        """
        print(file=self.outfile)

class TitlePrinter(Printer):
    """Printer that formats titles.

    Args:
        outfile: The output file.
        levels: The formatting associated with different header levels.
        kwargs: Additional keyword arguments passed to the print function.
    """
    def __init__(
            self, outfile,
            levels=(('=', '='), ('-', '-'), ('-', None), ('~', None)),
            **kwargs):
        super().__init__(outfile, **kwargs)
        self.levels = levels

    def __call__(self, *title, level=None, newline=True, **kwargs):
        title = ' '.join(title)

        if level is not None:
            if level >= len(self.levels):
                raise ValueError("Invalid level: {}".format(level))
            underline, overline = self.levels[level]
            if overline is True:
                overline = underline
            width = len(title)
            if overline:
                self._print(overline * width, **kwargs)

        self._print(title, **kwargs)

        if level is not None and underline:
            self._print(underline * width, **kwargs)

        if newline:
            self.newline()

class RowPrinter(Printer):
    """Priter that formats rows in a table.

    Args:
        outfile: The output file.
        colwidths: Column widths.
        justification: Column justifiations ('<' for left, '>' for right).
        indent: Column indents.
        pct: Whether floats should be formatted as percentages.
        default: Default value for None's.
        kwargs: Additional keyword arguments passed to the print function.

    colwidths, justification, and indent can be longer or shorter than the
    number of arguments; if shorter, the last value in the list is repeated;
    if longer, the list is truncated.
    """
    def __init__(
            self, outfile, colwidths=10, justification=('<', '>'), indent='',
            pct=False, default=0, **kwargs):
        super().__init__(outfile, **kwargs)
        self.colwidths, self.justification, self.indent = (
            (arg,) if isinstance(arg, typ) else tuple(arg)
            for arg, typ in zip(
                (colwidths, justification, indent),
                (int, str, str)))
        self.pct = pct
        self.default = default

    def print_rows(self, *rows, header=None, **kwargs):
        """Print multiple rows. Automatically figures out column widths.

        Args:
            rows: Rows to print.
            header: Header row.
            kwargs: Additional keyword arguments to self.__call__.
        """
        colwidths = tuple(sizeof(*x) for x in zip(*rows))
        if header:
            if isinstance(header[0], str):
                header_widths = (sizeof(h) for h in header)
                header_rows = [header]
            else:
                header_widths = (
                    max(sizeof(h) for h in header_part)
                    for header_part in header)
                header_rows = list(zip(*header))
            colwidths = tuple(
                max(h, c)
                for h, c in zip(header_widths, colwidths))
            for i, header_row in enumerate(header_rows, 1):
                self(*header_row, colwidths=colwidths, header=(i == len(header_rows)), **kwargs)
        for row in rows:
            self(*row, colwidths=colwidths)

    def __call__(
            self, *args, colwidths=None, extra_width=None, justification=None,
            extra_justification=None, indent=None, extra_indent=None,
            header=False, underline='-', pct=None, default=None, **kwargs):
        """Print a row.

        Args:
            args: Fields in the row.
            colwidths, justification, indent: Row-specific colwidths,
                justification, indent.
            extra_width, extra_justification, extra_indent: colwidth/
                justification/indent to use for extra fields.
            header: Whether this is a header row.
            underline: Whether to use an underline after the header row. Either
                a bool or a character.
            pct: Whether floating point values should be formatted as
                percentages.
            default: Default value.
            kwargs: Additional keyword arguments to pass to print.
        """
        ncols = len(args)
        if ncols == 0:
            self.newline()
            return

        if pct is None:
            pct = self.pct

        def adjust(arr, extra=None):
            """Adjust an array. If longer than the number of columns,
            truncate; if shorter, fill in by repeating the last element.
            """
            alen = len(arr)
            if alen == ncols:
                return arr
            elif alen > ncols:
                return arr[:ncols]
            else:
                return arr + ((extra or arr[-1],) * (ncols - alen))

        colwidths, justification, indent = (
            adjust(arr, extra)
            for arr, extra in zip(
                (
                    colwidths or self.colwidths,
                    justification or self.justification,
                    indent or self.indent),
                (extra_width, extra_justification, extra_indent)))

        # adjust colwidths if this is a header
        if header:
            colwidths = tuple(
                max(w, len(str(a)))
                for w, a in zip(colwidths, args))

        fmt_str = []
        fmt_args = []
        for i, (value, width, just, ind) in enumerate(
                zip(args, colwidths, justification, indent)):
            if value is None:
                value = default or self.default
            if isinstance(value, str):
                typ = 's'
                if len(value) > width:
                    value = truncate_string(value, width)
            elif isinstance(value, float):
                typ = ',.1' + ('%' if pct else 'f')
            else:
                typ = ',d'
            fmt_str.append(
                ind + '{' + str(i) + ':' + just + str(width-len(ind)) +
                typ + '}')
            fmt_args.append(value)

        fmt_str = ' '.join(fmt_str)
        self._print(fmt_str.format(*fmt_args), **kwargs)

        if header:
            sepline = ' '.join((underline * width) for width in colwidths)
            self._print(sepline, **kwargs)

class LegacyReportGenerator(BaseReportGenerator):
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == 'txt':
            with open_output(outfile, context_wrapper=True) as out:
                generate_report(summary, out)
        else:
            super().generate_from_template(fmt, summary, outfile, **kwargs)

def generate_report(summary, outfile):
    """Generate a report.

    Args:
        summary: The summary dict.
        outfile: The output file name/object.
    """

    print_summary_report(summary, outfile)
    if 'trim' in summary:
        print_trim_report(summary, outfile)
    if 'pre' in summary:
        print_pre_trim_report(summary, outfile)
    if 'post' in summary:
        print_post_trim_report(summary, outfile)

def print_summary_report(summary, outfile):
    """Print the top-level summary report.

    Args:
        summary: The summary dict.
        outfile: The output file object.
    """
    _print_title = TitlePrinter(outfile)
    _print = Printer(outfile)

    _print_title("Atropos", level=0)
    _print("Atropos version: {}".format(summary['version']))
    _print("Python version: {}".format(summary['python']))
    _print("Command line parameters: {} {}".format(
        summary['command'], " ".join(summary['options']['orig_args'])))
    _print()

    _print("Sample ID: {}".format(summary['sample_id']))
    _print("Input format: {}".format(summary['derived']['input_format']))
    _print("Input files:")
    for infile in summary['input']['input_names']:
        if infile is not None:
            _print(infile, indent=INDENT)
    _print()

    timing = summary['timing']
    total = summary["total_record_count"]
    wctime = ["Wallclock time: {:.2F} s".format(timing["wallclock"])]
    if total > 0:
        wctime.append("({0:.0F} us/read; {1:.2F} M reads/minute)".format(
            1E6 * timing["wallclock"] / total,
            total / timing["wallclock"] * 60 / 1E6))
    _print("Start time: {}".format(timing['start']))
    _print(*wctime)
    _print("CPU time (main process): {0:.2F} s".format(timing['cpu']))
    _print()

def print_trim_report(summary, outfile):
    """Print the trimming report.

    Args:
        summary: Summary dict.
        outfile: Open output stream.
    """
    paired = summary["options"]["paired"]
    pairs_or_reads = "Pairs" if paired else "Reads"
    total_bp = sum(summary['total_bp_counts'])
    max_width = len(str(total_bp))
    # account for commas
    max_width += (max_width // 3)

    _print_title = TitlePrinter(outfile)
    _print = RowPrinter(outfile, (35, max_width))

    total = summary["total_record_count"]
    if total == 0:
        _print_error = Printer(outfile)
        _print_error(
            "No reads processed! Either your input file is empty or you "
            "used the wrong -f/--format parameter.")
        return

    modifiers, filters, formatters = (
        summary['trim'][key]
        for key in ('modifiers', 'filters', 'formatters'))
    adapter_cutter = None
    error_corrector = None
    for modifier_dict in modifiers.values():
        if adapter_cutter is None and 'adapters' in modifier_dict:
            adapter_cutter = modifier_dict
            break
        if error_corrector is None and 'bp_corrected' in modifier_dict:
            error_corrector = modifier_dict
    correction_enabled = summary["options"]["correct_mismatches"]
    corrected = None
    trimmers = []
    for name, mod in modifiers.items():
        if 'bp_trimmed' in mod:
            trimmers.append((name, mod))
        if correction_enabled and 'records_corrected' in mod:
            corrected = mod

    _print_title("Trimming", level=1)
    _print(pairs_or_reads, 'records', 'fraction', header=True)
    _print(
        "Total {} processed:".format('read pairs' if paired else 'reads'),
        total)
    if adapter_cutter:
        if paired:
            for read in range(2):
                _print(
                    "Read {} with adapter:".format(read+1),
                    adapter_cutter['records_with_adapters'][read],
                    adapter_cutter['fraction_records_with_adapters'][read],
                    indent=(INDENT, ''), pct=True)
        else:
            _print(
                "Reads with adapters:",
                adapter_cutter['records_with_adapters'][0],
                adapter_cutter['fraction_records_with_adapters'][0],
                pct=True)

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
        "{} written (passing filters):".format(pairs_or_reads),
        formatters['records_written'], formatters['fraction_records_written'],
        pct=True)

    if corrected:
        _print(
            "Pairs corrected:", corrected['records_corrected'],
            corrected['fraction_records_corrected'], pct=True)

    _print()
    _print("Base pairs", 'bp', 'fraction', header=True)

    _print("Total bp processed:", total_bp)
    if paired:
        for read in range(2):
            _print(
                "Read {}:".format(read+1), summary['total_bp_counts'][read],
                indent=(INDENT, ''))

    def _print_bp(title, data, key, default=0):
        if paired:
            _print(
                title, data['total_{}'.format(key)],
                data['fraction_total_{}'.format(key)], pct=True)
            for read in range(2):
                _print(
                    "Read {}:".format(read+1),
                    data[key][read],
                    data['fraction_{}'.format(key)][read],
                    indent=(INDENT, ''), pct=True, default=default)
        else:
            _print(
                title, data[key][0], data['fraction_{}'.format(key)][0],
                pct=True, default=default)

    for name, mod in trimmers:
        _print_bp(mod['desc'], mod, 'bp_trimmed')

    _print_bp("Total bp written (filtered):", formatters, 'bp_written')

    if error_corrector:
        _print_bp("Total bp corrected:", error_corrector, 'bp_corrected')

    if adapter_cutter:
        _print()
        adapters = adapter_cutter['adapters']
        print_adapter_report(adapters, outfile, paired, total, max_width)

def print_adapter_report(adapters, outfile, paired, total_records, max_width):
    """Print details for a adapters.

    Args:
        adapters: Sequence of adapter info dicts.
        outfile: The output file
        paired: Whether the data is paired-end.

        max_width: Max column width.
    """
    adapter_lenghts = []
    for pair in adapters:
        if pair:
            for adapter in pair.values():
                if adapter['where']['name'] == 'linked':
                    adapter_lenghts.append(3 + len(
                        adapter['front_sequence'] +
                        adapter['back_sequence']))
                else:
                    adapter_lenghts.append(len(adapter['sequence']))
    max_seq_len = max(adapter_lenghts)

    _print = Printer(outfile)
    _print_title = TitlePrinter(outfile)
    _print_adj = RowPrinter(outfile, (12, 5), pct=True, indent=(INDENT, ''))

    seq_printer = RowPrinter(
        outfile, (max_seq_len, 14, 3, max_width), ('<', '<', '>'))
    hist_printer = RowPrinter(outfile, justification=('>', '>', '>', '>', '<'))

    def print_error_ranges(adapter_length, error_rate):
        """Print max number of errors for ranges of adapter match lengths.
        """
        _print("No. of allowed errors:")
        prev = 0
        for errors in range(1, int(error_rate * adapter_length) + 1):
            range_start = int(errors / error_rate)
            _print(
                "{0}-{1} bp: {2};".format(prev, range_start - 1, errors - 1),
                end=' ')
            prev = range_start
        if prev == adapter_length:
            _print("{0} bp: {1}".format(
                adapter_length, int(error_rate * adapter_length)))
        else:
            _print("{0}-{1} bp: {2}".format(
                prev, adapter_length, int(error_rate * adapter_length)))
        _print()

    def print_histogram(
            data, adapter_length, num_reads, error_rate, errors,
            match_probabilities):
        """Print a histogram. Also, print the no. of reads expected to be
        trimmed by chance (assuming a uniform distribution of nucleotides in
        the reads).

        Args:
            data: A dictionary mapping lengths of trimmed sequences to their
                respective frequency
            adapter_length: adapter length
            num_reads: total no. of reads.
            error_rate: Max error rate.
            errors: Histogram of actual numbers of errors.
            match_probabilities: List of random match proabilities for each
                position in the adapter seqeunce.
        """
        hist = []
        hist_errors = []
        for length, count in data.items():
            # when length surpasses adapter_length, the probability does not
            # increase anymore
            estimated = (
                num_reads * match_probabilities[min(length, adapter_length)])
            hist.append([
                length, count, estimated,
                int(error_rate * min(length, adapter_length))])
            hist_errors.append(errors['rows'][length])

        col_sizes = [len(str(max(col))) for col in zip(*hist_errors)]

        def _format_hist_errors(errs):
            trailing = True
            hist_str = []
            for i, e in reversed(list(enumerate(errs))):
                if trailing and e == 0:
                    pass
                else:
                    trailing = False
                    hist_str.append(('{:<' + str(col_sizes[i]) + 'd}').format(e))
            return ' '.join(reversed(hist_str))

        for i, errs in enumerate(hist_errors):
            hist[i].append(_format_hist_errors(errs))

        error_header = ' '.join(
            ('{:<' + str(e) + 'd}').format(i)
            for i, e in enumerate(col_sizes))

        hist_printer.print_rows(
            *hist,
            header=(
                ("length",""), ("count",""), ("expect",""), ("max.err",""),
                ("error counts", error_header)))
        hist_printer.newline()

    def print_adjacent_bases(bases):
        """Print a summary of the bases preceding removed adapter sequences.
        Print a warning if one of the bases is overrepresented and there are
        at least 20 preceding bases available.

        Return:
            True if a warning was printed.
        """
        total = sum(bases.values())
        if total == 0:
            return False
        _print("Bases preceding removed adapters:")
        warnbase = None
        for base in ['A', 'C', 'G', 'T', '']:
            base_label = base if base != '' else 'none/other'
            fraction = 1.0 * bases[base] / total
            _print_adj(base_label, fraction)
            if fraction > 0.8 and base != '':
                warnbase = base_label
        if total >= 20 and warnbase is not None:
            _print('WARNING:')
            _print("\n".join(INDENTED.wrap(
                'The adapter is preceded by "{0}" extremely often. The '
                'provided adapter sequence may be incomplete. To fix the '
                'problem, add "{0}" to the beginning of the adapter '
                'sequence.'.format(warnbase))))
            _print()
            return True
        _print()
        return False

    warning = False
    for pair in range(2 if paired else 1):
        if adapters[pair] is None:
            continue

        header = "Adapter {}"
        if paired:
            header = ("First read: " if pair == 0 else "Second read: ") + header

        for name, adapter in adapters[pair].items():
            _print_title(header.format(name), level=1)

            where_name = adapter["where"]["name"]
            if where_name == "linked":
                front_len, back_len = [
                    len(adapter[s])
                    for s in ('front_sequence', 'back_sequence')]
                seq_printer.print_rows(
                    (
                        "{}...{}".format(
                            adapter["front_sequence"],
                            adapter["back_sequence"]),
                        "linked",
                        "{}+{}".format(front_len, back_len),
                        adapter["total_front"],
                        adapter["total_back"]),
                    header=(
                        "Sequence", "Type", "Length", "Trimmed (x)",
                        "Half matches (x)"))
            else:
                seq_len = len(adapter["sequence"])
                seq_printer.print_rows(
                    (
                        adapter["sequence"], adapter["where"]["desc"], seq_len,
                        adapter["total"]),
                    header=(
                        "Sequence", "Type", "Length", "Trimmed (x)"))

            _print()

            if adapter["total"] == 0:
                continue

            if where_name == "anywhere":
                _print(
                    adapter["total_front"],
                    "times, it overlapped the 5' end of a read")
                _print(
                    adapter["total_back"],
                    "times, it overlapped the 3' end or was within the read")
                _print()
                print_error_ranges(seq_len, adapter["max_error_rate"])
                _print("Overview of removed sequences (5'):")
                print_histogram(
                    adapter["lengths_front"], seq_len, total_records,
                    adapter["max_error_rate"], adapter["errors_front"],
                    adapter['match_probabilities'])
                _print()
                _print("Overview of removed sequences (3' or within):")
                print_histogram(
                    adapter["lengths_back"], seq_len, total_records,
                    adapter["max_error_rate"], adapter["errors_back"],
                    adapter['match_probabilities'])

            elif where_name == "linked":
                print_error_ranges(front_len, adapter["front_max_error_rate"])
                print_error_ranges(back_len, adapter["back_max_error_rate"])
                _print("Overview of removed sequences at 5' end:")
                print_histogram(
                    adapter["front_lengths_front"], front_len, total_records,
                    adapter["front_max_error_rate"],
                    adapter["front_errors_front"],
                    adapter['front_match_probabilities'])
                _print()
                _print("Overview of removed sequences at 3' end:")
                print_histogram(
                    adapter["back_lengths_back"], back_len, total_records,
                    adapter["back_max_error_rate"], adapter["back_errors_back"],
                    adapter['back_match_probabilities'])

            elif where_name in ("front", "prefix"):
                print_error_ranges(seq_len, adapter["max_error_rate"])
                _print("Overview of removed sequences:")
                print_histogram(
                    adapter["lengths_front"], seq_len, total_records,
                    adapter["max_error_rate"], adapter["errors_front"],
                    adapter['match_probabilities'])

            elif where_name in ("back", "suffix"):
                print_error_ranges(seq_len, adapter["max_error_rate"])
                warning = warning or print_adjacent_bases(
                    adapter["adjacent_bases"])
                _print("Overview of removed sequences:")
                print_histogram(
                    adapter["lengths_back"], seq_len, total_records,
                    adapter["max_error_rate"], adapter["errors_back"],
                    adapter['match_probabilities'])

    if warning:
        _print('WARNING:')
        _print("\n".join(INDENTED.wrap(
            'One or more of your adapter sequences may be incomplete. '
            'Please see the detailed output above.')))

def print_pre_trim_report(summary, outfile):
    """Print pre-trimming stats.

    Args:
        summary: The summary dict.
        outfile: The output file.
    """
    pre = summary['pre']
    _print_title = TitlePrinter(outfile)
    _print = Printer(outfile)
    _print_title("Pre-trimming stats", level=1)
    for source, data in pre.items():
        _print_title("Source", level=3, newline=False)
        # TODO: When multi-file input is supported, this code will need to
        # get summary['input']['input_names'][source]
        for read, src in enumerate(summary['input']['input_names'], 1):
            if src is None:
                continue
            _print("Read {}: {}".format(read, src))
        _print()
        print_stats_report(data, outfile)

def print_post_trim_report(summary, outfile):
    """Print post-trimming stats.

    Args:
        summary: The summary dict.
        outfile: The output file.
    """
    post = summary['post']
    _print_title = TitlePrinter(outfile)
    _print = Printer(outfile)
    _print_title("Post-trimming stats", level=1)
    for dest, stats in post.items():
        _print_title("Destination: {}".format(dest), level=2)
        for source, data in stats.items():
            _print_title("Source", level=3, newline=False)
            # TODO: When multi-file input is supported, this code will need to
            # get summary['input']['input_names'][source]
            for read, src in enumerate(summary['input']['input_names'], 1):
                if src is None:
                    continue
                _print("Read {}: {}".format(read, src))
            _print()
            print_stats_report(data, outfile)


class StatsPrinter(metaclass=ABCMeta):
    def __init__(self, data, outfile):
        self._data = data
        self._title_printer = TitlePrinter(outfile)
        max_count = self._max_count()
        max_width = len(str(max_count))
        # add space for commas and column separation
        max_width += (max_width // 3) + 1
        self._printer = RowPrinter(outfile, (35, max_width))

    @abstractmethod
    def _max_count(self):
        pass

    def _print_histogram(self, title, hist1, hist2=None):
        self._title_printer(title, level=2)
        if hist1 is None:
            self._printer("No Data")
            return
        if hist2:
            hist = (
                (key, hist1.get(key, 0), hist2.get(key, 0))
                for key in sorted(set(hist1.keys()) | set(hist2.keys())))
        else:
            hist = sorted(hist1.items(), key=lambda x: x[0])
        for histbin in hist:
            self._printer(*histbin)

    def _print_base_histogram(self, title, hist, extra_width=4, index_name='Pos'):
        self._title_printer(title, level=2)
        if hist is None:
            self._printer("No Data")
            return
        self._printer(
            index_name, *hist['columns'], header=True, extra_width=extra_width)
        for pos, row in hist['rows'].items():
            total_count = sum(row)
            base_pcts = (
                round(count * 100 / total_count, 1)
                for count in row)
            self._printer(pos, *base_pcts, extra_width=extra_width)

    def _print_tile_histogram(self, title, hist):
        if hist is None:
            self._title_printer(title, level=2)
            self._printer("No Data")
            return
        ncol = len(hist['columns'])
        max_tile_width = max(
            4, len(str(math.ceil(self._data['read1']['counts'] / ncol)))) + 1
        self._print_base_histogram(
            title, hist, extra_width=max_tile_width, index_name='Tile')

    def _print_tile_base_histogram(self, title, hist):
        """Print a histogram of position x tile, with values as the median
        base quality.
        """
        self._title_printer(title, level=2)
        if hist is None:
            self._printer("No Data")
            return
        quals = hist['columns']
        tiles = hist['columns2']
        ncol = len(tiles)
        max_tile_width = max(
            4, len(str(math.ceil(self._data['read1']['counts'] / ncol)))) + 1
        self._printer('Pos', *tiles, header=True, extra_width=max_tile_width)
        for pos, tiles in hist['rows'].items():
            # compute the weighted median for each tile at each position
            self._printer(
                pos,
                *(
                    weighted_median(quals, tile_counts)
                    for tile_counts in tiles.values()),
                extra_width=max_tile_width)

    @abstractmethod
    def print_header(self):
        pass

    @abstractmethod
    def print_counts(self):
        pass

    @abstractmethod
    def print_histogram(self, title, key1, key2):
        pass

    @abstractmethod
    def print_tile_histograms(self, title, key):
        pass

    @abstractmethod
    def print_base_histograms(self, title, key):
        pass

    @abstractmethod
    def print_tile_base_histograms(self, title, key):
        pass


class SingleEndStatsPrinter(StatsPrinter):
    def _max_count(self):
        return self._data['read1']['counts']

    def print_header(self):
        self._printer('', 'Read1', header=True)

    def print_counts(self):
        self._printer("Reads:", self._data['read1']['counts'])
        self._printer()

    def print_histogram(self, title, key1, key2):
        if key1 in self._data['read1']:
            self._print_histogram(title, self._data['read1'][key1][key2])
            self._printer()

    def print_tile_histograms(self, title, key):
        if key in self._data['read1']:
            self._print_tile_histogram(
                "Read 1 {}".format(title),
                self._data['read1'][key])
            self._printer()

    def print_base_histograms(self, title, key):
        if key in self._data['read1']:
            self._print_base_histogram(
                "Read 1 {}".format(title),
                self._data['read1'][key])
            self._printer()

    def print_tile_base_histograms(self, title, key):
        if key in self._data['read1']:
            self._print_tile_base_histogram(
                "Read 1 {}".format(title),
                self._data['read1'][key])


class PairedEndStatsPrinter(StatsPrinter):
    def _max_count(self):
        return max(self._data['read1']['counts'], self._data['read2']['counts'])

    def print_header(self):
        self._printer('', 'Read1', 'Read2', header=True)

    def print_counts(self):
        self._printer(
            "Read pairs:",
            self._data['read1']['counts'],
            self._data['read2']['counts'])
        self._printer()

    def print_histogram(self, title, key1, key2):
        if key1 in self._data:
            self._print_histogram(
                title,
                self._data['read1'][key1][key2],
                self._data['read2'][key1][key2])
            self._printer()

    def print_tile_histograms(self, title, key):
        if 'tile_sequence_qualities' in self._data['read1']:
            self._print_tile_histogram(
                "Read 1 {}".format(title),
                self._data['read1'][key])
            self._printer()
            self._print_tile_histogram(
                "Read 2 {}".format(title),
                self._data['read2'][key])
            self._printer()

    def print_base_histograms(self, title, key):
        if key in self._data:
            self._print_base_histogram(
                "Read 1 {}".format(title),
                self._data['read1'][key])
            self._printer()
            self._print_base_histogram(
                "Read 2 {}".format(title),
                self._data['read2'][key])
            self._printer()

    def print_tile_base_histograms(self, title, key):
        if key in self._data['read1']:
            self._print_tile_base_histogram(
                "Read 1 {}".format(title),
                self._data['read1'][key])
            self._printer()
            self._print_tile_base_histogram(
                "Read 2 {}".format(title),
                self._data['read2'][key])
            self._printer()


def print_stats_report(data, outfile):
    """Print stats.

    Args:
        data: The stats dict.
        outfile: The output file.
    """
    paired = 'read2' in data
    if paired:
        printer = PairedEndStatsPrinter(data, outfile)
    else:
        printer = SingleEndStatsPrinter(data, outfile)

    printer.print_header()

    # Sequence-level stats
    printer.print_counts()
    printer.print_histogram("Sequence lengths:", "lengths", "hist")
    printer.print_histogram("Sequence qualities:", "qualities", "hist")
    printer.print_histogram("Sequence GC content (%)", "gc", "hist")
    printer.print_tile_histograms(
        "per-tile sequence qualities (%)", 'tile_sequence_qualities')

    # Base-level stats
    printer.print_base_histograms("base qualities (%)", 'base_qualities')
    printer.print_base_histograms("base composition (%)", 'bases')
    printer.print_tile_base_histograms(
        "per-tile base qualities (%)", 'tile_base_qualities')


def sizeof(*x, seps=True, prec=1):
    """Returns the largest string size of all objects in x, where x is a
    sequence of string or numeric values.

    Args:
        *x: The objects to test.
        seps: Whether to include separators (,.) in the size.
        prec: Precision of float values.
    """
    if isinstance(x[0], str):
        return max(len(s) for s in x)
    else:
        if isinstance(x[0], int):
            numlen = len(str(max(x)))
            if seps:
                numlen += (numlen // 3)
        elif isinstance(x[0], float):
            numlen = len(str(round(max(x), prec)))
            if seps:
                numlen += ((numlen - prec - 1) // 3)
        else:
            raise ValueError("Unexpected data type: {}".format(x[0].__class__))
        return numlen
