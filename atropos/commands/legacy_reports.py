from argparse import Namespace

import math
import textwrap
from typing import Any, Dict, IO, Optional, Sequence, Tuple, TypeVar, Union, cast

from atropos.commands.reports import BaseReportGenerator, BaseReportWriter, ReportWriter
from atropos.utils import truncate_string
from atropos.utils.statistics import weighted_median


INDENT = " " * 2
PARAGRAPH = textwrap.TextWrapper()
INDENTED = textwrap.TextWrapper(initial_indent=INDENT, subsequent_indent=INDENT)


class Printer:
    """
    Manages printing to a file.

    Args:
        outfile: The output file.
        kwargs: Additional keyword arguments passed to the print function.
    """

    def __init__(self, outfile: IO, indent: str = "", **kwargs):
        self.outfile = outfile
        self.indent = indent
        self.print_args = kwargs

    def __call__(self, *args, indent: Union[int, str] = None, **kwargs):
        if isinstance(indent, int):
            indent_str = self.indent * cast(int, indent)
        else:
            indent_str = cast(str, indent) or self.indent

        if indent_str:
            self._print(indent_str, end="")

        self._print(*args, **kwargs)

    def _print(self, *args, **kwargs):
        if self.print_args:
            print_args = self.print_args.copy()
            print_args.update(kwargs)
        else:
            print_args = kwargs

        print(*args, file=self.outfile, **print_args)

    def newline(self):
        """
        Prints a newline.
        """
        print(file=self.outfile)


class TitlePrinter(Printer):
    """
    Printer that formats titles.

    Args:
        outfile: The output file.
        levels: The formatting associated with different header levels.
        kwargs: Additional keyword arguments passed to the print function.
    """

    def __init__(
        self,
        outfile: IO,
        levels: Sequence[Tuple[str, Optional[str]]] = (
            ("=", "="),
            ("-", "-"),
            ("-", None),
            ("~", None),
        ),
        **kwargs,
    ):
        super().__init__(outfile, **kwargs)
        self.levels = levels

    def __call__(
        self, *title: str, level: Optional[int] = None, newline: bool = True, **kwargs
    ):
        title = " ".join(title)
        underline = None
        width = 0

        if level is not None:
            if level >= len(self.levels):
                raise ValueError(f"Invalid level: {level}")

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
    """
    Printer that formats rows in a table.
    """

    def __init__(
        self,
        outfile: IO,
        colwidths: Union[int, Tuple[int, ...]] = 10,
        justification: Union[str, Tuple[str, ...]] = ("<", ">"),
        indent: Union[str, Tuple[str, ...]] = "",
        pct: bool = False,
        default: int = 0,
        **kwargs,
    ):
        """
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
        super().__init__(outfile, **kwargs)
        self.colwidths, self.justification, self.indent = (
            (arg,) if isinstance(arg, typ) else tuple(arg)
            for arg, typ in zip((colwidths, justification, indent), (int, str, str))
        )
        self.pct = pct
        self.default = default

    def print_rows(
        self,
        *rows,
        header: Optional[Sequence[Union[str, Tuple[str, str]]]] = None,
        **kwargs,
    ):
        """
        Prints multiple rows. Automatically figures out column widths.

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
                    max(sizeof(h) for h in header_part) for header_part in header
                )
                header_rows = list(zip(*header))

            colwidths = tuple(max(h, c) for h, c in zip(header_widths, colwidths))

            for i, header_row in enumerate(header_rows, 1):
                self(
                    *header_row,
                    colwidths=colwidths,
                    header=(i == len(header_rows)),
                    **kwargs,
                )

        for row in rows:
            self(*row, colwidths=colwidths)

    def __call__(
        self,
        *args,
        colwidths: Optional[Tuple[int, ...]] = None,
        extra_width: Optional[int] = None,
        justification: Optional[Tuple[str, ...]] = None,
        extra_justification: Optional[str] = None,
        indent: Optional[Tuple[str, ...]] = None,
        extra_indent: Optional[str] = None,
        header: bool = False,
        underline: str = "-",
        pct: Optional[bool] = None,
        default: Optional[Any] = None,
        **kwargs,
    ):
        """
        Prints a row.

        Args:
            args: Fields in the row.
            colwidths: Row-specific colwidths.
            justification: Row-specific justifications.
            indent: Row-specific indents.
            extra_width: colwidth to use for extra fields.
            extra_justification: Justification to use for extra fields.
            extra_indent: Indent to use for extra fields.
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

        T = TypeVar("T")

        def adjust(tup: Tuple[T, ...], extra: Optional[T] = None) -> Tuple[T, ...]:
            """
            Adjust a tuple. If longer than the number of columns, truncate; if
            shorter, fill in by repeating the last element.
            """
            tlen = len(tup)
            if tlen == ncols:
                return tup
            elif tlen > ncols:
                return tup[:ncols]
            else:
                return tup + ((extra or tup[-1],) * (ncols - tlen))

        colwidths, justification, indent = (
            adjust(arr, extra)
            for arr, extra in zip(
                (
                    colwidths or self.colwidths,
                    justification or self.justification,
                    indent or self.indent,
                ),
                (extra_width, extra_justification, extra_indent),
            )
        )

        # adjust colwidths if this is a header
        if header:
            colwidths = tuple(max(w, len(str(a))) for w, a in zip(colwidths, args))

        fmt_str = []
        fmt_args = []

        for i, (value, width, just, ind) in enumerate(
            zip(args, colwidths, justification, indent)
        ):
            if value is None:
                value = default or self.default

            if isinstance(value, str):
                typ = "s"
                if len(value) > width:
                    value = truncate_string(value, width)
            elif isinstance(value, float):
                typ = ",.1" + ("%" if pct else "f")
            else:
                typ = ",d"

            fmt_str.append(
                ind + "{" + str(i) + ":" + just + str(width - len(ind)) + typ + "}"
            )
            fmt_args.append(value)

        self._print(" ".join(fmt_str).format(*fmt_args), **kwargs)

        if header:
            sepline = " ".join((underline * width) for width in colwidths)
            self._print(sepline, **kwargs)


class LegacyTextReportWriter(BaseReportWriter):
    def __init__(self, options: Namespace):
        super().__init__("txt", options)

    def serialize(self, summary: dict, stream: IO) -> None:
        self._print_summary_report(summary, stream)
        if "trim" in summary:
            self._print_trim_report(summary, stream)
        if "pre" in summary:
            self._print_pre_trim_report(summary, stream)
        if "post" in summary:
            self._print_post_trim_report(summary, stream)

    @classmethod
    def _print_summary_report(cls, summary: dict, stream: IO):
        """
        Prints the top-level summary report.

        Args:
            summary: The summary dict.
            stream: The output file object.
        """
        _print_title = TitlePrinter(stream)
        _print = Printer(stream)
        _print_title("Atropos", level=0)
        _print(f"Atropos version: {summary['version']}")
        _print(f"Python version: {summary['python']}")
        _print(
            f"Command line parameters: {summary['command']} "
            f"{' '.join(summary['options']['orig_args'])}"
        )
        _print()
        _print(f"Sample ID: {summary['sample_id']}")
        _print(f"Input format: {summary['derived']['input_format']}")
        _print("Input files:")

        for infile in summary["input"]["input_names"]:
            if infile is not None:
                _print(infile, indent=INDENT)

        _print()

        timing = summary["timing"]
        total = summary["total_record_count"]
        wctime = [f"Wallclock time: {timing['wallclock']:.2F} s"]

        if total > 0:
            wctime.append(
                f"({1e6 * timing['wallclock'] / total:.0F} us/read; "
                f"{total / timing['wallclock'] * 60 / 1e6:.2F} M reads/minute)"
            )

        _print(f"Start time: {timing['start']}")
        _print(*wctime)
        _print(f"CPU time (main process): {timing['cpu']:.2F} s")
        _print()

    @classmethod
    def _print_trim_report(cls, summary: dict, stream: IO):
        """
        Prints the trimming report.

        Args:
            summary: Summary dict.
            stream: Open output stream.
        """
        total_bp = sum(summary["total_bp_counts"])
        max_width = len(str(total_bp))
        # account for commas
        max_width += max_width // 3
        _print = RowPrinter(stream, (35, max_width))
        total = summary["total_record_count"]

        if total == 0:
            _print(
                "No reads processed! Either your input file is empty or you "
                "used the wrong -f/--format parameter."
            )
            return

        modifiers, filters, formatters = (
            summary["trim"][key] for key in ("modifiers", "filters", "formatters")
        )

        adapter_cutter = None
        error_corrector = None

        for modifier_dict in modifiers.values():
            if adapter_cutter is None and "adapters" in modifier_dict:
                adapter_cutter = modifier_dict
                break

            if error_corrector is None and "bp_corrected" in modifier_dict:
                error_corrector = modifier_dict

        correction_enabled = summary["options"]["correct_mismatches"]
        corrected = None
        trimmers = []

        for name, mod in modifiers.items():
            if "bp_trimmed" in mod:
                trimmers.append((name, mod))

            if correction_enabled and "records_corrected" in mod:
                corrected = mod

        _print_title = TitlePrinter(stream)
        paired = summary["options"]["paired"]
        pairs_or_reads = "Pairs" if paired else "Reads"

        _print_title("Trimming", level=1)
        _print(pairs_or_reads, "records", "fraction", header=True)
        _print(f"Total {'read pairs' if paired else 'reads'} processed:", total)

        if adapter_cutter:
            if paired:
                for read in range(2):
                    _print(
                        f"Read {read + 1} with adapter:",
                        adapter_cutter["records_with_adapters"][read],
                        adapter_cutter["fraction_records_with_adapters"][read],
                        indent=(INDENT, ""),
                        pct=True,
                    )
            else:
                _print(
                    "Reads with adapters:",
                    adapter_cutter["records_with_adapters"][0],
                    adapter_cutter["fraction_records_with_adapters"][0],
                    pct=True,
                )

        def _print_filter(_name, sep):
            if _name in filters:
                _print(
                    f"{pairs_or_reads} {sep} {_name.replace('_', ' ')}:",
                    filters[_name]["records_filtered"],
                    filters[_name]["fraction_records_filtered"],
                    pct=True,
                )

        _print_filter("too_short", "that were")
        _print_filter("too_long", "that were")
        _print_filter("too_many_n", "with")
        _print(
            f"{pairs_or_reads} written (passing filters):",
            formatters["records_written"],
            formatters["fraction_records_written"],
            pct=True,
        )

        if corrected:
            _print(
                "Pairs corrected:",
                corrected["records_corrected"],
                corrected["fraction_records_corrected"],
                pct=True,
            )

        _print()
        _print("Base pairs", "bp", "fraction", header=True)
        _print("Total bp processed:", total_bp)

        if paired:
            for read in range(2):
                _print(
                    f"Read {read + 1}:",
                    summary["total_bp_counts"][read],
                    indent=(INDENT, ""),
                )

        def _print_bp(title, data, key, default=0):
            if paired:
                _print(
                    title,
                    data[f"total_{key}"],
                    data[f"fraction_total_{key}"],
                    pct=True,
                )

                for _read in range(2):
                    _print(
                        f"Read {_read + 1}:",
                        data[key][_read],
                        data[f"fraction_{key}"][_read],
                        indent=(INDENT, ""),
                        pct=True,
                        default=default,
                    )
            else:
                _print(
                    title,
                    data[key][0],
                    data[f"fraction_{key}"][0],
                    pct=True,
                    default=default,
                )

        for name, mod in trimmers:
            _print_bp(mod["desc"], mod, "bp_trimmed")

        _print_bp("Total bp written (filtered):", formatters, "bp_written")

        if error_corrector:
            _print_bp("Total bp corrected:", error_corrector, "bp_corrected")

        if adapter_cutter:
            _print()
            adapters = adapter_cutter["adapters"]
            cls._print_adapter_report(adapters, stream, paired, total, max_width)

    @classmethod
    def _print_adapter_report(
        cls,
        adapters: Sequence[dict],
        stream: IO,
        paired: bool,
        total_records: int,
        max_width: int,
    ):
        """
        Prints details for a adapters.

        Args:
            adapters: Sequence of adapter info dicts.
            stream: The output file
            paired: Whether the data is paired-end.

            max_width: Max column width.
        """
        adapter_lenghts = []

        for pair in adapters:
            if pair:
                for adapter in pair.values():
                    if adapter["where"]["name"] == "linked":
                        adapter_lenghts.append(
                            3
                            + len(adapter["front_sequence"] + adapter["back_sequence"])
                        )
                    else:
                        adapter_lenghts.append(len(adapter["sequence"]))

        max_seq_len = max(adapter_lenghts)
        _print = Printer(stream)
        _print_title = TitlePrinter(stream)
        _print_adj = RowPrinter(stream, (12, 5), pct=True, indent=(INDENT, ""))
        seq_printer = RowPrinter(
            stream, (max_seq_len, 14, 3, max_width), ("<", "<", ">")
        )
        hist_printer = RowPrinter(stream, justification=(">", ">", ">", ">", "<"))

        def print_error_ranges(adapter_length: int, error_rate: float):
            """
            Prints max number of errors for ranges of adapter match lengths.
            """
            _print("No. of allowed errors:")

            prev = 0

            for errors in range(1, int(error_rate * adapter_length) + 1):
                range_start = int(errors / error_rate)
                _print(f"{prev}-{range_start - 1} bp: {errors - 1};", end=" ")
                prev = range_start

            if prev == adapter_length:
                _print(f"{adapter_length} bp: {int(error_rate * adapter_length)}")
            else:
                _print("{prev}-{adapter_length} bp: {int(error_rate * adapter_length)}")

            _print()

        def print_histogram(
            data: dict,
            adapter_length: int,
            num_reads: int,
            error_rate: float,
            errors: dict,
            match_probabilities: Sequence[float],
        ):
            """
            Prints a histogram. Also, print the no. of reads expected to be
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
                estimated = num_reads * match_probabilities[min(length, adapter_length)]
                hist.append(
                    [
                        length,
                        count,
                        estimated,
                        int(error_rate * min(length, adapter_length)),
                    ]
                )
                hist_errors.append(errors["rows"][length])

            col_sizes = [len(str(max(col))) for col in zip(*hist_errors)]

            def _format_hist_errors(_errs):
                trailing = True
                hist_str = []

                for e_ctr, e in reversed(list(enumerate(_errs))):
                    if trailing and e == 0:
                        pass
                    else:
                        trailing = False
                        hist_str.append(
                            ("{:<" + str(col_sizes[e_ctr]) + "d}").format(e)
                        )

                return " ".join(reversed(hist_str))

            for i, errs in enumerate(hist_errors):
                hist[i].append(_format_hist_errors(errs))

            error_header = " ".join(
                ("{:<" + str(e) + "d}").format(i) for i, e in enumerate(col_sizes)
            )

            hist_printer.print_rows(
                *hist,
                header=(
                    ("length", ""),
                    ("count", ""),
                    ("expect", ""),
                    ("max.err", ""),
                    ("error counts", error_header),
                ),
            )
            hist_printer.newline()

        def print_adjacent_bases(bases: Dict[str, int]):
            """
            Prints a summary of the bases preceding removed adapter sequences.
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

            for base in ["A", "C", "G", "T", ""]:
                base_label = base if base != "" else "none/other"
                fraction = bases[base] / total
                _print_adj(base_label, fraction)
                if fraction > 0.8 and base != "":
                    warnbase = base_label

            if total >= 20 and warnbase is not None:
                _print("WARNING:")
                _print(
                    "\n".join(
                        INDENTED.wrap(
                            f"The adapter is preceded by '{warnbase}'extremely often. "
                            f"The provided adapter sequence may be incomplete. To fix "
                            f"the problem, add '{warnbase}' to the beginning of the "
                            f"adapter sequence."
                        )
                    )
                )
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
                front_len = 0
                back_len = 0
                seq_len = 0

                if where_name == "linked":
                    front_len, back_len = [
                        len(adapter[s]) for s in ("front_sequence", "back_sequence")
                    ]
                    seq_printer.print_rows(
                        (
                            f"{adapter['front_sequence']}...{adapter['back_sequence']}",
                            "linked",
                            f"{front_len}+{back_len}",
                            adapter["total_front"],
                            adapter["total_back"],
                        ),
                        header=(
                            "Sequence",
                            "Type",
                            "Length",
                            "Trimmed (x)",
                            "Half matches (x)",
                        ),
                    )
                else:
                    seq_len = len(adapter["sequence"])
                    seq_printer.print_rows(
                        (
                            adapter["sequence"],
                            adapter["where"]["desc"],
                            seq_len,
                            adapter["total"],
                        ),
                        header=("Sequence", "Type", "Length", "Trimmed (x)"),
                    )

                _print()

                if adapter["total"] == 0:
                    continue

                if where_name == "anywhere":
                    _print(
                        adapter["total_front"],
                        "times, it overlapped the 5' end of a read",
                    )
                    _print(
                        adapter["total_back"],
                        "times, it overlapped the 3' end or was within the read",
                    )
                    _print()
                    print_error_ranges(seq_len, adapter["max_error_rate"])
                    _print("Overview of removed sequences (5'):")
                    print_histogram(
                        adapter["lengths_front"],
                        seq_len,
                        total_records,
                        adapter["max_error_rate"],
                        adapter["errors_front"],
                        adapter["match_probabilities"],
                    )
                    _print()
                    _print("Overview of removed sequences (3' or within):")
                    print_histogram(
                        adapter["lengths_back"],
                        seq_len,
                        total_records,
                        adapter["max_error_rate"],
                        adapter["errors_back"],
                        adapter["match_probabilities"],
                    )
                elif where_name == "linked":
                    print_error_ranges(front_len, adapter["front_max_error_rate"])
                    print_error_ranges(back_len, adapter["back_max_error_rate"])
                    _print("Overview of removed sequences at 5' end:")
                    print_histogram(
                        adapter["front_lengths_front"],
                        front_len,
                        total_records,
                        adapter["front_max_error_rate"],
                        adapter["front_errors_front"],
                        adapter["front_match_probabilities"],
                    )
                    _print()
                    _print("Overview of removed sequences at 3' end:")
                    print_histogram(
                        adapter["back_lengths_back"],
                        back_len,
                        total_records,
                        adapter["back_max_error_rate"],
                        adapter["back_errors_back"],
                        adapter["back_match_probabilities"],
                    )
                elif where_name in ("front", "prefix"):
                    print_error_ranges(seq_len, adapter["max_error_rate"])
                    _print("Overview of removed sequences:")
                    print_histogram(
                        adapter["lengths_front"],
                        seq_len,
                        total_records,
                        adapter["max_error_rate"],
                        adapter["errors_front"],
                        adapter["match_probabilities"],
                    )
                elif where_name in ("back", "suffix"):
                    print_error_ranges(seq_len, adapter["max_error_rate"])
                    warning = warning or print_adjacent_bases(adapter["adjacent_bases"])
                    _print("Overview of removed sequences:")
                    print_histogram(
                        adapter["lengths_back"],
                        seq_len,
                        total_records,
                        adapter["max_error_rate"],
                        adapter["errors_back"],
                        adapter["match_probabilities"],
                    )

        if warning:
            _print("WARNING:")
            _print(
                "\n".join(
                    INDENTED.wrap(
                        "One or more of your adapter sequences may be incomplete. "
                        "Please see the detailed output above."
                    )
                )
            )

    @classmethod
    def _print_pre_trim_report(cls, summary: dict, stream: IO):
        """
        Prints pre-trimming metrics.

        Args:
            summary: The summary dict.
            stream: The output file.
        """
        pre = summary["pre"]
        _print_title = TitlePrinter(stream)
        _print = Printer(stream)
        _print_title("Pre-trimming metrics", level=1)

        for source, data in pre.items():
            _print_title("Source", level=3, newline=False)

            # TODO: When multi-file input is supported, this code will need to
            #  get summary['input']['input_names'][source]
            for read, src in enumerate(summary["input"]["input_names"], 1):
                if src is None:
                    continue

                _print(f"Read {read}: {src}")

            _print()

            cls._print_metrics_report(data, stream)

    @classmethod
    def _print_post_trim_report(cls, summary: dict, outfile: IO):
        """
        Prints post-trimming metrics.

        Args:
            summary: The summary dict.
            outfile: The output file.
        """
        post = summary["post"]
        _print_title = TitlePrinter(outfile)
        _print = Printer(outfile)
        _print_title("Post-trimming metrics", level=1)

        for dest, metrics in post.items():
            _print_title(f"Destination: {dest}", level=2)

            for source, data in metrics.items():
                _print_title("Source", level=3, newline=False)

                # TODO: When multi-file input is supported, this code will need to
                #  get summary['input']['input_names'][source]

                for read, src in enumerate(summary["input"]["input_names"], 1):
                    if src is None:
                        continue

                    _print(f"Read {read}: {src}")

                _print()

                cls._print_metrics_report(data, outfile)

    @classmethod
    def _print_metrics_report(cls, data, outfile):
        """
        Prints metrics.

        Args:
            data: The metrics dict.
            outfile: The output file.
        """
        paired = "read2" in data
        max_count = data["read1"]["counts"]

        if paired:
            max_count = max(max_count, data["read2"]["counts"])

        max_width = len(str(max_count))
        # add space for commas and column separation
        max_width += (max_width // 3) + 1
        _print_title = TitlePrinter(outfile)
        _print = RowPrinter(outfile, (35, max_width))

        def _print_histogram(title, hist1, hist2=None):
            _print_title(title, level=2)

            if hist1 is None:
                _print("No Data")
                return

            if hist2:
                hist = (
                    (key, hist1.get(key, 0), hist2.get(key, 0))
                    for key in sorted(set(hist1.keys()) | set(hist2.keys()))
                )
            else:
                hist = sorted(hist1.items(), key=lambda x: x[0])

            for histbin in hist:
                _print(*histbin)

        def _print_base_histogram(title, hist, extra_width=4, index_name="Pos"):
            _print_title(title, level=2)

            if hist is None:
                _print("No Data")
                return

            _print(index_name, *hist["columns"], header=True, extra_width=extra_width)

            for pos, row in hist["rows"].items():
                total_count = sum(row)
                base_pcts = (round(count * 100 / total_count, 1) for count in row)
                _print(pos, *base_pcts, extra_width=extra_width)

        def _print_tile_histogram(title, hist):
            if hist is None:
                _print_title(title, level=2)
                _print("No Data")
                return

            ncol = len(hist["columns"])
            max_tile_width = (
                max(4, len(str(math.ceil(data["read1"]["counts"] / ncol)))) + 1
            )
            _print_base_histogram(
                title, hist, extra_width=max_tile_width, index_name="Tile"
            )

        def _print_tile_base_histogram(title: str, hist: dict):
            """Print a histogram of position x tile, with values as the median
            base quality.
            """
            _print_title(title, level=2)
            if hist is None:
                _print("No Data")
                return

            quals = hist["columns"]
            tiles = hist["columns2"]
            ncol = len(tiles)
            max_tile_width = (
                max(4, len(str(math.ceil(data["read1"]["counts"] / ncol)))) + 1
            )
            _print("Pos", *tiles, header=True, extra_width=max_tile_width)

            for pos, tiles in hist["rows"].items():
                # compute the weighted median for each tile at each position
                _print(
                    pos,
                    *(
                        weighted_median(quals, tile_counts)
                        for tile_counts in tiles.values()
                    ),
                    extra_width=max_tile_width,
                )

        _print("", "Read1", "Read2", header=True)

        # Sequence-level metrics
        _print(
            "Read pairs:" if paired else "Reads:",
            data["read1"]["counts"],
            data["read2"]["counts"],
        )
        _print()
        _print_histogram(
            "Sequence lengths:",
            data["read1"]["lengths"]["hist"],
            data["read2"]["lengths"]["hist"],
        )
        _print()
        if "qualities" in data["read1"]:
            _print_histogram(
                "Sequence qualities:",
                data["read1"]["qualities"]["hist"],
                data["read2"]["qualities"]["hist"],
            )
            _print()
        _print_histogram(
            "Sequence GC content (%)",
            data["read1"]["gc"]["hist"],
            data["read2"]["gc"]["hist"],
        )
        _print()
        if "tile_sequence_qualities" in data["read1"]:
            _print_tile_histogram(
                "Read 1 per-tile sequence qualities (%)",
                data["read1"]["tile_sequence_qualities"],
            )
            _print()
            _print_tile_histogram(
                "Read 2 per-tile sequence qualities (%)",
                data["read2"]["tile_sequence_qualities"],
            )
            _print()

        # Base-level metrics
        if "base_qualities" in data["read1"]:
            _print_base_histogram(
                "Read 1 base qualities (%)", data["read1"]["base_qualities"]
            )
            _print()
            _print_base_histogram(
                "Read 2 base qualities (%)", data["read2"]["base_qualities"]
            )
            _print()
        _print_base_histogram("Read 1 base composition (%)", data["read1"]["bases"])
        _print()
        _print_base_histogram("Read 2 base composition (%)", data["read2"]["bases"])
        _print()

        if "tile_base_qualities" in data["read1"]:
            _print_tile_base_histogram(
                "Read 1 per-tile base qualities (%)",
                data["read1"]["tile_base_qualities"],
            )
            _print()
            _print_tile_base_histogram(
                "Read 2 per-tile base qualities (%)",
                data["read2"]["tile_base_qualities"],
            )
            _print()


def sizeof(*x: Union[str, int, float], seps: bool = True, prec: int = 1):
    """
    Returns the largest string size of all objects in x, where x is a sequence of
    string or numeric values.

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
                numlen += numlen // 3
        elif isinstance(x[0], float):
            numlen = len(str(round(max(x), prec)))
            if seps:
                numlen += (numlen - prec - 1) // 3
        else:
            raise ValueError(f"Unexpected data type: {x[0].__class__}")

        return numlen


class LegacyReportGenerator(BaseReportGenerator):
    @classmethod
    def _create_report_writer(cls, fmt: str, options: Namespace) -> ReportWriter:
        if fmt == "txt":
            return LegacyTextReportWriter(options)
        else:
            return super()._create_report_writer(fmt, options)
