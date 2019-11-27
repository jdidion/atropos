"""Report generator for the detect command.
"""
from itertools import repeat
from pathlib import Path
from typing import IO

from xphyle import open_

from atropos.commands.reports import (
    BaseReportGenerator, ReportWriter, ReportWriterFactory
)
from atropos.commands.legacy_report import Printer, TitlePrinter
from atropos.io.seqio import FastaFormat


class FastaReportWriter(ReportWriter):
    """

    """
    def __init__(
        self,
        output_file: Path,
        options=None,
        **_
    ):
        super().__init__(output_file)
        self.union = False
        if options.fasta:
            self.perinput = False
            for arg in options.fasta:
                setattr(self, arg, True)
        else:
            self.perinput = True

    def serialize(self, summary: dict, stream: IO):
        names = summary["input"]["input_names"] or repeat(None)
        n_reads = summary["record_counts"][0]
        fasta_format = FastaFormat()
        union_records = []

        def format_match(_idx, _match, _records):
            name2 = [
                f"kmer_freq={_match['kmer_freq']}",
                f"kmer_freq_type={_match['kmer_freq_type']}",
            ]
            if _match["abundance"]:
                name2.append(f"abundance={_match['abundance']}")
                name2.append(f"abundance_frac={_match['abundance'] / n_reads}")
            if _match["contaminant_to_known_match_frac"]:
                name2.append(
                    f"contaminant_to_known_match_frac="
                    f"{_match['contaminant_to_known_match_frac']}"
                )
            if _match["is_known"]:
                name1 = _match["known_names"][0]
                name3 = []
                if len(_match["known_names"]) > 1:
                    name3 = [
                        f"other_names={'|'.join(_match['known_names'][1:])}"
                    ]
                if len(_match["known_seqs"]) > 1:
                    for seq in _match["known_seqs"]:
                        _records.append(
                            fasta_format.format_entry(
                                f"{name1}.{_idx} {';'.join(name2 + name3)}", seq
                            )
                        )
                else:
                    _records.append(
                        fasta_format.format_entry(
                            f"{name1} {';'.join(name2 + name3)}",
                            _match["known_seqs"][0],
                        )
                    )
            else:
                _records.append(
                    fasta_format.format_entry(
                        f"{_idx} {';'.join(name2)}", _match["longest_kmer"]
                    )
                )

        for i, (name, matches) in enumerate(zip(names, summary["detect"]["matches"])):
            records = []
            for idx, match in enumerate(matches, 1):
                format_match(idx, match, records)
            if self.union:
                union_records.extend(records)
            if self.perinput:
                path = self.output_file.with_suffix(f".{i}.fasta")
                with open_(path, "wt") as out:
                    out.write("".join(records))

        if self.union:
            stream.write("".join(union_records))


class ReportGenerator(BaseReportGenerator):
    def __init__(self, ):
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == "txt":
            with open_(outfile, "wt", context_wrapper=True) as out:
                generate_reports(out, summary, **kwargs)
        elif fmt == "fasta":
            generate_fasta(outfile, summary, **kwargs)
        else:
            super().generate_from_template(fmt, summary, outfile, **kwargs)


def generate_reports(outstream, summary):
    """Prints text reports for the results from one or a pair of detectors.
    """
    names = summary["input"]["input_names"] or repeat(None)
    n_reads = summary["record_counts"][0]
    for input_idx, (matches, name) in enumerate(
        zip(summary["detect"]["matches"], names), 1
    ):
        generate_detector_report(outstream, input_idx, n_reads, matches, name)


def generate_detector_report(outstream, input_idx, n_reads, matches, input_name=None):
    n_matches = len(matches)
    pad_size = len(str(n_matches))
    _print = Printer(outstream)
    _print_title = TitlePrinter(outstream)
    _print_indent = Printer(outstream, indent=" " * (pad_size + 2))
    _print.newline()
    _print_title(f"Input {input_idx}", level=0)
    if input_name:
        _print(f"File: {input_name}")
    _print(f"Detected {n_matches} adapters/contaminants:")
    if n_matches == 0:
        _print("Try increasing --max-reads")
        return

    for idx, match in enumerate(matches):
        _print(
            ("{:>" + str(pad_size) + "}. Longest kmer: {}").format(
                idx + 1, match["longest_kmer"]
            )
        )
        if match["longest_match"]:
            _print_indent(
                f"Longest matching sequence: {match['longest_match']}"
            )
        if match["is_known"]:
            _print_indent(
                "Name(s): {}".format(
                    f",\n{' ' * (pad_size + 11)}".join(match["known_names"])
                )
            )
            _print_indent(
                "Known sequence(s): {}".format(
                    f",\n{' ' * (pad_size + 11)}".join(match["known_seqs"])
                )
            )
            _print_indent(
                f"Known sequence K-mers that match detected contaminant: "
                f"{match['known_to_contaminant_match_frac']:.2%}"
            )
        if match["abundance"]:
            _print_indent(
                f"Abundance (full-length) in {n_reads} reads: {match['abundance']} "
                f"({match['abundance'] / n_reads:.1%})"
            )
        if match["contaminant_to_known_match_frac"]:
            _print_indent(
                f"Detected contaminant kmers that match known sequence: "
                f"{match['contaminant_to_known_match_frac']:.2%}"
            )
        if match["kmer_freq_type"] == "frequency":
            _print_indent(f"Frequency of k-mers: {match['kmer_freq']:.2%}")
        else:
            _print_indent(f"Number of k-mer matches: {match['kmer_freq']}")
