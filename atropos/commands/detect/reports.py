from argparse import Namespace
from itertools import repeat
from pathlib import Path
from typing import IO, Optional, Sequence

from xphyle import open_

from atropos.commands.legacy_reports import Printer, TitlePrinter
from atropos.commands.reports import BaseReportGenerator, ReportWriter, BaseReportWriter
from atropos.io.formatters import FastaFormat


class DetectTextReportWriter(BaseReportWriter):
    def __init__(self, options: Namespace):
        super().__init__("txt", options)

    def serialize(self, summary: dict, stream: IO) -> None:
        names = summary["input"]["input_names"] or repeat(None)
        n_reads = summary["record_counts"][0]
        for input_idx, (matches, name) in enumerate(
            zip(summary["detect"]["matches"], names), 1
        ):
            self._generate_detector_report(stream, input_idx, n_reads, matches, name)

    @staticmethod
    def _generate_detector_report(
        stream: IO,
        input_idx: int,
        n_reads: int,
        matches: Sequence[dict],
        input_name: Optional[str] = None
    ):
        n_matches = len(matches)
        pad_size = len(str(n_matches))
        _print = Printer(stream)
        _print_title = TitlePrinter(stream)
        _print_indent = Printer(stream, indent=" " * (pad_size + 2))
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


class DetectFastaReportWriter(ReportWriter):
    """

    """
    def __init__(self, options: Namespace):
        super().__init__("fasta", options)
        self.union = False
        if options.fasta:
            self.perinput = False
            for arg in options.fasta:
                setattr(self, arg, True)
        else:
            self.perinput = True

    def write_report(self, summary: dict, output_file: Path):
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
                path = self._get_output_file(output_file, (f"{i}.fasta",))
                with open_(path, "wt") as out:
                    out.write("".join(records))

        if self.union:
            with open_(self._get_output_file(output_file), "wt") as stream:
                stream.write("".join(union_records))


class DetectReportGenerator(BaseReportGenerator):
    @classmethod
    def _create_report_writer(cls, fmt: str, options: Namespace) -> ReportWriter:
        if fmt == "txt":
            return DetectTextReportWriter(options)
        elif fmt == "fasta":
            return DetectFastaReportWriter(options)
        else:
            return super()._create_report_writer(fmt, options)
