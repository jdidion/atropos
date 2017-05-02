"""Report generator for the detect command.
"""
from itertools import repeat
from atropos.commands.reports import BaseReportGenerator
from atropos.commands.legacy_report import Printer, TitlePrinter
from atropos.io import open_output

class ReportGenerator(BaseReportGenerator):
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == 'txt':
            with open_output(outfile, context_wrapper=True) as out:
                generate_reports(out, summary)
        else:
            super().generate_from_template(fmt, summary, outfile, **kwargs)

def generate_reports(outstream, summary):
    """Prints text reports for the results from one or a pair of detectors.
    """
    names = summary['input']['input_names'] or repeat(None)
    n_reads = summary['record_counts'][0]
    for input_idx, (matches, name) in enumerate(zip(summary['detect']['matches'], names), 1):
        generate_detector_report(outstream, input_idx, n_reads, matches, name)

def generate_detector_report(outstream, input_idx, n_reads, matches, input_name=None):
    n_matches = len(matches)
    pad_size = len(str(n_matches))
    
    _print = Printer(outstream)
    _print_title = TitlePrinter(outstream)
    _print_indent = Printer(outstream, indent=' ' * (pad_size + 2))
    
    _print.newline()
    _print_title("Input {}".format(input_idx), level=0)
    
    if input_name:
        _print("File: {}".format(input_name))
    
    _print("Detected {} adapters/contaminants:".format(n_matches))
    
    if n_matches == 0:
        _print("Try increasing --max-reads")
        return
    
    for idx, match in enumerate(matches):
        _print(
            ("{:>" + str(pad_size) + "}. Longest kmer: {}").format(idx+1, match.seq))
        if match.longest_match:
            _print_indent("Longest matching sequence: {}".format(
                match.longest_match[0]))
        if match.is_known:
            _print_indent("Name(s): {}".format(
                ",\n{}".format(' ' * (pad_size + 11)).join(match.names)))
            _print_indent("Known sequence(s): {}".format(
                ",\n{}".format(' ' * (pad_size + 11)).join(match.known_seqs)))
            _print_indent(
                "Known sequence K-mers that match detected contaminant: "
                "{:.2%}".format(match.match_frac))
        if match.abundance:
            _print_indent("Abundance (full-length) in {} reads: {} ({:.1%})".format(
                n_reads, match.abundance, match.abundance / n_reads))
        if match.match_frac2:
            _print_indent(
                "Detected contaminant kmers that match known sequence: "
                "{:.2%}".format(match.match_frac2))
        if match.count_is_frequency:
            _print_indent("Frequency of k-mers: {:.2%}".format(match.count))
        else:
            _print_indent("Number of k-mer matches: {}".format(match.count))
