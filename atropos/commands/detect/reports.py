"""Report generator for the detect command.
"""
from itertools import repeat
from atropos.commands.reports import BaseReportGenerator
from atropos.commands.legacy_report import Printer, TitlePrinter
from atropos.io import open_output
from atropos.io.seqio import FastaFormat

class ReportGenerator(BaseReportGenerator):
    def get_report_args(self, fmt, options):
        if fmt == 'fasta':
            if options.fasta:
                return dict((opt, True) for opt in options.fasta)
            else:
                return dict(perinput=True)
        return {}
    
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == 'txt':
            with open_output(outfile, context_wrapper=True) as out:
                generate_reports(out, summary, **kwargs)
        elif fmt == 'fasta':
            generate_fasta(outfile, summary, **kwargs)
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
            ("{:>" + str(pad_size) + "}. Longest kmer: {}").format(
                idx+1, match['longest_kmer']))
        if match['longest_match']:
            _print_indent("Longest matching sequence: {}".format(
                match['longest_match']))
        if match['is_known']:
            _print_indent("Name(s): {}".format(
                ",\n{}".format(' ' * (pad_size + 11)).join(match['known_names'])))
            _print_indent("Known sequence(s): {}".format(
                ",\n{}".format(' ' * (pad_size + 11)).join(match['known_seqs'])))
            _print_indent(
                "Known sequence K-mers that match detected contaminant: "
                "{:.2%}".format(match['known_to_contaminant_match_frac']))
        if match['abundance']:
            _print_indent("Abundance (full-length) in {} reads: {} ({:.1%})".format(
                n_reads, match['abundance'], match['abundance'] / n_reads))
        if match['contaminant_to_known_match_frac']:
            _print_indent(
                "Detected contaminant kmers that match known sequence: "
                "{:.2%}".format(match['contaminant_to_known_match_frac']))
        if match['kmer_freq_type'] == 'frequency':
            _print_indent("Frequency of k-mers: {:.2%}".format(match['kmer_freq']))
        else:
            _print_indent("Number of k-mer matches: {}".format(match['kmer_freq']))

def generate_fasta(outfile, summary, union=False, perinput=False):
    names = summary['input']['input_names'] or repeat(None)
    n_reads = summary['record_counts'][0]
    fasta_format = FastaFormat()
    if union:
        union_records = []
    if perinput:
        if outfile.endswith('.fasta'):
            name_prefix = outfile[:-6]
        elif outfile.endswith('.fa'):
            name_prefix = outfile[:-3]
        else:
            name_prefix = outfile
    
    def format_match(idx, match, records):
        name2 = [
            "kmer_freq={}".format(match['kmer_freq']),
            "kmer_freq_type={}".format(match["kmer_freq_type"])
        ]
        if match['abundance']:
            name2.append("abundance={}".format(match['abundance']))
            name2.append("abundance_frac={}".format(match['abundance'] / n_reads))
        if match['contaminant_to_known_match_frac']:
            name2.append("contaminant_to_known_match_frac={}".format(
                match["contaminant_to_known_match_frac"]))
            
        if match['is_known']:
            name = match['known_names'][0]
            name3 = []
            if len(match['known_names']) > 1:
                name3 = ["other_names={}".format('|'.join(match['known_names'][1:]))]
            if len(match['known_seqs']) > 1:
                for seq in match['known_seqs']:
                    records.append(fasta_format.format_entry(
                        "{}.{} {}".format(name, idx, ";".join(name2 + name3)),
                        seq))
            else:
                records.append(fasta_format.format_entry(
                    "{} {}".format(name, ";".join(name2 + name3)),
                    match['known_seqs'][0]))
        else:
            records.append(fasta_format.format_entry(
                "{} {}".format(idx, ";".join(name2)), 
                match['longest_kmer']))
    
    for i, (name, matches) in enumerate(zip(names, summary['detect']['matches'])):
        records = []
        for idx, match in enumerate(matches, 1):
            format_match(idx, match, records)
        if union:
            union_records.extend(records)
        if perinput:
            with open_output("{}.{}.fasta".format(name_prefix, i), 'wt') as out:
                out.write("".join(records))
    
    if union:
        with open_output(outfile, 'wt') as union_out:
            union_out.write("".join(union_records))
