"""Report generator for the error command.

TODO: move reporting functionality out of the ErrorEstimator class.
"""
from itertools import repeat
from atropos.commands.reports import BaseReportGenerator
from atropos.io import open_output
from atropos.commands.legacy_report import Printer, TitlePrinter

class ReportGenerator(BaseReportGenerator):
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == 'txt':
            with open_output(outfile, context_wrapper=True) as out:
                generate_reports(out, summary)
        else:
            super().generate_from_template(fmt, summary, outfile, **kwargs)

def generate_reports(outstream, summary):
    names = summary['input']['input_names'] or repeat(None)
    estimates = summary['errorrate']['estimate']
    
    _print = Printer(outstream)
    _print_title = TitlePrinter(outstream)
    
    input_idx = 0
    for input_idx, (estimate, details, name) in enumerate(zip(
            estimates, summary['errorrate']['details'], names), 1):
        generate_estimator_report(
            outstream, input_idx, estimate, details, _print, _print_title, name)
    
    if input_idx > 1:
        _print.newline()
        _print_title("Overall", level=0)
        total_lens = summary['errorrate']['total_len']
        overall_err = (
            sum(err * total_len for err, total_len in zip(estimates, total_lens)) / 
            sum(total_lens))
        print("Error rate: {:.2%}".format(overall_err), file=outstream)

def generate_estimator_report(
        outstream, input_idx, estimate, details, _print, _print_title, 
        input_name=None):
    _print_indent = Printer(outstream, indent='  ')
    
    _print.newline()
    _print_title("Input {}".format(input_idx), level=0)
    
    if input_name:
        _print("File: {}".format(input_name))
    
    _print("Error rate: {:.2%}".format(estimate))
    if details:
        _print("Details:\n")
        per_read = details['per_read']
        per_cycle = details['per_cycle']
        
        _print_indent("StdErr: {:.2%}".format(per_read['standard error']))
        _print_indent("Per-cycle rates:")
        for cycle in per_cycle:
            _print_indent(
                "Cycle: {}, Error: {:.2%}, StdErr: {:.2%}".format(*cycle), 
                indent=2)
