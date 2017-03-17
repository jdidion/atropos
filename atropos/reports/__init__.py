import importlib
import os
import sys

def generate_reports(summary, report_file, report_formats=None):
    """Generate report(s) from a summary.
    
    Args:
        report_file: File name (if generating a single )
    """
    if report_file in ('-', 'stdout', 'stderr'):
        report_file = (sys.stderr if report_file == 'stderr' else sys.stdout,)
        if not report_formats:
            report_formats = ('txt',)
        elif len(report_formats) > 1:
            report_file = report_file * len(report_formats)
    else:
        file_parts = os.path.splitext(report_file)
        if not report_formats:
            report_formats = (file_parts[1][1:] if file_parts[1] else 'txt',)
        if len(report_formats) == 1:
            report_files = (report_file,)
        else:
            report_files = ('{}.{}'.format(report_file, fmt) for fmt in report_formats)
    
    for fmt, outfile in zip(report_formats, report_files):
        #try:
        mod = importlib.import_module("atropos.reports.{}".format(fmt))
        mod.generate_report(summary, outfile)
        #except:
        #    from atropos.reports.jinja import generate_report
        #    generate_report(fmt, summary, outfile)
