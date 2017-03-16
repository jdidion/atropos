import importlib
import os

def generate_reports(summary, report_file, report_formats=None, quiet=False):
    """Generate report(s) from a summary.
    """
    if not report_file:
        return
    
    file_parts = os.path.splitext(report_file)
    
    if not report_formats:
        fmt = file_parts[1]
        report_formats = (fmt[1:] if fmt else 'txt',)
    
    if len(report_formats) == 1:
        report_files = (report_file,)
    else:
        report_files = ('{}.{}'.format(file_parts[0], fmt) for fmt in report_formats)
    
    for fmt, outfile in zip(report_formats, report_files):
        #try:
        mod = importlib.import_module("atropos.reports.{}".format(fmt))
        mod.generate_report(summary, outfile, quiet)
        #except:
        #    from atropos.reports.jinja import generate_report
        #    generate_report(fmt, summary, outfile)
