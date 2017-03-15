import importlib
import os

def generate_reports(options, summary):
    if not options.report_file:
        return None
    report_formats = options.report_formats
    if not report_formats:
        fmt = os.path.splitext(options.report_file)[1]
        report_formats = (fmt[1:] if fmt else 'txt',)
    for fmt, outfile in zip(report_formats, options.report_file):
        #try:
        mod = importlib.import_module("atropos.reports.{}".format(report_type))
        mod.generate_report(summary, outfile)
        #except:
        #    from atropos.reports.jinja import generate_report
        #    generate_report(fmt, summary, outfile)
