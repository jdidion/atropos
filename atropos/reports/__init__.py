import importlib
import os
import sys
from atropos.io import STDOUT, STDERR

TEXT_SERIALIZERS = ['json', 'yaml']
BINARY_SERIALIZERS = ['pickle']

def serialize(obj, fmt, mode, outfile, **kwargs):
    mod = importlib.import_module(fmt)
    with open(outfile, 'w' + mode) as o:
        mod.dump(obj, o, **kwargs)

def generate_reports(summary, report_file, report_formats=None, **kwargs):
    """Generate report(s) from a summary.
    
    Args:
        report_file: File name (if generating a single )
    """
    if report_file in (STDOUT, STDERR):
        report_files = (sys.stderr if report_file == STDERR else sys.stdout,)
        if not report_formats:
            report_formats = ('txt',)
        elif len(report_formats) > 1:
            report_files = report_file * len(report_formats)
    else:
        file_parts = os.path.splitext(report_file)
        if not report_formats:
            report_formats = (file_parts[1][1:] if file_parts[1] else 'txt',)
        if len(report_formats) == 1:
            report_files = (report_file,)
        else:
            report_files = ('{}.{}'.format(report_file, fmt) for fmt in report_formats)
    
    for fmt, outfile in zip(report_formats, report_files):
        if fmt in BINARY_SERIALIZERS:
            serialize(summary, fmt, 'b', outfile, **kwargs)
        elif fmt in TEXT_SERIALIZERS:
            # need to simplify some aspects of the summary to make it
            # compatible with JSON serialization
            serialize(simplify(summary), fmt, 't', outfile, **kwargs)
        else:
            try:
                mod = importlib.import_module("atropos.reports.{}".format(fmt))
                mod.generate_report(summary, outfile, **kwargs)
            except:
                import atropos.reports.jinja
                atropos.reports.jinja.generate_report(fmt, summary, outfile, **kwargs)

def simplify(summary):
    """
    1. JSON does not allow non-string map keys.
    """
    from collections import OrderedDict
    
    def _recurse(dest, src):
        for key, value in src.items():
            key = str(key)
            if isinstance(value, dict):
                dest[key] = OrderedDict()
                _recurse(dest[key], value)
            else:
                dest[key] = value
    
    simplified = OrderedDict()
    _recurse(simplified, summary)
    return simplified
