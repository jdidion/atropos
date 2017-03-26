"""Interface to report generation.
"""
import importlib
import os
import platform
import sys
from atropos.io import STDOUT, STDERR
from atropos.util import MergingDict

TEXT_SERIALIZERS = ['json', 'yaml']
BINARY_SERIALIZERS = ['pickle']

def create_summary(version, command, args, options):
    """Create the initial summary dict.
    
    Args:
        version: Program version.
        command: Atropos subcommand.
        args: Other Atropos command line arguments.
        options: Parsed command line options.
    
    Returns:
        The summary dict.
    """
    summary = MergingDict()
    summary['program'] = 'Atropos'
    summary['version'] = version
    summary['python'] = platform.python_version()
    summary['command_line'] = [command] + list(args)
    summary['options'] = options.__dict__.copy()
    summary['sample_id'] = options.sample_id
    return summary

def serialize(obj, fmt, mode, outfile, **kwargs):
    """Serialize a summary dict to a file.
    
    Args:
        obj: The summary dict.
        fmt: The serialization format (e.g. json, yaml).
        mode: The file mode (b=binary, t=text).
        outfile: The output file.
        kwargs: Additional arguments to pass to the `dump` method.
    """
    mod = importlib.import_module(fmt)
    with open(outfile, 'w' + mode) as stream:
        mod.dump(obj, stream, **kwargs)

def generate_reports(summary, report_file, report_formats=None, **kwargs):
    """Generate report(s) from a summary.
    
    Args:
        summary: The summary dict.
        report_file: File name (if generating a single report) or file name
            prefix.
        report_formats: Sequence of formats to generate.
        kwargs: Additional arguments to pass to report generator function(s).
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
            report_files = (
                '{}.{}'.format(report_file, fmt)
                for fmt in report_formats)
    
    for fmt, outfile in zip(report_formats, report_files):
        if fmt in BINARY_SERIALIZERS:
            serialize(summary, fmt, 'b', outfile, **kwargs)
        elif fmt in TEXT_SERIALIZERS:
            # need to simplify some aspects of the summary to make it
            # compatible with JSON serialization
            serialize(simplify(summary), fmt, 't', outfile, **kwargs)
        else:
            #try:
            mod = importlib.import_module("atropos.reports.{}".format(fmt))
            mod.generate_report(summary, outfile, **kwargs)
            #except:
            #    import atropos.reports.jinja
            #    atropos.reports.jinja.generate_report(
            #        fmt, summary, outfile, **kwargs)

def simplify(summary):
    """Simplify a summary dict for serialization to an external format (e.g.
    JSON, YAML).
    
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
