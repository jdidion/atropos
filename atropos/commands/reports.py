"""Interface to report generation.
"""
import importlib
import os
from atropos.io import STDOUT, STDERR, open_output
from atropos.io.seqio import PAIRED

SERIALIZERS = dict(
    json='t',
    yaml='t',
    pickle='b'
)

class BaseReportGenerator(object):
    """Base class for command report generators.
    
    Args:
        options: Command-line options.
    """
    def __init__(self, options):
        report_file = options.report_file
        report_formats = options.report_formats
        if report_file in (STDOUT, STDERR):
            self.report_formats = report_formats or ('txt',)
            self.report_files = (report_file,) * len(self.report_formats)
        else:
            file_parts = os.path.splitext(report_file)
            self.report_formats = (
                report_formats or
                (file_parts[1][1:] if file_parts[1] else 'txt',))
            if len(self.report_formats) == 1:
                self.report_files = (report_file,)
            else:
                self.report_files = (
                    '{}.{}'.format(report_file, fmt)
                    for fmt in self.report_formats)
        self.report_args = tuple(
            self.get_report_args(fmt, options)
            for fmt in self.report_formats)
    
    def get_report_args(self, fmt, options):
        """Returns report-specific options dict.
        """
        return {}
    
    def generate_reports(self, summary):
        """Generate report(s) from a summary.
        
        Args:
            summary: The summary dict.
            report_file: File name (if generating a single report) or file name
                prefix.
            report_formats: Sequence of formats to generate.
            kwargs: Additional arguments to pass to report generator
                function(s).
        """
        self.add_derived_data(summary)
        
        for fmt, outfile, kwargs in zip(
                self.report_formats, self.report_files, self.report_args):
            if fmt in SERIALIZERS:
                mode = SERIALIZERS[fmt]
                self.serialize(summary, fmt, mode, outfile, **kwargs)
            else:
                self.generate_text_report(fmt, summary, outfile, **kwargs)
    
    def add_derived_data(self, summary):
        """Get some additional fields that are useful in the report.
        """
        derived = {}
        derived['mean_sequence_lengths'] = tuple(
            None if bp is None else bp / summary['total_record_count']
            for bp in summary['total_bp_counts'])
        
        inp = summary['input']
        fmt = inp['file_format']
        if inp['input_read'] == PAIRED:
            fmt += ', Paired'
        else:
            fmt += ', Read {}'.format(inp['input_read'])
        if inp['colorspace']:
            fmt += ', Colorspace'
        if inp['interleaved']:
            fmt += ', Interleaved'
        if inp['delivers_qualities']:
            fmt += ', w/ Qualities'
        else:
            fmt += ', w/o Qualities'
        derived['input_format'] = fmt
        
        summary['derived'] = derived
    
    def serialize(self, obj, fmt, mode, outfile, **kwargs):
        """Serialize a summary dict to a file.
        
        Args:
            obj: The summary dict.
            fmt: The serialization format (e.g. json, yaml).
            mode: The file mode (b=binary, t=text).
            outfile: The output file.
            kwargs: Additional arguments to pass to the `dump` method.
        """
        mod = importlib.import_module(fmt)
        with open_output(outfile, 'w' + mode) as stream:
            mod.dump(obj, stream, **kwargs)
    
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        """Generate a text report. By default, calls generate_from_template.
        """
        self.generate_from_template(fmt, summary, outfile, **kwargs)
    
    def generate_from_template(
            self, fmt, summary, outfile, template_name=None,
            template_paths=None, template_globals=None):
        """Generate a report using Jinja2.
        
        Args:
            fmt: Report format. Must match the file extension on a discoverable
                template.
            summary: The summary dict.
            outfile: The output file name/prefix.
            template_name: A template name to use, rather than auto-discover.
            template_paths: Sequence of paths to search for templates.
            template_globals: Dict of additional globals to add to the template
                environment.
        """
        import jinja2
        
        if not template_name:
            template_name = 'template.{}'.format(fmt)
        if not template_paths:
            template_paths = []
        if hasattr(self, 'template_path'):
            template_paths.append(self.template_path)
        
        # Load the report template
        try:
            env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(template_paths))
            if template_globals:
                env.globals.update(template_globals)
            template = env.get_template(template_name)
        except:
            raise IOError(
                "Could not load template file '{}'".format(template_name))

        # Render the template
        report_output = template.render(summary=summary)
        
        # Write to output
        is_path = isinstance(outfile, str)
        if is_path:
            outfile = open(outfile, 'wt', encoding='utf-8')
        else:
            report_output = report_output.encode('utf-8')
        
        try:
            print(report_output, file=outfile)
        except IOError as err:
            raise IOError(
                "Could not print report to '{}' - {}".format(outfile, err))
        finally:
            if is_path:
                outfile.close()

def prettyprint_summary(summary, outfile='summary.dump.txt'):
    """Pretty-print the summary to a file. Mostly used for debugging.
    """
    from pprint import pprint
    with open(outfile, 'w') as out:
        pprint(summary, out)
