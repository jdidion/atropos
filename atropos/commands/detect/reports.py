"""Report generator for the detect command.

TODO: move reporting functionality out of the Detector class.
"""
from atropos.commands.reports import BaseReportGenerator
from atropos.io import open_output

class ReportGenerator(BaseReportGenerator):
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == 'txt':
            detector = summary['detect']['detector']
            with open_output(outfile) as out:
                detector.summarize(
                    out, summary['input']['input_names'],
                    summary['detect']['include'])
        else:
            super().generate_from_template(fmt, summary, outfile, **kwargs)
