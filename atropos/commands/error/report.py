"""Report generator for the error command.

TODO: move reporting functionality out of the ErrorEstimator class.
"""
from atropos.commands.reports import BaseReportGenerator
from atropos.io import open_output

class ReportGenerator(BaseReportGenerator):
    def generate_text_report(self, fmt, summary, outfile, **kwargs):
        if fmt == 'txt':
            estimator = summary['error_estimator']
            with open_output(outfile) as out:
                estimator.summarize(out, summary['input']['input_names'])
        else:
            super().generate_from_template(fmt, summary, outfile, **kwargs)
