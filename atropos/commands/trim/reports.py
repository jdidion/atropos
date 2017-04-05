"""Report generator for trim command.
"""
import os
from atropos.commands.legacy_report import LegacyReportGenerator

class ReportGenerator(LegacyReportGenerator):
    template_path = os.path.join(os.path.dirname(__file__), 'templates')
