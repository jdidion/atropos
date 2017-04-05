"""Command-line interface for the error command.
"""
from atropos.commands.cli import BaseCommandParser, writeable_file
from atropos.io import STDOUT

class CommandParser(BaseCommandParser):
    usage = """
atropos error-se input.fastq
atropos error -pe in1.fq -pe2 in2.fq
"""
    description = """
Estimate the sequencing error rate. This can help to determine
the quality of your data, and to decide the value for the max
error rate (-e) parameter. Normal error for an Illumina experiment
is around 1% or less. We recommend setting -e to 10x the empirical
error rate."""
    
    def add_command_options(self):
        parser = self.parser
        parser.set_defaults(
            max_reads=10000,
            counter_magnitude="K",
            report_formats=['txt'])
        parser.add_argument(
            "-a",
            "--algorithm",
            choices=('quality', 'shadow'), default="quality",
            help="Method for estimating error rates; quality = base qualities, "
                 "shadow = shadow regression. Be advised that the 'shadow' "
                 "method is incredibly slow.")
        parser.add_argument(
            "-m",
            "--max-bases",
            type=int, default=None,
            help="Maximum number of bases to use in the error calculation, "
                 "starting from the 5' end of the read.")
        parser.add_argument(
            "-o",
            "--output",
            type=writeable_file, default=STDOUT,
            help="File in which to write the summary of the estimated error "
                 "rates. (stdout)")
    
    def validate_command_options(self, options):
        options.report_file = options.output
