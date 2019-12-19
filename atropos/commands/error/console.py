from xphyle import STDOUT

from atropos.commands.console import (
    BaseCommandConsole, add_common_options, validate_common_options
)
from atropos.commands.error import ErrorCommand
from atropos.commands.error.reports import ErrorReportGenerator
from atropos.utils import classproperty
from atropos.utils.argparse import (
    AtroposArgumentParser,
    Namespace,
    writable_file
)


class ErrorCommandConsole(ErrorCommand, ErrorReportGenerator, BaseCommandConsole):
    @classproperty
    def description(cls) -> str:
        return """
        atropos error-se input.fastq
        atropos error -pe in1.fq -pe2 in2.fq
        """

    @classproperty
    def usage(cls) -> str:
        return """
        Estimate the sequencing error rate. Use this command to help determine
        the quality of your data, and to decide the value for the max
        error rate (-e) parameter. Normal error for an Illumina experiment
        is around 1% or less. We recommend setting -e to 10x the empirical
        error rate."""

    @classmethod
    def _add_arguments(cls, parser: AtroposArgumentParser) -> None:
        add_common_options(parser)
        cls._add_error_options(parser)

    @classmethod
    def _validate_options(
        cls, options: Namespace, parser: AtroposArgumentParser
    ) -> None:
        validate_common_options(options, parser)
        cls._validate_error_options(options, parser)

    @staticmethod
    def _add_error_options(parser: AtroposArgumentParser):
        parser.set_defaults(
            max_reads=10000, counter_magnitude="K", report_formats=['txt']
        )
        group = parser.add_group("Error Estimation")
        group.add_argument(
            "-a",
            "--algorithm",
            choices=('quality', 'shadow'),
            default="quality",
            help="Method for estimating error rates; quality = base qualities, "
            "shadow = shadow regression. Be advised that the 'shadow' "
            "method is incredibly slow.",
        )
        group.add_argument(
            "-m",
            "--max-bases",
            type=int,
            default=None,
            help="Maximum number of bases to use in the error calculation, "
            "starting from the 5' end of the read.",
        )
        group = parser.add_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writable_file,
            default=STDOUT,
            help="File in which to write the summary of the estimated error "
            "rates. (stdout)",
        )
        group.add_argument(
            "--output_formats",
            nargs="*",
            choices=("txt", "json", "yaml", "pickle"),
            default=None,
            metavar="FORMAT",
            dest="report_formats",
            help="Report type(s) to generate. If multiple, '--output' "
            "is treated as a prefix and the appropriate extensions are "
            "appended. If unspecified, the format is guessed from the "
            "file extension. Supported formats are: txt, json, yaml, "
            "pickle. See the documentation for a  full description of "
            "the structured output (json/yaml/pickle formats).",
        )

    @staticmethod
    def _validate_error_options(options: Namespace, parser: AtroposArgumentParser):
        options.report_file = options.output
