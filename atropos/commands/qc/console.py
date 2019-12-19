from atropos.commands.console import (
    BaseCommandConsole,
    add_common_options,
    validate_common_options,
    parse_metrics_args,
    configure_threads,
)
from atropos.commands.qc import QcCommand
from atropos.commands.legacy_reports import LegacyReportGenerator
from atropos.utils import classproperty
from atropos.utils.argparse import (
    AtroposArgumentParser,
    Namespace,
    positive,
    int_or_str,
    writable_file,
)


class QcCommandConsole(QcCommand, LegacyReportGenerator, BaseCommandConsole):
    @classproperty
    def description(cls) -> str:
        return """
        atropos qc -se input.fastq
        atropos qc -pe1 in1.fastq -pe2 in2.fastq
        """

    @classproperty
    def usage(cls) -> str:
        return """
        Compute read-level statistics. The output is identical to running the 'trim'
        command with '--metrics pre'. Use this command to get an idea of the quality of
        your raw data.
        """

    @classmethod
    def _add_arguments(cls, parser: AtroposArgumentParser) -> None:
        add_common_options(parser)
        cls._add_qc_options(parser)

    @classmethod
    def _validate_options(
        cls, options: Namespace, parser: AtroposArgumentParser
    ) -> None:
        validate_common_options(options, parser)
        cls._validate_qc_options(options, parser)

    @staticmethod
    def _add_qc_options(parser: AtroposArgumentParser):
        parser.set_defaults(action="qc", batch_size=None)
        group = parser.add_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writable_file,
            default="-",
            metavar="FILE",
            help="Write metrics to file rather than stdout.",
        )
        group = parser.add_group(
            "Report", title="Report content and formatting options"
        )
        group.add_argument(
            "--report-formats",
            nargs="*",
            choices=("txt", "json"),
            default=None,
            metavar="FORMAT",
            help="Report type(s) to generate. If multiple, '--output' "
            "is treated as a prefix and the appropriate extensions are "
            "appended. If unspecified, the format is guessed from the "
            "output file extension.",
        )
        group.add_argument(
            "--metrics",
            type=parse_metrics_args,
            default=None,
            help="Additional arguments for read statistic collection. E.g. "
            "'pre:tiles' means to also collect tile-level statistics "
            "(Illumina data only), and 'pre:tiles=<regexp>' means to use "
            "the specified regular expression to extract key portions of "
            "read names to collect the tile statistics.",
        )
        group = parser.add_group("Parallel", title="Parallel (multi-core) options")
        group.add_argument(
            "-T",
            "--threads",
            type=positive(int, True),
            default=None,
            metavar="THREADS",
            help="Number of threads to use for read trimming. Set to 0 to use "
            "max available threads. (Do not use multithreading)",
        )
        group.add_argument(
            "--process-timeout",
            type=positive(int, True),
            default=60,
            metavar="SECONDS",
            help="Number of seconds process should wait before escalating "
            "messages to ERROR level. (60)",
        )
        group.add_argument(
            "--read-queue-size",
            type=int_or_str,
            default=None,
            metavar="SIZE",
            help="Size of queue for batches of reads to be processed. "
            "(THREADS * 100)",
        )

    @staticmethod
    def _validate_qc_options(options: Namespace, parser: AtroposArgumentParser):
        options.report_file = options.output

        if options.threads is not None:
            threads = configure_threads(options, parser)
            if options.read_queue_size is None:
                options.read_queue_size = threads * 100
            elif 0 < options.read_queue_size < threads:
                parser.error("Read queue size must be >= than 'threads'")

        if options.batch_size is None:
            options.batch_size = 1000
