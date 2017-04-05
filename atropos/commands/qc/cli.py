"""Command line interface for the qc command.
"""
from atropos.commands.cli import (
    BaseCommandParser, parse_stat_args, configure_threads, positive, int_or_str,
    writeable_file)

class CommandParser(BaseCommandParser):
    name = 'qc'
    usage = """
atropos qc -se input.fastq
atropos qc -pe1 in1.fastq -pe2 in2.fastq
"""
    description = """
Compute read-level statistics. The output is identical to running the 'trim'
command with '--stats pre'. Use this command to get an idea of the quality of
your raw data.
"""

    def add_command_options(self):
        self.parser.set_defaults(
            action='qc',
            batch_size=None)
        
        group = self.add_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writeable_file, default="-", metavar="FILE",
            help="Write stats to file rather than stdout.")
        
        group = self.add_group(
            "Report", title="Report content and formatting options")
        group.add_argument(
            "--report-formats",
            nargs="*", choices=("txt", "json"), default=None, metavar="FORMAT",
            help="Report type(s) to generate. If multiple, '--output' "
                 "is treated as a prefix and the appropriate extensions are "
                 "appended. If unspecified, the format is guessed from the "
                 "output file extension.")
        group.add_argument(
            "--stats",
            type=parse_stat_args, default=None,
            help="Additional arguments for read statistic collection. E.g. "
                 "'pre:tiles' means to also collect tile-level statistics "
                 "(Illumina data only), and 'pre:tiles=<regexp>' means to use "
                 "the specified regular expression to extract key portions of "
                 "read names to collect the tile statistics.")
        
        group = self.add_group(
            "Parallel", title="Parallel (multi-core) options")
        group.add_argument(
            "-T",
            "--threads",
            type=positive(int, True), default=None, metavar="THREADS",
            help="Number of threads to use for read trimming. Set to 0 to use "
                 "max available threads. (Do not use multithreading)")
        group.add_argument(
            "--process-timeout",
            type=positive(int, True), default=60, metavar="SECONDS",
            help="Number of seconds process should wait before escalating "
                 "messages to ERROR level. (60)")
        group.add_argument(
            "--read-queue-size",
            type=int_or_str, default=None, metavar="SIZE",
            help="Size of queue for batches of reads to be processed. "
                 "(THREADS * 100)")
    
    def validate_command_options(self, options):
        options.report_file = options.output
        if options.threads is not None:
            threads = configure_threads(options, self.parser)
            if options.read_queue_size is None:
                options.read_queue_size = threads * 100
            elif (
                    options.read_queue_size > 0 and
                    options.read_queue_size < threads):
                self.parser.error("Read queue size must be >= than 'threads'")
        if options.batch_size is None:
            options.batch_size = 1000
