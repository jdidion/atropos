from abc import ABCMeta, abstractmethod
import copy
from multiprocessing import cpu_count
from pathlib import Path
import platform
import re
import textwrap
from typing import Dict, Sequence as SequenceType, Tuple, Type, Union, cast

from loguru import logger
from xphyle.paths import STDOUT, STDERR

from atropos import __version__
from atropos.commands import BaseCommand
from atropos.io import InputRead
from atropos.utils import (
    LOGGING_CONFIG, ReturnCode, classproperty
)
from atropos.utils.argparse import (
    ParagraphHelpFormatter
)
from atropos.utils.paths import splitext_compressed
from atropos.utils.argparse import (
    AtroposArgumentParser,
    Namespace,
    int_or_str,
    positive,
    probability,
    readable_file,
    writable_file,
)
from atropos.utils.ngs import ALPHABETS


class CommandConsole(metaclass=ABCMeta):
    @classproperty
    @abstractmethod
    def name(cls) -> str:
        """
        The command name.
        """

    @classmethod
    @abstractmethod
    def get_help(
        cls, fmt: str = "* {name}: {description}", wrap: int = 80, indent: int = 2
    ) -> str:
        """
        Returns a string to include in the command help.
        """

    @classmethod
    @abstractmethod
    def execute(cls, args: SequenceType[str] = ()) -> Tuple[ReturnCode, dict]:
        """
        Parse command line arguments, execute the command, and generate summary reports.

        Returns:
            A tuple (return_code, summary)
        """


class BaseCommandConsole(metaclass=ABCMeta):
    """
    Base class for Atropos commands.
    """

    @classproperty
    @abstractmethod
    def name(cls) -> str:
        """
        The command name.
        """

    @classproperty
    @abstractmethod
    def description(cls) -> str:
        """
        The command description.
        """

    @classproperty
    def usage(cls) -> str:
        """
        Command usage string.
        """
        return "atropos {command} [options]"

    @classproperty
    def preamble(cls) -> str:
        """
        Help message preamble.
        """
        return "Atropos version {version}"

    @classproperty
    def details(cls) -> str:
        """
        Help message details.
        """
        return ""

    @classmethod
    def get_description(cls, **kwargs) -> str:
        description = "{}\n\n{}\n\n{}".format(
            *(part.strip() for part in (cls.preamble, cls.description, cls.details))
        )
        return description.format(**kwargs)

    @classmethod
    def get_help(
        cls, fmt: str = "* {name}: {description}", wrap: int = 80, indent: int = 2
    ) -> str:
        """
        Returns a string to include in the command help.

        Args:
            fmt: The help string format.
            wrap: The number of characters at which to wrap lines.
            indent: The number of spaces to indent each line.
        """
        helpstr = fmt.format(name=cls.name, description=cls.description.strip())

        if wrap:
            return "\n".join(
                textwrap.wrap(
                    re.sub(r"\s+", " ", helpstr), wrap, subsequent_indent=" " * indent
                )
            )
        else:
            return helpstr

    @classmethod
    def execute(cls, args: SequenceType[str] = ()) -> Tuple[ReturnCode, dict]:
        """
        Parse command line arguments, execute the command, and generate summary reports.

        Returns:
            The tuple (return_code, summary)
        """
        options = cls._parse_args(args)
        command = cls._create_command(options)
        retcode, summary = command.run()
        if retcode is ReturnCode.SUCCESS and options.report_file:
            logger.debug(f"Writing report to {options.report_file}")
            cls.generate_reports(summary, options)
        else:
            logger.debug("Not generating report file")
        return retcode, summary

    @classmethod
    def _parse_args(cls, args: SequenceType[str]) -> Namespace:
        """
        Parses command line arguments.

        Args:
            args: command-line arguments

        Returns:
            A `Namespace`
        """
        parser = cls._create_argument_parser()
        cls._add_arguments(parser)

        options = parser.parse_args(args)
        options.orig_args = copy.copy(args)

        cls._setup_logging(options)
        cls._validate_options(options, parser)

        return options

    @classmethod
    def _create_command(cls, options: Namespace) -> BaseCommand:
        if issubclass(cls, BaseCommand):
            return cast(Type[BaseCommand], cls)(options)
        else:
            raise NotImplementedError()

    @classmethod
    def _create_argument_parser(cls) -> AtroposArgumentParser:
        """
        Creates an `ArgumentParser` to parse command line arguments to this command.

        Returns:
            An `ArgumentParser`
        """
        format_args = dict(name=cls.name, version=__version__)
        return AtroposArgumentParser(
            prog=f"atropos {cls.name}",
            usage=cls.usage.format(**format_args),
            description=cls.get_description(**format_args),
            formatter_class=ParagraphHelpFormatter,
        )

    @classmethod
    @abstractmethod
    def _add_arguments(cls, parser: AtroposArgumentParser) -> None:
        """
        Adds arguments to the ArgumentParser.

        Args:
            parser:
        """

    @classmethod
    def _validate_options(
        cls, options: Namespace, parser: AtroposArgumentParser
    ) -> None:
        """
        Validates command-line options.

        Args:
            options:
        """

    @classmethod
    def _setup_logging(cls, options: Namespace) -> None:
        """
        Setup logging and print an introductory message.

        Logging setup is only done if there are not already any handlers (can
        happen when this function is being called externally such as from unit
        tests).

        Args:
            options: Command-line options
        """
        if not LOGGING_CONFIG.was_setup:
            level = options.log_level or ("ERROR" if options.quiet else "INFO")

            if options.log_file is None:
                # Due to backwards compatibility, log messages are sent to
                # standard output instead of standard error if the -o option is
                # used.
                if options.output is STDERR:
                    log_path = STDOUT
                else:
                    log_path = STDERR
            else:
                log_path = options.log_file

            LOGGING_CONFIG.setup(log_path, level)

        logger.info(
            f"This is Atropos {__version__} with Python {platform.python_version()}"
        )

    @classmethod
    @abstractmethod
    def generate_reports(cls, summary: dict, options: Namespace) -> None:
        """
        Generates summary report(s) for this command.

        Args:
            summary:
            options:
        """


def add_common_options(parser: AtroposArgumentParser) -> None:
    """
    Adds common arguments to the parser.
    """
    parser.set_defaults(
        orig_args=None,
        paired=False,
        default_outfile=STDOUT,
        report_file=None,
        report_formats=None,
        batch_size=1000,
        counter_magnitude="M",
        ngstream_reader=None,
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Print debugging information. (no)",
    )
    parser.add_argument(
        "--progress",
        nargs="?",
        choices=("bar", "msg"),
        default=False,
        help="Show progress. An optional argument determines what type of progress "
        "meter to use: bar = progress bar; msg = status message. Ohterwise the default "
        "progress meter is shown. (now)"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Print only error messages. (no)",
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARN", "ERROR"),
        default=None,
        help="Logging level. (ERROR when --quiet else INFO)",
    )
    parser.add_argument(
        "--log-file",
        type=writable_file,
        default=None,
        metavar="FILE",
        help="File to write logging info. (stdout)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
        help="Show version information and exit.",
    )
    group = parser.add_group("Input")
    group.add_argument(
        "-pe1",
        "--input1",
        type=readable_file,
        default=None,
        metavar="FILE1",
        help="The first input file.",
    )
    group.add_argument(
        "-pe2",
        "--input2",
        type=readable_file,
        default=None,
        metavar="FILE2",
        help="The second input file.",
    )
    group.add_argument(
        "-l",
        "--interleaved-input",
        type=readable_file,
        default=None,
        metavar="FILE",
        help="Interleaved input file.",
    )
    group.add_argument(
        "-se",
        "--single-input",
        type=readable_file,
        default=None,
        metavar="FILE",
        help="A single-end read file.",
    )
    group.add_argument(
        "--single-input-read",
        type=int,
        dest="input_read",
        choices=(1, 2),
        default=None,
        help="When treating an interleaved FASTQ or paired-end SAM/BAM "
        "file as single-end, this option specifies which of the two "
        "reads to process. (both reads used)",
    )
    group.add_argument(
        "-sq",
        "--single-quals",
        type=readable_file,
        default=None,
        metavar="FILE",
        help="A single-end qual file.",
    )
    group.add_argument(
        "--accession",
        default=None,
        metavar="ACCN",
        help="Accession to stream using a supported protocol. Should be of the form "
        "<protocol>:<accession>, e.g. 'sra:SRR000066'. If no protocol is specified, "
        "it assumed to be an SRA accession."
    )
    group.add_argument(
        "--query",
        default=None,
        metavar="URL",
        help="Query URL for a supported protocol. Should be of the form "
        "'<protocol>+http://...'. If no protocol is specified, it is assumed to be "
        "htsget."
    )
    group.add_argument(
        "-f",
        "--input-format",
        choices=("fasta", "fastq", "sra-fastq", "sam", "bam"),
        default=None,
        help="Input file format. Ignored when reading csfasta/qual files. "
        "(auto-detect from file name extension)",
    )
    group.add_argument(
        "-Q",
        "--quality-base",
        type=positive(),
        default=33,
        help="Assume that quality values in FASTQ are encoded as "
        "ascii(quality + QUALITY_BASE). This needs to be set to 64 "
        "for some old Illumina FASTQ files. (33)",
    )
    group.add_argument(
        "-c",
        "--colorspace",
        action="store_true",
        default=False,
        help="Enable colorspace mode: Also trim the color that is adjacent "
        "to the found adapter. (no)",
    )
    group.add_argument(
        "--max-reads",
        type=int_or_str,
        default=None,
        metavar="N",
        help="Maximum number of reads/pairs to process (no max)",
    )
    group.add_argument(
        "--subsample",
        type=probability,
        default=None,
        metavar="PROB",
        help="Subsample a fraction of reads. (no)",
    )
    group.add_argument(
        "--subsample-seed",
        type=int,
        default=None,
        metavar="SEED",
        help="The seed to use for the pseudorandom number generator. Using "
        "the same seed will result in the same subsampling of reads.",
    )
    group.add_argument(
        "--batch-size",
        type=int_or_str,
        metavar="SIZE",
        help="Number of records to process in each batch. (1000)",
    )
    group.add_argument(
        "-D",
        "--sample-id",
        default=None,
        metavar="ID",
        help="Optional sample ID. Added to the summary output.",
    )
    group.add_argument(
        "--alphabet",
        default=None,
        metavar="NAME",
        choices=tuple(ALPHABETS.keys()),
        help="Specify a sequence alphabet to use for validating inputs. "
        "Currently, only 'dna' is supported. (no validation)",
    )


def validate_common_options(options: Namespace, parser: AtroposArgumentParser) -> None:
    """
    Validates arguments to common options.
    """
    # Find out which 'mode' we need to use.
    # TODO: unit tests for SRA streaming
    if options.accession or options.query:
        if options.input_format is None:
            options.input_format = "fastq"
        elif options.input_format not in ("fastq", "sra-fastq", "sam", "bam"):
            raise ValueError(
                f"Invalid file format for accession/query: {options.input_format}"
            )

        if options.accession:
            accession = options.accession

            if ":" in accession:
                protocol, accession = accession.split(":", 1)
            else:
                protocol = "sra"

            logger.debug(
                f"Opening reader for {protocol} Accession {options.accession}"
            )
        else:
            accession = options.query

            if "+" in accession:
                protocol, accession = accession.split("+", 1)
            else:
                protocol = "htsget"

            logger.debug(
                f"Opening reader for {protocol} query {options.accession}"
            )

        try:
            import ngstream

            options.ngstream_reader = ngstream.open(
                accession, protocol, batch_size=options.batch_size or 1000
            )
            options.ngstream_reader.start()
            options.paired = options.ngstream_reader.paired
        except Exception:
            logger.exception(
                f"Error while fetching accession {options.sra_accession} from SRA"
            )
            parser.error(f"Unable to read from accession {options.sra_accession}")
    elif options.single_input:
        if options.input1 or options.input2 or options.interleaved_input:
            parser.error("Cannot use -se together with -pe1, -pe2, or -l")

        options.paired = False
        options.input1 = options.single_input
        options.input2 = options.single_quals
    elif options.interleaved_input and options.input_read:
        options.input1 = options.interleaved_input
        options.paired = False
    else:
        if not options.interleaved_input and not (options.input1 and options.input2):
            parser.error(
                "Both '-pe1' and '-pe2' are required for paired-end "
                "trimming. If this is an interleaved file, use '-l' "
                "instead."
            )

        options.paired = True

    if options.input_read is None:
        options.input_read = InputRead.PAIRED if options.paired else InputRead.SINGLE

    # Set sample ID from the input file name(s)
    if options.sample_id is None:
        if options.ngstream_reader:
            options.sample_id = options.ngstream_reader.name
        else:
            fname = Path(options.input1 or options.interleaved_input).name
            name = splitext_compressed(fname)[0]

            if options.input2:
                name2 = splitext_compressed(Path(options.input2).name)[0]
                if name != name2:
                    for i in range(min(len(name), len(name2))):
                        if name[i] != name2[i]:
                            name = name[:i]
                            break

            if name.endswith("."):
                name = name[:-1]

            options.sample_id = name

    if options.quiet:
        options.progress = False
    else:
        if options.progress is None:
            options.progress = True

        if options.progress and options.output == STDERR:
            logger.warning("Progress bar may corrupt output written to STDERR")

    if options.report_file in (STDOUT, STDERR) and options.quiet:
        logger.warning("Quiet mode - report will not be written to stdout")
        options.report_file = None


def configure_threads(options, parser) -> int:
    """
    Determines the number of threads to use from the command-line options. Updates
    the value in options, and returns the number of threads.
    """
    if options.debug:
        parser.error("Cannot use debug mode with multiple threads")

    threads = options.threads

    if threads <= 0:
        threads = cpu_count()
    elif threads == 1:
        parser.error("--threads must be >= 2")

    options.threads = threads

    return threads


def parse_metrics_args(args_str) -> Dict[str, Union[bool, str]]:
    """
    Parses the optional value to the '--stat' option.
    """
    args = {}

    for arg in args_str.split(";"):
        arg_parts = arg.split("=")
        if len(arg_parts) == 1:
            args[arg_parts[0]] = True
        else:
            args[arg_parts[0]] = arg_parts[1]

    return args
