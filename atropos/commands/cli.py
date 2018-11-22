"""Base class and functions for implementing command-line interfaces for
Atropos subcommands.
"""
from argparse import ArgumentParser, ArgumentError, HelpFormatter
import copy
import logging
from multiprocessing import cpu_count
import operator
import os
import platform
import re
import sys
import textwrap
import urllib
from atropos import __version__
from atropos.io import STDOUT, STDERR, resolve_path, check_path, check_writeable
from atropos.io.compression import splitext_compressed
from atropos.io.seqio import SINGLE, PAIRED
from atropos.util import MAGNITUDE, ALPHABETS

class BaseCommandParser(object):
    """Base class for Atropos sub-commands.

    Subclasses must define name, description, and usage members.
    """
    preamble = "Atropos version {version}"
    usage = "atropos {command} [options]"
    description = ''
    details = ''

    def __init__(self):
        self.groups = {}
        self.create_parser()
        self.add_common_options()
        self.add_command_options()

    def parse(self, args):
        """Parse args using the conifgured ArgumentParser.

        Args:
            args: Command line arguments.

        Returns:
            A Namespace-like object.
        """
        options = self.parser.parse_args(args)
        options.orig_args = copy.copy(args)
        self.setup_logging(options)
        self.validate_common_options(options)
        self.validate_command_options(options)
        return options

    def create_parser(self):
        """Create the ArgumentParser.
        """
        format_args = dict(
            name=self.name,
            version=__version__)
        self.parser = ArgumentParser(
            prog="atropos {}".format(self.name),
            usage=self.usage.format(**format_args),
            description=self.get_description(**format_args),
            formatter_class=ParagraphHelpFormatter)

    def get_description(self, **kwargs):
        description = "{}\n\n{}\n\n{}".format(*(
            part.strip()
            for part in (self.preamble, self.description, self.details)))
        return description.format(**kwargs)

    def add_group(
            self, name, title=None, description=None, mutex=False,
            required=False):
        """Add a group to the parser. The group will be stored under `name`
        and can later be retrieved via `get_group`.
        """
        if name in self.groups:
            raise ValueError("Group already exists: {}".format(name))
        if mutex:
            group = self.parser.add_mutually_exclusive_group(required)
        else:
            group = self.parser.add_argument_group(title or name, description)
        self.groups[name] = group
        return group

    def get_group(self, name):
        """If a group has already been created with `name`, return the group,
        otherwise create a new group with that name.
        """
        if name in self.groups:
            return self.groups[name]
        else:
            return self.add_group(name)

    def add_common_options(self):
        """Add common arguments to the parser.
        """
        self.parser.set_defaults(
            orig_args=None,
            paired=False,
            default_outfile=STDOUT,
            report_file=None,
            report_formats=None,
            batch_size=1000,
            counter_magnitude="M",
            sra_reader=None)
        self.parser.add_argument(
            "--debug",
            action='store_true', default=False,
            help="Print debugging information. (no)")
        self.parser.add_argument(
            "--progress",
            choices=('bar', 'msg'), default=None,
            help="Show progress. bar = show progress bar; msg = show a status "
                 "message. (no)")
        self.parser.add_argument(
            "--quiet",
            action='store_true', default=False,
            help="Print only error messages. (no)")
        self.parser.add_argument(
            "--log-level",
            choices=('DEBUG', 'INFO', 'WARN', 'ERROR'), default=None,
            help="Logging level. (ERROR when --quiet else INFO)")
        self.parser.add_argument(
            "--log-file",
            type=writeable_file, default=None, metavar="FILE",
            help="File to write logging info. (stdout)")
        self.parser.add_argument(
            '--version', action='version', version=__version__,
            help="Show version information and exit.")

        group = self.add_group("Input")
        group.add_argument(
            "-pe1",
            "--input1",
            type=readable_file, default=None, metavar="FILE1",
            help="The first input file.")
        group.add_argument(
            "-pe2",
            "--input2",
            type=readable_file, default=None, metavar="FILE2",
            help="The second input file.")
        group.add_argument(
            "-l",
            "--interleaved-input",
            type=readable_file, default=None, metavar="FILE",
            help="Interleaved input file.")
        group.add_argument(
            "-se",
            "--single-input",
            type=readable_file, default=None, metavar="FILE",
            help="A single-end read file.")
        group.add_argument(
            "--single-input-read",
            type=int, dest='input_read', choices=(1, 2), default=None,
            help="When treating an interleaved FASTQ or paired-end SAM/BAM "
                 "file as single-end, this option specifies which of the two "
                 "reads to process. (both reads used)")
        group.add_argument(
            "-sq",
            "--single-quals",
            type=readable_file, default=None, metavar="FILE",
            help="A single-end qual file.")
        group.add_argument(
            "-sra",
            "--sra-accession",
            default=None, metavar="ACCN",
            help="Accession to stream from SRA (requires optional NGS "
                 "dependency to be installed).")
        group.add_argument(
            "-f",
            "--format",
            choices=('fasta','fastq','sra-fastq','sam','bam'), default=None,
            help="Input file format. Ignored when reading csfasta/qual files. "
                 "(auto-detect from file name extension)")
        group.add_argument(
            "-Q",
            "--quality-base",
            type=positive(), default=33,
            help="Assume that quality values in FASTQ are encoded as "
                 "ascii(quality + QUALITY_BASE). This needs to be set to 64 "
                 "for some old Illumina FASTQ files. (33)")
        group.add_argument(
            "-c",
            "--colorspace",
            action='store_true', default=False,
            help="Enable colorspace mode: Also trim the color that is adjacent "
                 "to the found adapter. (no)")
        group.add_argument(
            "--max-reads",
            type=int_or_str, default=None, metavar="N",
            help="Maximum number of reads/pairs to process (no max)")
        group.add_argument(
            "--subsample",
            type=probability, default=None, metavar="PROB",
            help="Subsample a fraction of reads. (no)")
        group.add_argument(
            "--subsample-seed",
            type=int, default=None, metavar="SEED",
            help="The seed to use for the pseudorandom number generator. Using "
                 "the same seed will result in the same subsampling of reads.")
        group.add_argument(
            "--batch-size",
            type=int_or_str, metavar="SIZE",
            help="Number of records to process in each batch. (1000)")
        group.add_argument(
            "-D",
            "--sample-id",
            default=None, metavar="ID",
            help="Optional sample ID. Added to the summary output.")
        group.add_argument(
            "--alphabet",
            default=None, metavar="NAME", choices=tuple(ALPHABETS.keys()),
            help="Specify a sequence alphabet to use for validating inputs. "
                 "Currently, only 'dna' is supported. (no validation)")

    def add_command_options(self):
        """Add command-specific options. At the very least,
        "-o, --output" is required.
        """
        raise NotImplementedError()

    def setup_logging(self, options):
        """Setup logging and print an introductory message.

        Logging setup is only done if there are not already any handlers (can
        happen when this function is being called externally such as from unit
        tests).
        """
        if not logging.root.handlers:
            level = options.log_level or ("ERROR" if options.quiet else "INFO")
            level = getattr(logging, level)
            if options.log_file is None:
                # Due to backwards compatibility, log messages are sent to
                # standard output instead of standard error if the -o option is
                # used.
                stream = sys.stdout
                if options.output in (None, STDOUT, STDERR):
                    stream = sys.stderr
                handler = logging.StreamHandler(stream)
            else:
                handler = logging.FileHandler(options.log_file)
            handler.setFormatter(logging.Formatter(
                '%(asctime)s %(levelname)s: %(message)s'))
            handler.setLevel(level)
            logging.getLogger().setLevel(level)
            logging.getLogger().addHandler(handler)

        logger = logging.getLogger()
        logger.info(
            "This is Atropos %s with Python %s", __version__,
            platform.python_version())

    def validate_common_options(self, options):
        """Validate arguments to common options.
        """
        parser = self.parser

        # Find out which 'mode' we need to use.
        # TODO: unit tests for SRA streaming
        # TODO: add srastream to pypi
        if options.sra_accession:
            if options.format not in ('fastq', 'sam', 'bam', None):
                raise ValueError(
                    "Invalid file format for SRA accession: {}".format(
                        options.format))
            options.format = 'fastq'
            logging.getLogger().debug(
                "Opening reader for SRA Accession {}".format(
                    options.sra_accession))
            try:
                from srastream import SraReader
                options.sra_reader = SraReader(
                    options.sra_accession,
                    batch_size=options.batch_size or 1000)
                options.sra_reader.start()
                options.paired = options.sra_reader.paired
            except Exception:
                logging.getLogger().exception(
                    "Error while fetching accession {} from SRA".format(
                        options.sra_accession))
                parser.error("Unable to read from accession {}".format(
                    options.sra_accession))
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
            if (not options.interleaved_input and (
                    not options.input1 or not options.input2)):
                parser.error(
                    "Both '-pe1' and '-pe2' are required for paired-end "
                    "trimming. If this is an interleaved file, use '-l' "
                    "instead.")
            options.paired = True

        if options.input_read is None:
            options.input_read = PAIRED if options.paired else SINGLE

        # Set sample ID from the input file name(s)
        if options.sample_id is None:
            if options.sra_reader:
                options.sample_id = options.sra_reader.name
            else:
                fname = os.path.basename(
                    options.input1 or options.interleaved_input)
                name = splitext_compressed(fname)[0]
                if options.input2:
                    name2 = splitext_compressed(os.path.basename(options.input2))[0]
                    if name != name2:
                        for i in range(min(len(name), len(name2))):
                            if name[i] != name2[i]:
                                name = name[:i]
                                break
                if name.endswith('.'):
                    name = name[:-1]
                options.sample_id = name

        if options.quiet:
            options.progress = None
        elif options.progress and options.output == STDERR:
            logging.getLogger().warning(
                "Progress bar may corrupt output written to STDERR")

        if options.report_file in (STDOUT, STDERR) and options.quiet:
            logging.getLogger().warning(
                "Quiet mode - report will not be written to stdout")
            options.report_file = None

    def validate_command_options(self, options):
        """Validate command-specific options.
        """
        pass

# Extensions to argparse

class ParagraphHelpFormatter(HelpFormatter):
    """HelpFormatter that wraps text in paragraphs.
    """
    def _fill_text(self, text, width, indent):
        text = re.sub('[ \t]{2,}', ' ', text)
        paragraphs = [
            textwrap.fill(
                p, width, initial_indent=indent, subsequent_indent=indent)
            for p in re.split("\n\n", text)]
        return "\n\n".join(paragraphs)

class TypeWithArgs(object): # pylint: disable=no-member
    """A data type that requires additional static arguments.
    """
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, string):
        return self._do_call(string, *self.args, **self.kwargs) or string

    def _do_call(self, string, *args, **kwargs):
        """Convert and/or validate `string`.
        """
        raise NotImplementedError()

class CompositeType(object):
    """A composite of multiple data types.
    """
    def __init__(self, *types):
        self.types = types

    def __call__(self, string):
        result = string
        for datatype in self.types:
            result = datatype(result)
        return result

class ComparisonValidator(TypeWithArgs):
    """Validator that compares an argument against an expected value.
    """
    def _do_call(self, lhs, rhs, oper, expected=True):
        assert oper(lhs, rhs) == expected, "{}({}, {}) != {}".format(
            oper, lhs, rhs, expected)

class CharList(object):
    """Parses a string into a list of characters and ensures they are all in
    the choices tuple.
    """
    def __init__(self, choices):
        self.choices = set(choices)

    def __call__(self, string):
        chars = list(string)
        assert all(char in self.choices for char in chars)
        return chars

class Delimited(TypeWithArgs):
    """Splits a string argument using a delimiter.
    """
    def _do_call(
            self, string, delim=",", data_type=None, choices=None,
            min_len=None, max_len=None):
        if isinstance(string, str):
            vals = string.split(delim) if delim else (string,)
        else:
            vals = string

        if vals[0] == "*" and choices is not None:
            vals = choices

        if data_type:
            vals = [data_type(v) for v in vals]

        if min_len and len(vals) < min_len:
            raise ArgumentError(
                self, "there must be at least {} values".format(min_len))

        if max_len and len(vals) > max_len:
            raise ArgumentError(
                self, "there can be at most {} values".format(max_len))

        return vals

ACCESS = dict(
    r=os.R_OK,
    rU=os.R_OK,
    rb=os.R_OK,
    w=os.W_OK,
    wb=os.W_OK,
    x=os.X_OK
)

class AccessiblePath(TypeWithArgs):
    """Test that a path is accessible.
    """
    def _do_call(self, path, type_, mode):
        if type_ == 'f' and path in (STDOUT, STDERR):
            return path
        if 'w' in mode:
            return check_writeable(path, type_)
        else:
            return check_path(path, type_, ACCESS[mode])

class ReadwriteableFile(object):
    """Validator for a file argument that must be both readable and writeable.
    """
    def __init__(self):
        self.read_type = AccessiblePath('f', 'r')
        self.write_type = AccessiblePath('f', 'w')

    def __call__(self, string):
        path = string
        if os.path.exists(path):
            path = self.read_type(path)
        path = self.write_type(path)
        return path

def existing_path(path):
    """Test that a path exists."""
    if path == STDOUT:
        return path
    return resolve_path(path)

readable_file = CompositeType(existing_path, AccessiblePath('f', 'r'))
"""Test that a file exists and is readable."""

writeable_file = AccessiblePath('f', 'w')
"""Test that a file 1) exists and is writelable, or 2) does not exist but
is in a writeable directory.
"""

readwriteable_file = ReadwriteableFile()
"""Test that a file is both readable and writeable."""

def readable_url(url):
    """Validator for a URL that must be readable.

    Args:
        url: The URL to validate.
    """
    parsed = urllib.parse.urlparse(url)
    scheme = parsed.scheme or 'file'
    if scheme == 'file':
        filename = readable_file(parsed.path)
        return 'file:' + filename
    else:
        return url

str_list = Delimited(data_type=str)
"""Comma-delimited list of strings."""

INT_OR_STR_RE = re.compile(r"([\d\.]+)([KkMmGg]?)")

def int_or_str(arg):
    """Similar to int(), but accepts K, M, and G abbreviations.
    """
    if arg is None or isinstance(arg, int):
        return arg
    elif isinstance(arg, str):
        match = INT_OR_STR_RE.match(arg.upper())
        num, mult = match.groups()
        if mult:
            return int(float(num) * MAGNITUDE[mult])
        else:
            return int(num)
    else:
        raise ValueError("Unsupported type {}".format(arg))

def positive(type_=int, inclusive=False):
    """Test that a number is greater than (or equal to, if ``inclusive=True``)
    zero.
    """
    oper = operator.ge if inclusive else operator.gt
    return CompositeType(type_, ComparisonValidator(0, oper))

def between(min_val=None, max_val=None, type_=int):
    """Returns a CompositeType that validates `min_val <= x <= max_val`.
    """
    return CompositeType(
        type_,
        ComparisonValidator(min_val, operator.ge),
        ComparisonValidator(max_val, operator.le))

probability = between(0, 1, float)
"""A float between 0-1 (inclusive)."""

def configure_threads(options, parser):
    """Determine the number of threads to use from the command-line options.
    Updates the value in options, and returns the number of threads.
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

def parse_stat_args(args_str):
    """Parse the optional value to the '--stat' option.
    """
    args = {}
    for arg in args_str.split(';'):
        arg_parts = arg.split('=')
        if len(arg_parts) == 1:
            args[arg_parts[0]] = True
        else:
            args[arg_parts[0]] = arg_parts[1]
    return args
