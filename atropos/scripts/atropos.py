#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-spaces all;

from argparse import ArgumentParser, HelpFormatter
from collections import OrderedDict
import copy
import logging
import operator
import os
import platform
import re
import sys
import time
import textwrap

# Print a helpful error message if the extension modules cannot be imported.
from atropos import *
check_importability()

from atropos import __version__
import atropos.commands
from atropos.util import MAGNITUDE
from atropos.xopen import STDOUT, STDERR, resolve_path, check_path, check_writeable

# Extensions to argparse

class ParagraphHelpFormatter(HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = re.sub('[ \t]{2,}', ' ', text)
        paragraphs = [
            textwrap.fill(p, width, initial_indent=indent, subsequent_indent=indent)
            for p in re.split("\n\n", text)
        ]
        return "\n\n".join(paragraphs)

class TypeWithArgs(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, string):
        return self._do_call(string, *self.args, **self.kwargs) or string

class CompositeType(object):
    def __init__(self, *types):
        self.types = types

    def __call__(self, string):
        result = string
        for t in self.types:
            result = t(result)
        return result

class ComparisonValidator(TypeWithArgs):
    def _do_call(self, lhs, rhs, oper, expected=True):
        assert oper(lhs, rhs) == expected, "{}({}, {}) != {}".format(oper, lhs, rhs, expected)

class char_list(object):
    def __init__(self, choices):
        self.choices = set(choices)
    
    def __call__(self, string):
        l = list(string)
        assert all(c in self.choices for c in l)
        return l

class delimited(TypeWithArgs):
    """Splits a string argument using a delimiter."""
    def _do_call(self, string, delim=",", data_type=None, choices=None, min_len=None, max_len=None):
        if isinstance(string, str):
            vals = string.split(delim) if delim else (string,)
        else:
            vals = string
        
        if vals[0] == "*" and choices is not None:
            vals = choices
        
        if data_type:
            vals = [data_type(v) for v in vals]
        
        if min_len and len(vals) < min_len:
            raise ArgumentError(self, "there must be at least {} values".format(min_len))
        
        if max_len and len(vals) > max_len:
            raise ArgumentError(self, "there can be at most {} values".format(max_len))
        
        return vals

ACCESS = dict(r=os.R_OK, rU=os.R_OK, rb=os.R_OK, w=os.W_OK, wb=os.W_OK, x=os.X_OK)

class accessible_path(TypeWithArgs):
    """Test that a path is accessible"""
    def _do_call(self, path, type_, mode):
        if type_ == 'f' and path in (STDOUT, STDERR):
            return path
        if 'w' in mode:
            return check_writeable(path, type_)
        else:
            return check_path(path, type_, ACCESS[mode])

def existing_path(path):
    """Test that a path exists."""
    if path == STDOUT:
        return path
    return resolve_path(path)

readable_file = CompositeType(existing_path, accessible_path('f', 'r'))
"""Test that a file exists and is readable."""

writeable_file = accessible_path('f', 'w')
"""
Test that a file 1) exists and is writelable, or 2) does not exist but
is in a writeable directory.
"""

class _readwriteable_file(object):
    def __init__(self):
        self.r = accessible_path('f', 'r')
        self.w = accessible_path('f', 'w')
    
    def __call__(self, string):
        path = string
        if os.path.exists(path):
            path = self.r(path)
        path = self.w(path)
        return path
readwriteable_file = _readwriteable_file()
"""Test that a file is both readable and writeable."""

import urllib
def readable_url(url):
    parsed = urllib.parse.urlparse(url)
    scheme = parsed.scheme or 'file'
    if scheme == 'file':
        filename = readable_file(parsed.path)
        return 'file:' + filename
    else:
        return url

str_list = delimited(data_type=str)
"""Comma-delimited list of strings."""

int_or_str_re = re.compile("([\d\.]+)([KkMmGg]?)")
def int_or_str(x):
    """Similar to int(), but accepts K, M, and G abbreviations"""
    if x is None or isinstance(x, int):
        return x
    elif isinstance(x, str):
        match = int_or_str_re.match(x.upper())
        num, mult = match.groups()
        if mult:
            return int(float(num) * MAGNITUDE[mult])
        else:
            return int(num)
    else:
        raise Exception("Unsupported type {}".format(x))

def positive(type_=int, inclusive=False):
    """Test that a number is greater than (or equal to, if ``inclusive=True``) zero."""
    return ge(0, type_) if inclusive else gt(0, type_)

def gt(x, type_=int):
    """Test that a number is greater than another number."""
    return CompositeType(type_, ComparisonValidator(x, operator.gt))

def ge(x, type_=int):
    """Test that a number is greater than another number."""
    return CompositeType(type_, ComparisonValidator(x, operator.ge))

def between(min_val=None, max_val=None, type_=int):
    return CompositeType(
        type_,
        ComparisonValidator(min_val, operator.ge),
        ComparisonValidator(max_val, operator.le))

probability = between(0, 1, float)
"""A float between 0-1 (inclusive)."""

# Commands

COMMANDS = OrderedDict()

class Command(object):
    def __init__(self, args):
        self.orig_args = copy.copy(args)
        self.create_parser()
        self.add_common_options()
        self.add_command_options()
        self.options = self.parser.parse_args(args)
        self.setup_logging()
        self.validate_common_options()
        self.validate_command_options()
    
    def create_parser(self):
        self.parser = ArgumentParser(
            usage=self.usage,
            description=self.description.format(version=__version__),
            formatter_class=ParagraphHelpFormatter)
    
    def add_common_options(self):
        # Add common options
        self.parser.set_defaults(
            paired=False,
            default_outfile=STDOUT,
            batch_size=1000)
        self.parser.add_argument(
            "--debug",
            action='store_true', default=False,
            help="Print debugging information. (no)")
        self.parser.add_argument(
            "--progress",
            choices=('bar', 'msg'), default=None,
            help="Show progress. bar = show progress bar; msg = show a status message. (no)")
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
        
        group = self.parser.add_argument_group("Input")
        group.add_argument(
            "-pe1",
            "--input1",
            type=readable_file, default=None, metavar="FILE1",
            help="The first (and possibly only) input file.")
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
            type=int, choices=(1, 2), default=None,
            help="When treating an interleaved FASTQ or paired-end SAM/BAM "
                "file as single-end, this option specifies which of the two "
                "reads to process. (both reads used)")
        group.add_argument(
            "-sq",
            "--single-quals",
            type=readable_file, default=None, metavar="FILE",
            help="A single-end qual file.")
        group.add_argument(
            "-f",
            "--format",
            choices=('fasta','fastq','sra-fastq','sam','bam'), default=None,
            help="Input file format. Ignored when reading csfasta/qual files. (auto-detect from "
                 "file name extension)")
        group.add_argument(
            "-c",
            "--colorspace",
            action='store_true', default=False,
            help="Enable colorspace mode: Also trim the color that is adjacent to "
                 "the found adapter. (no)")
        group.add_argument(
            "--max-reads",
            type=int_or_str, default=None, metavar="N",
            help="Maximum number of reads/pairs to process (no max)")
        group.add_argument(
            "--subsample",
            type=probability, default=None, metavar="PROB",
            help="Subsample a fraction of reads. (no)")
        group.add_argument(
            "--batch-size",
            type=int_or_str, metavar="SIZE",
            help="Number of records to process in each batch. (5000)")
    
    def add_command_options(self):
        """
        Add command-specific options. At the very least,
        "-o, --output" is required.
        """
        raise NotImplemented()
    
    def setup_logging(self):
        # Setup logging only if there are not already any handlers (can happen when
        # this function is being called externally such as from unit tests)
        if not logging.root.handlers:
            level = self.options.log_level or ("ERROR" if self.options.quiet else "INFO")
            level = getattr(logging, level)
            if self.options.log_file is None:
                # Due to backwards compatibility, log messages are sent to standard output
                # instead of standard error if the -o option is used.
                handler = logging.StreamHandler(
                    sys.stdout if self.options.output not in (None, STDOUT, STDERR) else sys.stderr)
            else:
                handler = logging.FileHandler(options.log_file)
            handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
            handler.setLevel(level)
            logging.getLogger().setLevel(level)
            logging.getLogger().addHandler(handler)
        
        logger = logging.getLogger()
        logger.info("This is Atropos %s with Python %s", __version__, platform.python_version())
        logger.info("Command line parameters: %s %s", self.name, " ".join(self.orig_args))
    
    def validate_common_options(self):
        options = self.options
        parser = self.parser
        
        # Find out which 'mode' we need to use.
        if options.single_input:
            if options.input1 or options.input2 or options.interleaved_input:
                parser.error("Cannot use -se together with -pe1, -pe2, or -l")
            options.paired = False
            options.input1 = options.single_input
            options.input2 = options.single_quals
        elif options.interleaved_input and options.single_input_read:
            options.input1 = options.interleaved_input
            options.paired = False
        else:
            if not options.interleaved_input and (not options.input1 or not options.input2):
                parser.error("Both '-pe1' and '-pe2' are required for paired-end trimming. If this is an "
                    "interleaved file, use '-l' instead.")
            options.paired = True
        
        if options.quiet:
            options.progress = None
        elif options.progress and options.output == STDERR:
            parser.warn("Progress bar may corrupt output written to STDERR")
    
    def validate_command_options(self):
        pass
    
    def execute(self):
        cmd = getattr(atropos.commands, self.name)
        cmd(self.options, self.parser)

class TrimCommand(Command):
    """trim sequencing reads."""
    name = "trim"
    usage = """
atropos trim -a ADAPTER [options] [-o output.fastq] -se input.fastq
atropos trim -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq -pe1 in1.fastq -pe2 in2.fastq
"""
    description = """
Atropos version {version}

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. The reverse complement is *not* automatically
searched. All reads from input.fastq will be written to output.fastq with the
adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter
sequences can be given (use further -a options), but only the best-matching
adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.
"""
    
    def add_command_options(self):
        parser = self.parser
        parser.set_defaults(
            zero_cap=None,
            action='trim',
            batch_size=5000,
            known_adapter=None)
        
        group = parser.add_argument_group("Finding adapters",
            description="Parameters -a, -g, -b specify adapters to be removed from "
                "each read (or from the first read in a pair if data is paired). "
                "If specified multiple times, only the best matching adapter is "
                "trimmed (but see the --times option). When the special notation "
                "'file:FILE' is used, adapter sequences are read from the given "
                "FASTA file. When the --adapter-file option is used, adapters can "
                "be specified by name rather than sequence.")
        group.add_argument(
            "-a",
            "--adapter",
            action="append", default=[], metavar="ADAPTER", dest="adapters",
            help="Sequence of an adapter ligated to the 3' end (paired data: of the "
                "first read). The adapter and subsequent bases are trimmed. If a "
                "'$' character is appended ('anchoring'), the adapter is only "
                "found if it is a suffix of the read. (none)")
        group.add_argument(
            "-g",
            "--front",
            action="append", default=[], metavar="ADAPTER",
            help="Sequence of an adapter ligated to the 5' end (paired data: of the "
                "first read). The adapter and any preceding bases are trimmed. "
                "Partial matches at the 5' end are allowed. If a '^' character is "
                "prepended ('anchoring'), the adapter is only found if it is a "
                "prefix of the read. (none)")
        group.add_argument(
            "-b",
            "--anywhere",
            action="append", default=[], metavar="ADAPTER",
            help="Sequence of an adapter that may be ligated to the 5' or 3' end "
                "(paired data: of the first read). Both types of matches as "
                "described under -a und -g are allowed. If the first base of the "
                "read is part of the match, the behavior is as with -g, otherwise "
                "as with -a. This option is mostly for rescuing failed library "
                "preparations - do not use if you know which end your adapter was "
                "ligated to! (none)")
        group.add_argument(
            "-F",
            "--known-adapters-file",
            type=readable_file, action="append", default=None,
            help="Path or URL of a FASTA file containing adapter sequences.")
        group.add_argument(
            "--no-default-adapters",
            action="store_false", dest="default_adapters", default=True,
            help="Don't fetch the default adapter list (which is currently stored as a) "
                "GitHub gist).")
        group.add_argument(
            "--adapter-cache-file",
            type=readwriteable_file, default='.adapters',
            help="File where adapter sequences will be cached, unless "
                "--no-cache-adapters is set.")
        group.add_argument(
            "--no-cache-adapters",
            action="store_false", dest="cache_adapters", default=True,
            help="Don't cache adapters list as '.adapters' in the working directory.")
        group.add_argument(
            "--no-trim",
            action='store_const', dest='action', const=None,
            help="Match and redirect reads to output/untrimmed-output as usual, "
                "but do not remove adapters. (no)")
        group.add_argument(
            "--mask-adapter",
            action='store_const', dest='action', const='mask',
            help="Mask adapters with 'N' characters instead of trimming them. (no)")
        
        ## Arguments specific to the choice of aligner
        
        group.add_argument("--aligner", choices=('adapter', 'insert'), default='adapter',
            help="Which alignment algorithm to use for identifying adapters. Currently, "
                 "you can choose between the semi-global alignment strategy used in Cutdapt "
                 "('adapter') or the more accurate insert-based alignment algorithm ('insert'). "
                 "Note that insert-based alignment can only be used with paired-end reads "
                 "containing 3' adapters. New algorithms are being implemented and the default "
                 "is likely to change. (adapter)")
        
        # TODO: all the different matching options are pretty confusing. Either explain their
        # usage better in the docs or find a way to simplify the choices.
        
        # Arguments for adapter match
        group.add_argument(
            "-e",
            "--error-rate",
            type=probability, default=None,
            help="Maximum allowed error rate for adapter match (no. of errors divided by the length "
                "of the matching region). (0.1)")
        group.add_argument(
            "--indel-cost",
            type=positive(int, True), default=None, metavar="COST",
            help="Integer cost of insertions and deletions during adapter match. Substitutions always "
                 "have a cost of 1. (1)")
        group.add_argument(
            "--no-indels",
            action='store_false', dest='indels', default=True,
            help="Allow only mismatches in alignments. (allow both mismatches and indels)")
        group.add_argument(
            "-n",
            "--times",
            type=positive(int, False), default=1, metavar="COUNT",
            help="Remove up to COUNT adapters from each read. (1)")
        group.add_argument(
            "--match-read-wildcards",
            action="store_true", default=False,
            help="Interpret IUPAC wildcards in reads. (no)")
        group.add_argument(
            "-N",
            "--no-match-adapter-wildcards",
            action="store_false", dest='match_adapter_wildcards', default=True,
            help="Do not interpret IUPAC wildcards in adapters. (no)")
        group.add_argument(
            "-O",
            "--overlap",
            type=positive(int, False), default=None, metavar="MINLENGTH",
            help="If the overlap between the read and the adapter is shorter than "
                "MINLENGTH, the read is not modified. Reduces the no. of bases "
                "trimmed due to random adapter matches. (3)")
        group.add_argument(
            "--adapter-max-rmp",
            type=probability, default=None, metavar="PROB",
            help="If no minimum overlap (-O) is specified, then adapters are only matched "
                 "when the probabilty of observing k out of n matching bases is <= PROB. (1E-6)")
        
        # Arguments for insert match
        group.add_argument(
            "--insert-max-rmp",
            type=probability, default=1E-6, metavar="PROB",
            help="Overlapping inserts only match when the probablity of observing k of n "
                 "matching bases is <= PROB. (1E-6)")
        group.add_argument(
            "--insert-match-error-rate",
            type=probability, default=None,
            help="Maximum allowed error rate for insert match (no. of errors divided by the length "
                "of the matching region). (0.2)")
        group.add_argument(
            "--insert-match-adapter-error-rate",
            type=probability, default=None,
            help="Maximum allowed error rate for matching adapters after successful insert match "
                 "(no. of errors divided by the length of the matching region). (0.2)")
        
        # Arguments for merging and error correction
        # TODO: add RMP parameter for MergeOverlap
        group.add_argument(
            "-R",
            "--merge-overlapping",
            action="store_true", default=False,
            help="Merge read pairs that overlap into a single sequence. This is experimental. (no)")
        group.add_argument(
            "--merge-min-overlap",
            type=positive(float, True), default=0.9,
            help="The minimum overlap between reads required for merging. If this number is (0,1.0], "
                 "it specifies the minimum length as the fraction of the length of the *shorter* read "
                 "in the pair; otherwise it specifies the minimum number of overlapping base pairs ("
                 "with an absolute minimum of 2 bp). (0.9)")
        group.add_argument(
            "--merge-error-rate",
            type=probability, default=None,
            help="The maximum error rate for merging. (0.2)")
        group.add_argument(
            "--correct-mismatches",
            choices=("liberal", "conservative", "N"), default=None,
            help="How to handle mismatches while aligning/merging. 'Liberal' and 'conservative' error "
                 "correction both involve setting the base to the one with the best quality. They differ "
                 "only when the qualities are equal -- liberal means set it to the base from the read with "
                 "the overall best median base quality, while conservative means to leave it unchanged. "
                 "'N' means to set the base to N. If exactly one base is ambiguous, the non-ambiguous base "
                 "is always used. (no error correction)")
        
        group = parser.add_argument_group("Additional read modifications")
        group.add_argument(
            "--op-order",
            type=char_list(choices=('A','C','G','Q')), default="CGQA",
            help="The order in which trimming operations are be applied. This is a string of "
                 "1-4 of the following characters: A = adapter trimming; C = cutting "
                 "(unconditional); G = NextSeq trimming; Q = quality trimming. The default is "
                 "'CGQA' to maintain compatibility with Cutadapt; however, this is likely to "
                 "change to 'GACQ' in the near future.")
        group.add_argument(
            "-u",
            "--cut",
            type=int, action='append', default=[], metavar="LENGTH",
            help="Remove bases from each read (first read only if paired). "
                "If LENGTH is positive, remove bases from the beginning. "
                "If LENGTH is negative, remove bases from the end. "
                "Can be used twice if LENGTHs have different signs. (no)")
        group.add_argument(
            "-q",
            "--quality-cutoff",
            type=delimited(data_type=positive(int, True), min_len=1, max_len=2),
            default=None, metavar="[5'CUTOFF,]3'CUTOFF",
            help="Trim low-quality bases from 5' and/or 3' ends of each read before "
                "adapter removal. Applied to both reads if data is paired. If one "
                "value is given, only the 3' end is trimmed. If two "
                "comma-separated cutoffs are given, the 5' end is trimmed with "
                "the first cutoff, the 3' end with the second. (no)")
        group.add_argument(
            "-i",
            "--cut-min",
            type=int, action='append', default=[], metavar="LENGTH",
            help="Similar to -u, except that cutting is done AFTER adapter trimming, and only "
                "if a minimum of LENGTH bases was not already removed. (no)")
        group.add_argument(
            "--nextseq-trim",
            type=positive(), default=None, metavar="3'CUTOFF",
            help="NextSeq-specific quality trimming (each read). Trims also dark "
                "cycles appearing as high-quality G bases (EXPERIMENTAL). (no)")
        group.add_argument(
            "--quality-base",
            type=positive(), default=33,
            help="Assume that quality values in FASTQ are encoded as ascii(quality "
                "+ QUALITY_BASE). This needs to be set to 64 for some old Illumina "
                "FASTQ files. (33)")
        group.add_argument(
            "--trim-n",
             action='store_true', default=False,
            help="Trim N's on ends of reads. (no)")
        group.add_argument(
            "-x",
            "--prefix",
            default='',
            help="Add this prefix to read names. Use {name} to insert the name of "
                 "the matching adapter. (no)")
        group.add_argument(
            "-y",
            "--suffix",
            default='',
            help="Add this suffix to read names; can also include {name}. (no)")
        group.add_argument(
            "--strip-suffix",
            action='append', default=[],
            help="Remove this suffix from read names if present. Can be given "
                 "multiple times. (no)")
        group.add_argument(
            "--length-tag",
            metavar="TAG",
            help="Search for TAG followed by a decimal number in the description "
                "field of the read. Replace the decimal number with the correct "
                "length of the trimmed read. For example, use --length-tag 'length=' "
                "to correct fields like 'length=123'. (no)")
        
        group = parser.add_argument_group("Filtering of processed reads")
        group.add_argument(
            "--discard-trimmed", "--discard",
            action='store_true', default=False,
            help="Discard reads that contain an adapter. Also use -O to avoid "
                "discarding too many randomly matching reads! (no)")
        group.add_argument(
            "--discard-untrimmed", "--trimmed-only",
            action='store_true', default=False,
            help="Discard reads that do not contain the adapter. (no)")
        group.add_argument(
            "-m",
            "--minimum-length",
            type=positive(int, True), default=None, metavar="LENGTH",
            help="Discard trimmed reads that are shorter than LENGTH. Reads that "
                "are too short even before adapter removal are also discarded. In "
                "colorspace, an initial primer is not counted. (0)")
        group.add_argument(
            "-M",
            "--maximum-length",
            type=positive(int, True), default=sys.maxsize, metavar="LENGTH",
            help="Discard trimmed reads that are longer than LENGTH. "
                "Reads that are too long even before adapter removal "
                "are also discarded. In colorspace, an initial primer "
                "is not counted. (no limit)")
        group.add_argument(
            "--max-n",
            type=positive(float, True), default=None, metavar="COUNT",
            help="Discard reads with too many N bases. If COUNT is an integer, it "
                "is treated as the absolute number of N bases. If it is between 0 "
                "and 1, it is treated as the proportion of N's allowed in a read. "
                "(no)")
        
        group = parser.add_argument_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writeable_file, metavar="FILE",
            help="Write trimmed reads to FILE. FASTQ or FASTA format is chosen "
                "depending on input. The summary report is sent to standard output. "
                "Use '{name}' in FILE to demultiplex reads into multiple "
                "files. (write to standard output)")
        group.add_argument(
            "--info-file",
            type=writeable_file, metavar="FILE",
            help="Write information about each read and its adapter matches into FILE. "
                 "See the documentation for the file format. (no)")
        group.add_argument(
            "-r",
            "--rest-file",
            type=writeable_file, metavar="FILE",
            help="When the adapter matches in the middle of a read, write the "
                "rest (after the adapter) into FILE. (no)")
        group.add_argument(
            "--wildcard-file",
            type=writeable_file, metavar="FILE",
            help="When the adapter has N bases (wildcards), write adapter bases "
                "matching wildcard positions to FILE. When there are indels in the "
                "alignment, this will often not be accurate. (no)")
        group.add_argument(
            "--too-short-output",
            type=writeable_file, metavar="FILE",
            help="Write reads that are too short (according to length specified by "
                "-m) to FILE. (no - too short reads are discarded)")
        group.add_argument(
            "--too-long-output",
            type=writeable_file, metavar="FILE",
            help="Write reads that are too long (according to length specified by "
                "-M) to FILE. (no - too long reads are discarded)")
        group.add_argument(
            "--untrimmed-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write reads that do not contain the adapter to FILE. (no - "
                "untrimmed reads are written to default output)")
        group.add_argument(
            "--merged-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write reads that have been merged to this file. (merged reads "
                "are discarded)")
        group.add_argument(
            "--report-file",
            type=writeable_file, default=None, metavar="FILE",
            help="Write report to file rather than stdout/stderr. (no)")
        
        group = parser.add_argument_group("Colorspace options")
        group.add_argument(
            "-d",
            "--double-encode",
            action='store_true', default=False,
            help="Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N). (no)")
        group.add_argument(
            "-t",
            "--trim-primer",
            action='store_true', default=False,
            help="Trim primer base and the first color (which is the transition "
                "to the first nucleotide). (no)")
        group.add_argument(
            "--strip-f3",
            action='store_true', default=False,
            help="Strip the _F3 suffix of read names. (no)")
        group.add_argument(
            "--maq", "--bwa",
            action='store_true', default=False,
            help="MAQ- and BWA-compatible colorspace output. This enables -c, -d, "
                "-t, --strip-f3 and -y '/1'. (no)")
        group.add_argument(
            "--no-zero-cap",
            dest='zero_cap', action='store_false',
            help="Do not change negative quality values to zero in colorspace "
                "data. By default, they are since many tools have problems with "
                "negative qualities. (no)")
        group.add_argument(
            "-z",
            "--zero-cap",
            action='store_true',
            help="Change negative quality values to zero. This is enabled "
                "by default when -c/--colorspace is also enabled. Use the above option "
                "to disable it. (no)")
        
        group = parser.add_argument_group("Paired-end options", description="The "
            "-A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but "
            "are applied to the second read in each pair.")
        group.add_argument(
            "-A",
            action='append', dest='adapters2', default=[], metavar='ADAPTER',
            help="3' adapter to be removed from second read in a pair. (no)")
        group.add_argument(
            "-G",
            action='append', dest='front2', default=[], metavar='ADAPTER',
            help="5' adapter to be removed from second read in a pair. (no)")
        group.add_argument(
            "-B",
            action='append', dest='anywhere2',  default=[], metavar='ADAPTER',
            help="5'/3 adapter to be removed from second read in a pair. (no)")
        group.add_argument(
            "-U",
            type=int, action='append', dest='cut2', default=[], metavar="LENGTH",
            help="Remove LENGTH bases from second read in a pair (see --cut). (no)")
        group.add_argument(
            "-I",
            "--cut-min2",
            type=int, action='append', default=[], metavar="LENGTH",
            help="Similar to -U, except that cutting is done AFTER adapter trimming, "
                 "and only if a minimum of LENGTH bases was not already removed "
                 "(see --cut-min). (no)")
        group.add_argument(
            "-p",
            "--paired-output",
            type=writeable_file, metavar="FILE",
            help="Write second read in a pair to FILE. (no)")
        group.add_argument(
            "-L",
            "--interleaved-output",
            type=writeable_file, metavar="FILE",
            help="Write output to interleaved file.")
        # Setting the default for pair_filter to None allows us to find out whether
        # the option was used at all.
        group.add_argument(
            "--pair-filter",
            choices=("any", "both"), default=None, metavar='(any|both)',
            help="Which of the reads in a paired-end read have to match the "
                "filtering criterion in order for it to be filtered. (any)")
        group.add_argument(
            "--untrimmed-paired-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write second read in a pair to this FILE when no adapter "
                "was found in the first read. Use this option together with "
                "--untrimmed-output when trimming paired-end reads. (no - output "
                "to same file as trimmed reads)")
        group.add_argument(
            "--too-short-paired-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write second read in a pair to this file if pair is too short. "
                "Use together with --too-short-output. (no - too short reads are "
                "discarded)")
        group.add_argument(
            "--too-long-paired-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write second read in a pair to this file if pair is too long. "
                "Use together with --too-long-output. (no - too long reads are "
                "discarded)")
        
        group = parser.add_argument_group("Method-specific options")
        group = group.add_mutually_exclusive_group()
        group.add_argument(
            "--bisulfite",
            default=False, metavar="METHOD",
            help="Set default option values for bisulfite-treated data. The argument specifies the "
                 "type of bisulfite library (rrbs, non-directional, non-directional-rrbs, truseq, "
                 "epignome, or swift) or custom parameters for trimming: '<read1>[;<read2>]' where "
                 "trimming parameters for each read are: '<5' cut>,<3' cut>,<include trimmed>,<only trimmed>' "
                 "where 'include trimmed' is 1 or 0 for whether or not the bases already trimmed "
                 "during/prior to adapter trimming should be counted towards the total bases to be "
                 "cut and 'only trimmed' is 1 or 0 for whether or not only trimmed reads should be "
                 "further cut. (no)")
        group.add_argument(
            "--mirna",
            action="store_true", default=False,
            help="Set default option values for miRNA data. (no)")
        
        group = parser.add_argument_group("Parallel (multi-core) options")
        group.add_argument(
            "-T",
            "--threads",
            type=positive(int, True), default=None, metavar="THREADS",
            help="Number of threads to use for read trimming. Set to 0 to use "
                 "max available threads. (Do not use multithreading)")
        group.add_argument(
            "--no-writer-process",
            action="store_false", dest="writer_process", default=True,
            help="Do not use a writer process; instead, each worker thread writes its "
                 "own output to a file with a '.N' suffix. (no)")
        group.add_argument(
            "--preserve-order",
            action="store_true", default=False,
            help="Preserve order of reads in input files (ignored if "
                 "--no-writer-process is set). (no)")
        group.add_argument(
            "--process-timeout",
            type=positive(int, True), default=60, metavar="SECONDS",
            help="Number of seconds process should wait before escalating messages "
                 "to ERROR level. (60)")
        group.add_argument(
            "--read-queue-size",
            type=int_or_str, default=None, metavar="SIZE",
            help="Size of queue for batches of reads to be processed. (THREADS * 100)")
        group.add_argument(
            "--result-queue-size",
            type=int_or_str, default=None, metavar="SIZE",
            help="Size of queue for batches of results to be written. (THREADS * 100)")
        group.add_argument(
            "--compression",
            choices=("worker", "writer"), default=None,
            help="Where data compression should be performed. Defaults to 'writer' if "
                 "system-level compression can be used and (1 < threads < 8), otherwise "
                 "defaults to 'worker'.")
    
    def validate_command_options(self):
        options = self.options
        parser = self.parser
        paired = options.paired
        
        if not paired:
            if not options.output:
                parser.error("An output file is required")
            if options.untrimmed_paired_output:
                parser.error("Option --untrimmed-paired-output can only be used when "
                    "trimming paired-end reads (with option -p).")
        else:
            if not options.interleaved_output:
                if not options.output:
                    parser.error("When you use -p or --paired-output, you must also "
                        "use the -o option.")
                if not options.paired_output:
                    parser.error("When paired-end trimming is enabled via -A/-G/-B/-U, "
                        "a second output file needs to be specified via -p (--paired-output).")
                if bool(options.untrimmed_output) != bool(options.untrimmed_paired_output):
                    parser.error("When trimming paired-end reads, you must use either none "
                        "or both of the --untrimmed-output/--untrimmed-paired-output options.")
                if options.too_short_output and not options.too_short_paired_output:
                    parser.error("When using --too-short-output with paired-end "
                        "reads, you also need to use --too-short-paired-output")
                if options.too_long_output and not options.too_long_paired_output:
                    parser.error("When using --too-long-output with paired-end "
                        "reads, you also need to use --too-long-paired-output")
            
            # Any of these options switch off legacy mode
            if (options.adapters2 or options.front2 or options.anywhere2 or options.cut2 or
                options.cut_min2 or options.interleaved_input or options.pair_filter or
                options.too_short_paired_output or options.too_long_paired_output):
                # Full paired-end trimming when both -p and -A/-G/-B/-U given
                # Read modifications (such as quality trimming) are applied also to second read.
                paired = 'both'
            else:
                # Modify first read only, keep second in sync (-p given, but not -A/-G/-B/-U).
                # This exists for backwards compatibility ('legacy mode').
                paired = 'first'
            
            options.paired = paired
        
        # If the user specifies a max rmp, that is used for determining the
        # minimum overlap and -O is set to 1, otherwise -O is set to the old
        # default of 3.
        # TODO: This is pretty confusing logic - need to simplify
        
        if options.aligner == 'adapter':
            if options.indels and options.indel_cost is None:
                options.indel_cost = 1
            if options.overlap is None:
                if options.adapter_max_rmp is None:
                    options.overlap = 3
                else:
                    options.overlap = 1
        elif options.aligner == 'insert':
            if paired != 'both':
                parser.error("Insert aligner only works with paired-end reads")
                # TODO: should also be checking that there is exactly one 3' adapter for each read
                # TODO: have the aligner tell us whether it can be used based on options?
            if options.indels and options.indel_cost is None:
                options.indel_cost = 3
            if options.overlap is None:
                options.overlap = 1
                if options.adapter_max_rmp is None:
                    options.adapter_max_rmp = 1E-6
            if options.insert_match_error_rate is None:
                options.insert_match_error_rate = options.error_rate or 0.2
            if options.insert_match_adapter_error_rate is None:
                options.insert_match_adapter_error_rate = options.insert_match_error_rate
                
        if options.merge_overlapping:
            if options.merged_output is None:
                parser.warn("--merge-output is not set; merged reads will be discarded")
            if options.merge_error_rate is None:
                options.merge_error_rate = options.error_rate or 0.2
        
        if options.mirna:
            if options.adapter is None and options.front is None and options.anywhere is None:
                options.adapter = ['TGGAATTCTCGG'] # illumina small RNA adapter
            if options.quality_cutoff is None:
                options.quality_cutoff = (20, 20)
            if options.minimum_length is None:
                options.minimum_length = 16
            if options.error_rate is None:
                options.error_rate = 0.12
        elif options.bisulfite:
            # TODO: set default adapter sequences
            # Jury is out on whether quality trimming helps. For aligners like
            # bwameth, it actually leads to worse results.
            #if options.quality_cutoff is None:
            #    options.quality_cutoff = "20,20"
            if options.bisulfite == "swift" and paired != "both":
                parser.error("Swift trimming is only compatible with paired-end reads")
            if options.bisulfite not in ("rrbs", "non-directional", "truseq", "epignome", "swift"):
                def parse_bisulfite_params(r):
                    try:
                        parts = [int(p) for p in r.split(",")]
                        assert len(parts) == 4
                        if parts[0] <= 0 and parts[1] <= 0:
                            return None
                        return dict(zip(
                            ("lengths", "count_trimmed", "only_trimmed"),
                            ((parts[0], -1 * parts[1]), (False,True)[parts[2]], (False,True)[parts[3]])
                        ))
                    except:
                        parser.error("Invalidate format for bisulfite parameters")
                
                temp = [parse_bisulfite_params(r) for r in options.bisulfite.split(";")]
                if paired == "both" and len(temp) == 1:
                    temp = [temp[0], temp[0]]
                elif paired != "both" and len(temp) > 1:
                    parser.error("Too many bisulfite parameters for single-end reads")
                options.bisulfite = temp
        
        if options.quality_cutoff:
            if all(c <= 0 for c in options.quality_cutoff):
                options.quality_cutoff = None
            elif len(options.quality_cutoff) == 1:
                options.quality_cutoff = [0] + options.quality_cutoff
        
        if options.pair_filter is None:
            options.pair_filter = 'any'

        if (options.discard_trimmed or options.discard_untrimmed) and (options.untrimmed_output is not None):
            parser.error("Only one of the --discard-trimmed, --discard-untrimmed "
                "and --untrimmed-output options can be used at the same time.")
        
        if options.output is not None and '{name}' in options.output:
            if options.discard_trimmed:
                parser.error("Do not use --discard-trimmed when demultiplexing.")
            if paired:
                parser.error("Demultiplexing not supported for paired-end files, yet.")
        
        if options.maq:
            options.colorspace = True
            options.double_encode = True
            options.trim_primer = True
            options.suffix = "/1"
        
        if options.strip_f3 or options.maq:
            options.strip_suffix.append('_F3')
        
        if options.zero_cap is None:
            options.zero_cap = options.colorspace
        
        if options.colorspace:
            if options.anywhere:
                parser.error("Using --anywhere with colorspace reads is currently not supported (if you "
                    "think this may be useful, contact the author).")
            if options.match_read_wildcards:
                parser.error('IUPAC wildcards not supported in colorspace')
            options.match_adapter_wildcards = False
        else:
            if options.trim_primer:
                parser.error("Trimming the primer makes only sense in colorspace.")
            if options.double_encode:
                parser.error("Double-encoding makes only sense in colorspace.")
                
        if options.error_rate is None:
            options.error_rate = 0.1
        
        if options.cut:
            if len(options.cut) > 2:
                parser.error("You cannot remove bases from more than two ends.")
            if len(options.cut) == 2 and options.cut[0] * options.cut[1] > 0:
                parser.error("You cannot remove bases from the same end twice.")
        
        if options.cut_min:
            if len(options.cut_min) > 2:
                parser.error("You cannot remove bases from more than two ends.")
            if len(options.cut_min) == 2 and options.cut_min[0] * options.cut_min[1] > 0:
                parser.error("You cannot remove bases from the same end twice.")
        
        if paired == 'both' and options.cut2:
            if len(options.cut2) > 2:
                parser.error("You cannot remove bases from more than two ends.")
            if len(options.cut2) == 2 and options.cut2[0] * options.cut2[1] > 0:
                parser.error("You cannot remove bases from the same end twice.")

        if paired == 'both' and options.cut_min2:
            if len(options.cut_min2) > 2:
                parser.error("You cannot remove bases from more than two ends.")
            if len(options.cut_min2) == 2 and options.cut_min2[0] * options.cut_min2[1] > 0:
                parser.error("You cannot remove bases from the same end twice.")
        
        if options.threads is not None:
            if options.debug:
                parser.error("Cannot use debug mode with multiple threads")
            
            threads = options.threads
            if threads <= 0:
                threads = cpu_count()
            elif threads == 1:
                parser.error("--threads must be >= 2")
            options.threads = threads
            
            if options.compression is None:
                # Our tests show that with 8 or more threads, worker compression is
                # more efficient.
                if options.writer_process and 2 < threads < 8:
                    from atropos.compression import can_use_system_compression
                    options.compression = "writer" if can_use_system_compression() else "worker"
                else:
                    options.compression = "worker"
            elif options.compression == "writer":
                if not options.writer_process:
                    parser.error("Writer compression and --no-writer-process are mutually exclusive")
                elif threads == 2:
                    logging.getLogger.warn("Writer compression requires > 2 threads; using worker compression instead")
                    options.compression = "worker"
            
            # Set queue sizes if necessary.
            # If we are using writer compression, the back-up will be in the result queue,
            # otherwise it will be in the read queue.
            if options.read_queue_size is None:
                options.read_queue_size = threads * (100 if options.compression == "writer" else 1000)
            elif options.read_queue_size > 0:
                assert options.read_queue_size >= threads
        
            if options.result_queue_size is None:
                options.result_queue_size = threads * (100 if options.compression == "worker" else 1000)
            elif options.result_queue_size > 0:
                assert options.result_queue_size >= threads

COMMANDS['trim'] = TrimCommand

class DetectCommand(Command):
    """detect adapter and other contaminant sequences."""
    name = "detect"
    usage = "atropos -se input.fastq detect\natropos -pe1 in1.fq -pe2 in2.fq detect"
    description = "Detect adapter sequences directly from read sequences."
    
    def add_command_options(self):
        parser = self.parser
        parser.set_defaults(max_reads=10000)
        parser.add_argument(
            "-d",
            "--detector",
            choices=('known', 'heuristic', 'khmer'), default=None,
            help="Which detector to use. (automatically choose based on other options)")
        parser.add_argument(
            "-o",
            "--output",
            type=writeable_file, default=STDOUT, metavar="FILE",
            help="File in which to write the summary of detected adapters. (stdout)")
        parser.add_argument(
            "-k",
            "--kmer-size",
            type=positive(), default=12,
            help="Size of k-mer used to scan reads for adapter sequences. (12)")
        parser.add_argument(
            "-m",
            "--max-adapters",
            type=positive(), default=None,
            help="The maximum number of candidate adapters to report. (report all)")
        
        group = parser.add_argument_group("Contaminants")
        group.add_argument(
            "-i",
            "--include-contaminants",
            choices=('all','known','unknown'), default='all',
            help="What conaminants to search for: all, only known adapters/contaminants ('known'), "
                 "or only unknown contaminants ('unknown'). (all)")
        group.add_argument(
            "-x",
            "--known-contaminant",
            action="append", dest='known_adapter', default=None,
            help="Pass known contaminants in on the commandline as 'name=sequence'. "
                 "Can be specified multiple times.")
        group.add_argument(
            "-F",
            "--known-contaminants-file",
            type=readable_url, action="append", dest='known_adapters_file', default=None,
            help="Points to FASTA File or URL with known contaminants.")
        group.add_argument(
            "--no-default-contaminants",
            action="store_false", dest="default_adapters", default=True,
            help="Don't fetch the default contaminant list (which is currently stored as a "
                "GitHub gist).")
        group.add_argument(
            "--contaminant-cache-file",
            type=readwriteable_file, dest='adapter_cache_file', default='.adapters',
            help="File where known contaminant sequences will be cached, unless "
                "--no-cache-contaminants is set.")
        group.add_argument(
            "--no-cache-contaminants",
            action="store_false", dest="cache_adapters", default=True,
            help="Don't cache contaminant list as '.contaminants' in the working directory.")

COMMANDS['detect'] = DetectCommand

class ErrorCommand(Command):
    """estimate the sequencing error rate."""
    name = "error"
    usage = "atropos -se input.fastq error\natropos -pe in1.fq -pe2 in2.fq error"
    description = """
Estimate the error rate from base qualities. This can help to
determine the quality of your data, and to decide the value for
the max error rate (-e) parameter. Normal error for an Illumina
experiment is around 1% or less. We recommend setting -e to
10 x the empirical error rate."""
    
    def add_command_options(self):
        parser = self.parser
        parser.set_defaults(max_reads=10000)
        parser.add_argument(
            "-o",
            "--output",
            type=writeable_file, default=STDOUT,
            help="File in which to write the summary of the estimated error rates. (stdout)")

COMMANDS['error'] = ErrorCommand

# Main

def main(cmdlineargs=None):
    """
    Main function that evaluates command-line parameters.
    
    :param cmdlineargs: iterable of command line arguments; `None` is equivalent to
    `sys.argv[1:]`
    """
    args = cmdlineargs or sys.argv[1:]
    if len(args) == 0 or any(h in args for h in ('-h', '--help')):
        print_subcommands()
        return
    
    if args[0][0] == '-':
        command_name = "trim"
    else:
        command_name = args[0]
        del args[0]
    
    command_class = COMMANDS[command_name]
    command = command_class(args)
    command.execute()

def print_subcommands():
    print("Atropos version {}\n".format(__version__))
    print("usage: atropos <command> [options]\n")
    print("commands:")
    for command_class in COMMANDS.values():
        print("  {}: {}".format(command_class.name, command_class.__doc__))
    print("\noptional arguments:\n  -h, --help  show this help message and exit")
    print("""
Use "atropos <command> --help" to see all options for a specific command.
See http://atropos.readthedocs.org/ for full documentation.

Atropos is a fork of Cutadapt 1.10 (
https://github.com/marcelm/cutadapt/tree/2f3cc0717aa9ff1e0326ea6bcb36b712950d4999)
by John Didion, "Atropos: sensitive, specific, and speedy trimming of NGS reads,"
in prep.

Cutadapt (https://github.com/marcelm/cutadapt) was developed by Marcel Martin,
"Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads,"
EMBnet Journal, 2011, 17(1):10-12.
""")

if __name__ == '__main__':
    main()
