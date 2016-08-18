#!/usr/bin/env python
# -*- coding: utf-8 -*-
# kate: word-wrap off; remove-trailing-spaces all;

"""
atropos version {version}

Atropos removes adapter sequences from high-throughput sequencing reads.

Usage:
    atropos -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    atropos -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. The reverse complement is *not* automatically
searched. All reads from input.fastq will be written to output.fastq with the
adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter
sequences can be given (use further -a options), but only the best-matching
adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

Use "atropos --help" to see all command-line options.
See http://atropos.readthedocs.org/ for full documentation.

Cutadapt (https://github.com/marcelm/cutadapt) was developed by Marcel Martin,
"Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads,"
EMBnet Journal, 2011, 17(1):10-12.

Atropos is a fork of Cutadapt 1.10 (
https://github.com/marcelm/cutadapt/tree/2f3cc0717aa9ff1e0326ea6bcb36b712950d4999)
by John Didion, "Atropos: sensitive, specific, and speedy trimming of NGS reads,"
in prep.
"""

# Print a helpful error message if the extension modules cannot be imported.
from atropos import *
check_importability()

from atropos import __version__
from atropos.seqio import (open_reader, UnknownFileType, Formatters, RestFormatter,
                            InfoFormatter, WildcardFormatter, FormatError, Writers)
from atropos.adapters import AdapterParser, BACK
from atropos.modifiers import *
from atropos.filters import *
from atropos.report import *

import argparse
from collections import defaultdict
import copy
import logging
import os
import platform
import sys
import time

def main(cmdlineargs=None, default_outfile="-"):
    """
    Main function that evaluates command-line parameters and iterates over all reads.
    
    :param cmdlineargs: iterable of command line arguments; `None` is equivalent to
    `sys.argv[1:]`
    :param default_outfile: the file to which trimmed reads are sent if the ``-o``
    parameter is not used.
    """
    orig_args = copy.copy(cmdlineargs or sys.argv)
    parser = get_argument_parser()
    options = parser.parse_args(args=cmdlineargs)
    
    # Setup logging only if there are not already any handlers (can happen when
    # this function is being called externally such as from unit tests)
    if not logging.root.handlers:
        level = options.log_level or ("ERROR" if options.quiet else "INFO")
        setup_logging(stdout=bool(options.output), level=level)
    
    paired = validate_options(options, parser)
    reader, qualities, has_qual_file = create_reader(options.inputs, options, paired, parser)
    modifiers, num_adapters = create_modifiers(options, paired, qualities, has_qual_file, parser)
    min_affected = 2 if options.pair_filter == 'both' else 1
    filters = create_filters(options, paired, min_affected)
    formatters, force_create = create_formatters(options, qualities, default_outfile)
    writers = Writers(force_create)
    
    logger = logging.getLogger()
    logger.info("This is atropos %s with Python %s", __version__, platform.python_version())
    logger.info("Command line parameters: %s", " ".join(orig_args))
    logger.info("Trimming %s adapter%s with at most %.1f%% errors in %s mode ...",
        num_adapters, 's' if num_adapters > 1 else '', options.error_rate * 100,
        { False: 'single-end', 'first': 'paired-end legacy', 'both': 'paired-end' }[paired])
    if paired == 'first' and (len(modifiers.get_modifiers(read=2)) > 0 or options.quality_cutoff):
        import textwrap
        logger.warning('\n'.join(textwrap.wrap('WARNING: Requested read '
            'modifications are applied only to the first '
            'read since backwards compatibility mode is enabled. '
            'To modify both reads, also use any of the -A/-B/-G/-U options. '
            'Use a dummy adapter sequence when necessary: -A XXX')))
    
    start_wallclock_time = time.time()
    start_cpu_time = time.clock()
    
    if options.threads is None:
        # Run single-threaded version
        rc, summary = run_serial(reader, modifiers, filters, formatters, writers)
    else:
        # Run multiprocessing version
        import atropos.multicore
        rc, summary = atropos.multicore.run_parallel(
            reader, modifiers, filters, formatters, writers, options.threads,
            options.process_timeout, options.preserve_order, options.read_queue_size,
            options.result_queue_size, not options.no_writer_process, options.compression)
    
    if rc != 0:
        sys.exit(rc)
    
    stop_wallclock_time = time.time()
    stop_cpu_time = time.clock()
    report = print_report(
        paired,
        options,
        stop_wallclock_time - start_wallclock_time,
        stop_cpu_time - start_cpu_time,
        summary,
        modifiers.get_trimmer_classes())
    
def get_argument_parser():
    parser = argparse.ArgumentParser(usage=__doc__.lstrip().format(version=__version__))
    parser.add_argument("inputs", nargs=argparse.REMAINDER)
    
    parser.add_argument("--debug", action='store_true', default=False,
        help="Print debugging information. (no)")
    parser.add_argument("-f", "--format",
        help="Input file format; can be either 'fasta', 'fastq' or 'sra-fastq'. "
            "Ignored when reading csfasta/qual files. (auto-detect "
            "from file name extension)")
    parser.add_argument("--max-reads", default=None,
        help="Maximum number of reads/pairs to process (no max)")
    parser.add_argument("--progress", default=None,
        help="Show progress. bar = show progress bar; msg = show a status "
             "message. (no)")
    parser.add_argument("--op-order", default="CGQA",
        help="The order in which trimming operations are be applied. This is a string of "
             "1-4 of the following characters: A = adapter trimming; C = cutting "
             "(unconditional); G = NextSeq trimming; Q = quality trimming. The default is "
             "'CGQA' to maintain compatibility with Cutadapt; however, this is likely to "
             "change to 'GACQ' in the near future.")
    
    group = parser.add_argument_group("Finding adapters:",
        description="Parameters -a, -g, -b specify adapters to be removed from "
            "each read (or from the first read in a pair if data is paired). "
            "If specified multiple times, only the best matching adapter is "
            "trimmed (but see the --times option). When the special notation "
            "'file:FILE' is used, adapter sequences are read from the given "
            "FASTA file.")
    group.add_argument("-a", "--adapter", action="append", default=[], metavar="ADAPTER",
        dest="adapters",
        help="Sequence of an adapter ligated to the 3' end (paired data: of the "
            "first read). The adapter and subsequent bases are trimmed. If a "
            "'$' character is appended ('anchoring'), the adapter is only "
            "found if it is a suffix of the read. (none)")
    group.add_argument("-g", "--front", action="append", default=[], metavar="ADAPTER",
        help="Sequence of an adapter ligated to the 5' end (paired data: of the "
            "first read). The adapter and any preceding bases are trimmed. "
            "Partial matches at the 5' end are allowed. If a '^' character is "
            "prepended ('anchoring'), the adapter is only found if it is a "
            "prefix of the read. (none)")
    group.add_argument("-b", "--anywhere", action="append", default=[], metavar="ADAPTER",
        help="Sequence of an adapter that may be ligated to the 5' or 3' end "
            "(paired data: of the first read). Both types of matches as "
            "described under -a und -g are allowed. If the first base of the "
            "read is part of the match, the behavior is as with -g, otherwise "
            "as with -a. This option is mostly for rescuing failed library "
            "preparations - do not use if you know which end your adapter was "
            "ligated to! (none)")
    group.add_argument("--aligner", choices=('adapter', 'insert'), default='adapter',
        help="Which alignment algorithm to use for identifying adapters. Currently, "
             "you can choose between the semi-global alignment strategy used in Cutdapt "
             "('adapter') or the more accurate insert-based alignment algorithm ('insert'). "
             "Note that insert-based alignment can only be used with paired-end reads "
             "containing 3' adapters. New algorithms are being implemented and the default "
             "is likely to change. (adapter)")
    group.add_argument("-e", "--error-rate", type=float, default=None,
        help="Maximum allowed error rate (no. of errors divided by the length "
            "of the matching region). (0.1)")
    group.add_argument("--no-indels", action='store_false', dest='indels', default=True,
        help="Allow only mismatches in alignments. "
            "(allow both mismatches and indels)")
    group.add_argument("-n", "--times", type=int, metavar="COUNT", default=1,
        help="Remove up to COUNT adapters from each read. (1)")
    group.add_argument("-O", "--overlap", type=int, metavar="MINLENGTH", default=3,
        help="If the overlap between the read and the adapter is shorter than "
            "MINLENGTH, the read is not modified. Reduces the no. of bases "
            "trimmed due to random adapter matches. (3)")
    group.add_argument("--match-read-wildcards", action="store_true", default=False,
        help="Interpret IUPAC wildcards in reads. (no)")
    group.add_argument("-N", "--no-match-adapter-wildcards", action="store_false",
        default=True, dest='match_adapter_wildcards',
        help="Do not interpret IUPAC wildcards in adapters. (no)")
    group.add_argument("--no-trim", dest='action', action='store_const', const=None,
        help="Match and redirect reads to output/untrimmed-output as usual, "
            "but do not remove adapters. (no)")
    group.add_argument("--mask-adapter", dest='action', action='store_const', const='mask',
        help="Mask adapters with 'N' characters instead of trimming them. (no)")
    group.add_argument("-R", "--merge-overlapping", action="store_true", default=False,
        help="Merge read pairs that overlap into a single sequence. This is an experimental "
             "and currently only works with read pairs where an adapter match was found at "
             "the 3' ends of both reads. The minimum required overlap is controled by "
             "--merge-min-overlap, and the mimimum error rate in the alignment is controlled "
             "by --error-rate. (no)")
    group.add_argument("--merge-min-overlap", type=float, default=0.9,
        help="The minimum overlap between reads required for merging. If this number is (0,1.0], "
             "it specifies the minimum length as the fraction of the length of the *shorter* read "
             "in the pair; otherwise it specifies the minimum number of overlapping base pairs ("
             "with an absolute minimum of 2 bp). (0.9)")
    # TODO: error correction:
    #group.add_argument("--merge-mismatches", choices=("best", "N", "random"), default="best",
    #	help="How to handle mismatches while merging; best: select the base with the best quality; "
    #		 "N: set the base to N; random: select one of the two bases at random. Note that if "
    #		 "exactly one base is ambiguous, the non-ambiguous base is always used.")
    
    group = parser.add_argument_group("Additional read modifications")
    group.add_argument("-u", "--cut", action='append', default=[], type=int, metavar="LENGTH",
        help="Remove bases from each read (first read only if paired). "
            "If LENGTH is positive, remove bases from the beginning. "
            "If LENGTH is negative, remove bases from the end. "
            "Can be used twice if LENGTHs have different signs. (no)")
    group.add_argument("-q", "--quality-cutoff", default=None, metavar="[5'CUTOFF,]3'CUTOFF",
        help="Trim low-quality bases from 5' and/or 3' ends of each read before "
            "adapter removal. Applied to both reads if data is paired. If one "
            "value is given, only the 3' end is trimmed. If two "
            "comma-separated cutoffs are given, the 5' end is trimmed with "
            "the first cutoff, the 3' end with the second. (no)")
    group.add_argument("-i", "--cut-min", action='append', default=[], type=int, metavar="LENGTH",
        help="Similar to -u, except that cutting is done AFTER adapter trimming, and only "
            "if a minimum of LENGTH bases was not already removed. (no)")
    group.add_argument("--nextseq-trim", type=int, default=None, metavar="3'CUTOFF",
        help="NextSeq-specific quality trimming (each read). Trims also dark "
            "cycles appearing as high-quality G bases (EXPERIMENTAL). (no)")
    group.add_argument("--quality-base", type=int, default=33,
        help="Assume that quality values in FASTQ are encoded as ascii(quality "
            "+ QUALITY_BASE). This needs to be set to 64 for some old Illumina "
            "FASTQ files. (33)")
    group.add_argument("--trim-n", action='store_true', default=False,
        help="Trim N's on ends of reads. (no)")
    group.add_argument("-x", "--prefix", default='',
        help="Add this prefix to read names. Use {name} to insert the name of "
             "the matching adapter. (no)")
    group.add_argument("-y", "--suffix", default='',
        help="Add this suffix to read names; can also include {name}. (no)")
    group.add_argument("--strip-suffix", action='append', default=[],
        help="Remove this suffix from read names if present. Can be given "
             "multiple times. (no)")
    group.add_argument("--length-tag", metavar="TAG",
        help="Search for TAG followed by a decimal number in the description "
            "field of the read. Replace the decimal number with the correct "
            "length of the trimmed read. For example, use --length-tag 'length=' "
            "to correct fields like 'length=123'. (no)")
    
    group = parser.add_argument_group("Filtering of processed reads")
    group.add_argument("--discard-trimmed", "--discard", action='store_true', default=False,
        help="Discard reads that contain an adapter. Also use -O to avoid "
            "discarding too many randomly matching reads! (no)")
    group.add_argument("--discard-untrimmed", "--trimmed-only", action='store_true', default=False,
        help="Discard reads that do not contain the adapter. (no)")
    group.add_argument("-m", "--minimum-length", type=int, default=None, metavar="LENGTH",
        help="Discard trimmed reads that are shorter than LENGTH. Reads that "
            "are too short even before adapter removal are also discarded. In "
            "colorspace, an initial primer is not counted. (0)")
    group.add_argument("-M", "--maximum-length", type=int, default=sys.maxsize, metavar="LENGTH",
        help="Discard trimmed reads that are longer than LENGTH. "
            "Reads that are too long even before adapter removal "
            "are also discarded. In colorspace, an initial primer "
            "is not counted. (no limit)")
    group.add_argument("--max-n", type=float, default=-1.0, metavar="COUNT",
        help="Discard reads with too many N bases. If COUNT is an integer, it "
            "is treated as the absolute number of N bases. If it is between 0 "
            "and 1, it is treated as the proportion of N's allowed in a read. "
            "(no)")
    
    group = parser.add_argument_group("Output")
    group.add_argument("--quiet", default=False, action='store_true',
        help="Print only error messages. (no)")
    group.add_argument("-o", "--output", metavar="FILE",
        help="Write trimmed reads to FILE. FASTQ or FASTA format is chosen "
            "depending on input. The summary report is sent to standard output. "
            "Use '{name}' in FILE to demultiplex reads into multiple "
            "files. (write to standard output)")
    group.add_argument("--info-file", metavar="FILE",
        help="Write information about each read and its adapter matches into FILE. "
             "See the documentation for the file format. (no)")
    group.add_argument("-r", "--rest-file", metavar="FILE",
        help="When the adapter matches in the middle of a read, write the "
            "rest (after the adapter) into FILE. (no)")
    group.add_argument("--wildcard-file", metavar="FILE",
        help="When the adapter has N bases (wildcards), write adapter bases "
            "matching wildcard positions to FILE. When there are indels in the "
            "alignment, this will often not be accurate. (no)")
    group.add_argument("--too-short-output", metavar="FILE",
        help="Write reads that are too short (according to length specified by "
        "-m) to FILE. (no - too short reads are discarded)")
    group.add_argument("--too-long-output", metavar="FILE",
        help="Write reads that are too long (according to length specified by "
        "-M) to FILE. (no - too long reads are discarded)")
    group.add_argument("--untrimmed-output", default=None, metavar="FILE",
        help="Write reads that do not contain the adapter to FILE. (no - "
             "untrimmed reads are written to default output)")
    group.add_argument("--merged-output", default=None, metavar="FILE",
        help="Write reads that have been merged to this file. By default, merged "
             "reads are written to --too-short-output if specified, and otherwise "
             "are discarded. (no - merged reads are written to default output)")
    group.add_argument("--report-file", default=None, metavar="FILE",
        help="Write report to file rather than stdout/stderr. (no)")
    group.add_argument("--log-level", default=None,
        choices=("DEBUG", "INFO", "WARN", "ERROR"),
        help="Logging level. (ERROR when --quiet else INFO)")
    
    group = parser.add_argument_group("Colorspace options")
    group.add_argument("-c", "--colorspace", action='store_true', default=False,
        help="Enable colorspace mode: Also trim the color that is adjacent to "
             "the found adapter. (no)")
    group.add_argument("-d", "--double-encode", action='store_true', default=False,
        help="Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N). (no)")
    group.add_argument("-t", "--trim-primer", action='store_true', default=False,
        help="Trim primer base and the first color (which is the transition "
            "to the first nucleotide). (no)")
    group.add_argument("--strip-f3", action='store_true', default=False,
        help="Strip the _F3 suffix of read names. (no)")
    group.add_argument("--maq", "--bwa", action='store_true', default=False,
        help="MAQ- and BWA-compatible colorspace output. This enables -c, -d, "
            "-t, --strip-f3 and -y '/1'. (no)")
    group.add_argument("--no-zero-cap", dest='zero_cap', action='store_false',
        help="Do not change negative quality values to zero in colorspace "
            "data. By default, they are since many tools have problems with "
            "negative qualities. (no)")
    group.add_argument("--zero-cap", "-z", action='store_true',
        help="Change negative quality values to zero. This is enabled "
        "by default when -c/--colorspace is also enabled. Use the above option "
        "to disable it. (no)")
    parser.set_defaults(zero_cap=None, action='trim')
    
    group = parser.add_argument_group("Paired-end options", description="The "
        "-A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but "
        "are applied to the second read in each pair.")
    group.add_argument("-A", dest='adapters2', action='append', default=[],
        metavar='ADAPTER',
        help="3' adapter to be removed from second read in a pair. (no)")
    group.add_argument("-G", dest='front2', action='append', default=[],
        metavar='ADAPTER',
        help="5' adapter to be removed from second read in a pair. (no)")
    group.add_argument("-B", dest='anywhere2', action='append', default=[],
        metavar='ADAPTER',
        help="5'/3 adapter to be removed from second read in a pair. (no)")
    group.add_argument("-U", dest='cut2', action='append', default=[], type=int,
        metavar="LENGTH",
        help="Remove LENGTH bases from second read in a pair (see --cut). (no)")
    group.add_argument("-I", "--cut-min2", action='append', default=[],
        type=int, metavar="LENGTH",
        help="Similar to -U, except that cutting is done AFTER adapter trimming, "
             "and only if a minimum of LENGTH bases was not already removed "
             "(see --cut-min). (no)")
    group.add_argument("-p", "--paired-output", metavar="FILE",
        help="Write second read in a pair to FILE. (no)")
    # Setting the default for pair_filter to None allows us to find out whether
    # the option was used at all.
    group.add_argument("--pair-filter", metavar='(any|both)', default=None,
        choices=("any", "both"),
        help="Which of the reads in a paired-end read have to match the "
            "filtering criterion in order for it to be filtered. (any)")
    group.add_argument("--interleaved", action='store_true', default=False,
        help="Read and write interleaved paired-end reads. (no)")
    group.add_argument("--untrimmed-paired-output", metavar="FILE",
        help="Write second read in a pair to this FILE when no adapter "
            "was found in the first read. Use this option together with "
            "--untrimmed-output when trimming paired-end reads. (no - output "
            "to same file as trimmed reads)")
    group.add_argument("--too-short-paired-output", metavar="FILE", default=None,
        help="Write second read in a pair to this file if pair is too short. "
            "Use together with --too-short-output. (no - too short reads are "
            "discarded)")
    group.add_argument("--too-long-paired-output", metavar="FILE", default=None,
        help="Write second read in a pair to this file if pair is too long. "
            "Use together with --too-long-output. (no - too long reads are "
            "discarded)")
    
    group = parser.add_argument_group("Method-specific options")
    group.add_argument("--bisulfite", default=False, metavar="METHOD",
        help="Set default option values for bisulfite-treated data. The argument specifies the "
             "type of bisulfite library (non-directional, truseq, epignome, or swift) or custom "
             "parameters for trimming: '<read1>[;<read2>]' where trimming parameters for each "
             "read are: '<5' cut>,<3' cut>,<include trimmed>,<only trimmed>' where 'include "
             "trimmed' is 1 or 0 for whether or not the bases already trimmed during/prior to adapter "
             "trimming should be counted towards the total bases to be cut and 'only trimmed' "
             "is 1 or 0 for whether or not only trimmed reads should be further cut. (no)")
    group.add_argument("--mirna", action="store_true", default=False,
        help="Set default option values for miRNA data. (no)")
    
    group = parser.add_argument_group("Parallel (multi-core) options")
    group.add_argument("-T", "--threads", type=int, default=None, metavar="THREADS",
        help="Number of threads to use for read trimming. Set to 0 to use "
             "max available threads. (Do not use multithreading)")
    group.add_argument("--no-writer-process", action="store_true", default=False,
        help="Do not use a writer process; instead, each worker thread writes its "
             "own output to a file with a '.N' suffix. (no)")
    group.add_argument("--preserve-order", action="store_true", default=False,
        help="Preserve order of reads in input files (ignored if "
             "--no-writer-process is set). (no)")
    group.add_argument("--process-timeout", type=int, default=60, metavar="SECONDS",
        help="Number of seconds process should wait before escalating messages "
             "to ERROR level. (60)")
    group.add_argument("--batch-size", type=int, default=5000, metavar="SIZE",
        help="Number of records to process in each batch. (5000)")
    group.add_argument("--read-queue-size", type=int, default=None, metavar="SIZE",
        help="Size of queue for batches of reads to be processed. (THREADS * 100)")
    group.add_argument("--result-queue-size", type=int, default=None, metavar="SIZE",
        help="Size of queue for batches of results to be written. (THREADS * 100)")
    group.add_argument("--compression", choices=("worker", "writer"), default=None,
        help="Where data compression should be performed. Defaults to 'writer' if "
             "system-level compression can be used and (1 < threads < 8), otherwise "
             "defaults to 'worker'.")
    
    return parser

def setup_logging(stdout=False, level="INFO"):
    """
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    level = getattr(logging, level)
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
    stream_handler.setLevel(level)
    logging.getLogger().setLevel(level)
    logging.getLogger().addHandler(stream_handler)

MAGNITUDE = dict(
    G=(1E9, "000000000"),
    M=(1E6, "000000"),
    K=(1E3, "000")
)

def validate_options(options, parser):
    if sum(int(opt is not False) for opt in (options.mirna, options.bisulfite)) > 1:
        parser.error("Cannot specify more than one method-specific option")
    
    if options.debug and options.threads is not None:
        parser.error("Cannot use debug mode with multiple threads")
    
    args = options.inputs
    if len(args) == 0:
        parser.error("At least one parameter needed: name of a FASTA or FASTQ file.")
    elif len(args) > 2:
        parser.error("Too many parameters: {}".format(args))
    if args[0].endswith('.qual'):
        parser.error("If a .qual file is given, it must be the second argument.")
    
    # Find out which 'mode' we need to use.
    # Default: single-read trimming (neither -p nor -A/-G/-B/-U/--interleaved given)
    paired = False
    if options.paired_output:
        # Modify first read only, keep second in sync (-p given, but not -A/-G/-B/-U).
        # This exists for backwards compatibility ('legacy mode').
        paired = 'first'
    # Any of these options switch off legacy mode
    if (options.adapters2 or options.front2 or options.anywhere2 or
        options.cut2 or options.cut_min2 or options.interleaved or options.pair_filter or
        options.too_short_paired_output or options.too_long_paired_output):
        # Full paired-end trimming when both -p and -A/-G/-B/-U given
        # Read modifications (such as quality trimming) are applied also to second read.
        paired = 'both'

    if paired and len(args) == 1 and not options.interleaved:
        parser.error("When paired-end trimming is enabled via -A/-G/-B/-U or -p, "
            "two input files are required.")
    if options.interleaved and len(args) != 1:
        parser.error("When reading interleaved files, only one input file may "
            "be given.")
    if not paired:
        if options.untrimmed_paired_output:
            parser.error("Option --untrimmed-paired-output can only be used when "
                "trimming paired-end reads (with option -p).")

    if paired:
        if not options.interleaved:
            if not options.paired_output:
                parser.error("When paired-end trimming is enabled via -A/-G/-B/-U, "
                    "a second output file needs to be specified via -p (--paired-output).")
            if not options.output:
                parser.error("When you use -p or --paired-output, you must also "
                    "use the -o option.")
            if bool(options.untrimmed_output) != bool(options.untrimmed_paired_output):
                parser.error("When trimming paired-end reads, you must use either none "
                    "or both of the --untrimmed-output/--untrimmed-paired-output options.")
            if options.too_short_output and not options.too_short_paired_output:
                parser.error("When using --too-short-output with paired-end "
                    "reads, you also need to use --too-short-paired-output")
            if options.too_long_output and not options.too_long_paired_output:
                parser.error("When using --too-long-output with paired-end "
                    "reads, you also need to use --too-long-paired-output")
    elif len(args) == 2:
        if options.format is not None:
            parser.error("If a pair of .fasta and .qual files is given, the -f/--format "
                "parameter cannot be used.")
    
    if options.aligner != 'adapter':
        if options.aligner == 'insert' and paired != 'both':
            parser.error("Insert aligner only works with paired-end reads")
            # TODO: should also be checking that there is exactly one 3' adapter for each read
            # TODO: have the aligner tell us whether it can be used based on options?
    
    if options.merge_overlapping and (
            paired != "both" or
            options.adapters is None or
            options.adapters2 is None or
            options.times != 1):
        parser.error("--merge-overlapping currently only works for paired-end reads with 3' adapters.")
    
    if options.format is not None and options.format.lower() not in ['fasta', 'fastq', 'sra-fastq']:
        parser.error("The input file format must be either 'fasta', 'fastq' or "
            "'sra-fastq' (not '{0}').".format(options.format))
    
    options.op_order = list(options.op_order)
    if not (1 <= len(options.op_order) <= 4 and all(o in ('A','C','G','Q') for o in options.op_order)):
        parser.error("--op-order must be a string of 1-4 characters consisting of A,C,G,Q")
    
    if options.mirna and options.bisulfite is not None:
        parser.error("Can only specify one of --mirna, --bisulfite")
    if options.mirna:
        if options.adapter is None and options.front is None and options.anywhere is None:
            options.adapter = ['TGGAATTCTCGG'] # illumina small RNA adapter
        if options.quality_cutoff is None:
            options.quality_cutoff = "20,20"
        if options.minimum_length is None:
            options.minimum_length = 16
        if options.error_rate is None:
            options.error_rate = 0.12
    elif options.bisulfite:
        if options.quality_cutoff is None:
            options.quality_cutoff = "20,20"
        if options.bisulfite not in ("non-directional", "truseq", "epignome", "swift"):
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
        elif not paired and options.bisulfite != "non-directional":
            parser.error("Single-end bisulfite libraries must be set to 'non-directional'.")
    
    if options.quality_cutoff is not None:
        cutoffs = options.quality_cutoff.split(',')
        if len(cutoffs) == 1:
            try:
                cutoffs_array = [0, int(cutoffs[0])]
            except ValueError as e:
                parser.error("Quality cutoff value not recognized: {0}".format(e))
        elif len(cutoffs) == 2:
            try:
                cutoffs_array = [int(cutoffs[0]), int(cutoffs[1])]
            except ValueError as e:
                parser.error("Quality cutoff value not recognized: {0}".format(e))
        else:
            parser.error("Expected one value or two values separated by comma for the quality cutoff")
        
        if all(x <= 0 for x in cutoffs_array):
            options.quality_cutoff = None
        else:
            options.quality_cutoff = cutoffs_array
        
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
    
    if options.trim_primer and not options.colorspace:
        parser.error("Trimming the primer makes only sense in colorspace.")
    
    if options.double_encode and not options.colorspace:
        parser.error("Double-encoding makes only sense in colorspace.")
    
    if options.anywhere and options.colorspace:
        parser.error("Using --anywhere with colorspace reads is currently not supported (if you "
            "think this may be useful, contact the author).")
    
    if options.error_rate is None:
        options.error_rate = 0.1
    elif not (0 <= options.error_rate <= 1.):
        parser.error("The maximum error rate must be between 0 and 1.")
    
    if options.overlap < 1:
        parser.error("The overlap must be at least 1.")

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
    
    if options.colorspace:
        if options.match_read_wildcards:
            parser.error('IUPAC wildcards not supported in colorspace')
        options.match_adapter_wildcards = False
    
    if options.quiet or options.output is None or options.output == "-":
        options.progress = None
    elif options.progress == "bar":
        try:
            import progressbar
        except:
            parser.error("The python-progressbar library is required for --progress=bar")
    
    # TODO: once we switch to argparse, int_or_str can be passed
    # as the argument type
    def int_or_str(x):
        if x is None or isinstance(x, int):
            return x
        elif isinstance(x, str):
            x = x.upper()
            for a, mag in MAGNITUDE.items():
                x = x.replace(a, mag[1])
            return int(x)
        else:
            raise Exception("Unsupported type {}".format(x))
    
    options.max_reads = int_or_str(options.max_reads)
    
    if options.threads is not None:
        threads = options.threads
        if threads <= 0:
            threads = cpu_count()
        if threads < 2:
            parser.error("At least two threads are required for multi-processing")
        options.threads = threads
        
        if options.compression is None:
            # Our tests show that with 8 or more threads, worker compression is
            # more efficient.
            if options.no_writer_process or threads == 2 or threads >= 8:
                options.compression = "worker"
            else:
                from atropos.compression import can_use_system_compression
                options.compression = "writer" if can_use_system_compression() else "worker"
        elif options.compression == "writer" and options.no_writer_process:
            parser.error("Writer compression and --no-writer-process are mutually exclusive")
        elif threads == 2 and options.compression == "writer":
            logging.getLogger.warn("Writer compression requires > 2 threads; using worker compression instead")
            options.compression = "worker"
        
        options.batch_size = int_or_str(options.batch_size)
        if options.process_timeout < 0:
            options.process_timeout = 0
        
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
    
    return paired

def create_reader(input_files, options, paired, parser, counter_magnitude="M"):
    input_filename = input_files[0]
    input_paired_filename = None
    quality_filename = None
    if len(input_files) > 1:
        if paired:
            input_paired_filename = input_files[1]
        else:
            quality_filename = input_files[1]
    
    try:
        reader = open_reader(
            input_filename, file2=input_paired_filename,
            qualfile=quality_filename, colorspace=options.colorspace,
            fileformat=options.format, interleaved=options.interleaved)
    except (UnknownFileType, IOError) as e:
        parser.error(e)
    
    qualities = reader.delivers_qualities
    
    # Wrap reader in batch iterator
    batch_size = options.batch_size or 1000
    reader = BatchIterator(reader, batch_size, options.max_reads)
    
    # Wrap iterator in progress bar
    if options.progress == "msg":
        reader = ProgressMessageReader(reader, batch_size)
    elif options.progress == "bar":
        reader = create_progress_reader(reader, options.max_reads, counter_magnitude)
    
    return (reader, qualities, quality_filename is not None)

class BatchIterator(object):
    def __init__(self, reader, size, max_reads=None):
        self.reader = reader
        self.iterable = enumerate(reader, 1)
        self.size = size
        self.max_reads = max_reads
        self.done = False
        self._empty_batch = [None] * size
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.done:
            raise StopIteration()
        
        try:
            read_index, record = next(self.iterable)
        except StopIteration:
            self.close()
            raise
        
        batch = copy.copy(self._empty_batch)
        batch[0] = record
        batch_index = 1
        max_size = self.size
        if self.max_reads:
            max_size = min(max_size, self.max_reads - read_index + 1)
        
        while batch_index < max_size:
            try:
                read_index, record = next(self.iterable)
                batch[batch_index] = record
                batch_index += 1
            except StopIteration:
                self.close()
                break
        
        if self.max_reads and read_index >= self.max_reads:
            self.close()
        
        if batch_index == self.size:
            return (batch_index, batch)
        else:
            return (batch_index, batch[0:batch_index])
    
    # py2x alias
    next = __next__
    
    def close(self):
        self.done = True
        self.reader.close()

class ProgressMessageReader(object):
    def __init__(self, iterable, batch_size, interval=1000000):
        self.iterable = iterable
        self.batch_size = batch_size
        self.interval = interval
        self.ctr = 0
    
    def __next__(self):
        value = self.iterable.next()
        if value:
            self.ctr += value[0]
            if self.ctr % self.interval < self.batch_size:
                now = time.time()
                logging.getLogger().info("Read {0} records in {1:.1f} seconds".format(self.ctr, now - self.start))
        return value
    
    next = __next__
    def __iter__(self):
        self.start = time.time()
        return self
    
    def close(self):
        logging.getLogger().info("Read a total of {} records".format(self.ctr))
        self.iterable.close()
        
def create_progress_reader(reader, max_reads=None, counter_magnitude="M"):
    try:
        return create_progressbar_reader(reader, max_reads, counter_magnitude)
    except:
        pass
    
    try:
        return create_tqdm_reader(reader)
    except:
        pass
    
    logging.getLogger().warn("No progress bar library available")
    return reader

def create_progressbar_reader(reader, max_reads=None, counter_magnitude="M"):
    import progressbar
    import progressbar.widgets
    import math

    class ProgressBarReader(progressbar.ProgressBar):
        def __init__(self, iterable, widgets, max_value=None):
            super(ProgressBarReader, self).__init__(
                widgets=widgets, max_value=max_value or progressbar.UnknownLength)
            self._iterable = iterable
            self.done = False
        
        def __next__(self):
            try:
                value = next(self._iterable)
                if self.start_time is None:
                    self.start()
                self.update(self.value + value[0])
                return value
            except StopIteration:
                self.close()
                raise
        
        def close(self):
            if not self.done:
                self.finish()
                self.done = True
            if not self._iterable.done:
                self._iterable.close()
    
    class MagCounter(progressbar.widgets.WidgetBase):
        def __init__(self, magnitude):
            suffix = ""
            if magnitude is None:
                div = 1.0
            else:
                div = float(MAGNITUDE[magnitude][0])
                suffix = magnitude
        
            self._format = lambda val: "{:.1f} {}".format(val / div, suffix)
    
        def __call__(self, progress, data):
            return self._format(data["value"])
        
    if max_reads:
        reader = ProgressBarReader(reader, [
            MagCounter(counter_magnitude), " Reads (", progressbar.Percentage(), ") ",
            progressbar.Timer(), " ", progressbar.Bar(), progressbar.AdaptiveETA()
        ], max_reads)
    else:
        reader = ProgressBarReader(reader, [
            MagCounter(counter_magnitude), " Reads", progressbar.Timer(),
            progressbar.AnimatedMarker()
        ])
    
    return reader

def create_tqdm_reader(reader):
    import tqdm
    return tqdm.tqdm(reader)

def create_modifiers(options, paired, qualities, has_qual_file, parser):
    adapter_parser = AdapterParser(
        colorspace=options.colorspace,
        max_error_rate=options.error_rate,
        min_overlap=options.overlap,
        read_wildcards=options.match_read_wildcards,
        adapter_wildcards=options.match_adapter_wildcards,
        indels=options.indels)

    try:
        adapters1 = adapter_parser.parse_multi(options.adapters, options.anywhere, options.front)
        adapters2 = adapter_parser.parse_multi(options.adapters2, options.anywhere2, options.front2)
    except IOError as e:
        if e.errno == errno.ENOENT:
            parser.error(e)
        raise
    except ValueError as e:
        parser.error(e)
    
    if not adapters1 and not adapters2 and not options.quality_cutoff and \
            options.nextseq_trim is None and \
            options.cut == [] and options.cut2 == [] and \
            options.cut_min == [] and options.cut_min2 == [] and \
            (options.minimum_length is None or options.minimum_length <= 0) and \
            options.maximum_length == sys.maxsize and \
            not has_qual_file and \
            options.max_n == -1 and not options.trim_n:
        parser.error("You need to provide at least one adapter sequence.")
    
    if options.aligner == 'insert' and (
            not adapters1 or len(adapters1) > 1 or adapters1[0].where != BACK or
            not adapters2 or len(adapters2) > 1 or adapters2[0].where != BACK):
        parser.error("Insert aligner requires a single 3' adapter for each read")
    
    if options.debug:
        for adapter in adapters1 + adapters2:
            adapter.enable_debug()
    
    merger = None
    if options.merge_overlapping:
        merger = MergeOverlapping(
            min_overlap=options.merge_min_overlap,
            error_rate=options.error_rate)
    
    modifiers = Modifiers(paired, merger)
    
    for op in options.op_order:
        if op == 'A' and (adapters1 or adapters2):
            # TODO: generalize this using some kind of factory class
            if options.aligner == 'insert':
                modifiers.add_modifier(InsertAdapterCutter,
                    adapter1=adapters1[0], adapter2=adapters2[0], action=options.action)
            else:
                modifiers.add_modifier_pair(AdapterCutter,
                    dict(adapters=adapters1, times=options.times, action=options.action),
                    dict(adapters=adapters2, times=options.times, action=options.action)
                )
        elif op == 'C' and (options.cut or options.cut2):
            modifiers.add_modifier_pair(UnconditionalCutter,
                dict(lengths=options.cut),
                dict(lengths=options.cut2)
            )
        elif op == 'G' and (options.nextseq_trim is not None):
            modifiers.add_modifier(NextseqQualityTrimmer,
                read=1, cutoff=options.nextseq_trim, base=options.quality_base)
        elif op == 'Q' and options.quality_cutoff:
            modifiers.add_modifier(QualityTrimmer,
                cutoff_front=options.quality_cutoff[0],
                cutoff_back=options.quality_cutoff[1],
                base=options.quality_base)
    
    if options.bisulfite:
        if isinstance(options.bisulfite, str):
            if options.bisulfite == "non_directional":
                modifiers.add_modifier(NonDirectionalBisulfiteTrimmer)
            elif options.bisulfite == "epignome":
                modifiers.add_modifier(EpiGnomeBisulfiteTrimmer)
            elif options.bisulfite == "truseq":
                # TODO: how to trim?
                pass
            elif options.bisulfite == "swift":
                modifiers.add_modifier(SwiftBisulfiteTrimmer)
        else:
            if options.bisulfite[0]:
                modifiers.add_modifier(MinCutter, read=1, **(options.bisulfite[0]))
            if len(options.bisulfite) > 1 and options.bisulfite[1]:
                modifiers.add_modifier(MinCutter, read=2, **(options.bisulfite[1]))
    
    if options.trim_n:
        modifiers.add_modifier(NEndTrimmer)
    
    if options.cut_min or options.cut_min2:
        modifiers.add_modifier_pair(MinCutter,
            dict(lengths=options.cut_min),
            dict(lengths=options.cut_min2)
        )
    
    if options.length_tag:
        modifiers.add_modifier(LengthTagModifier, length_tag=options.length_tag)
    
    if options.strip_suffix:
        modifiers.add_modifier(SuffixRemover, suffixes=options.strip_suffix)
    
    if options.prefix or options.suffix:
        modifiers.add_modifier(PrefixSuffixAdder, prefix=options.prefix, suffix=options.suffix)
    
    if options.double_encode:
        modifiers.add_modifier(DoubleEncoder)
    
    if options.zero_cap and qualities:
        modifiers.add_modifier(ZeroCapper, quality_base=options.quality_base)
    
    if options.trim_primer:
        modifiers.add_modifier(PrimerTrimmer)
    
    return (modifiers, len(adapters1) + len(adapters2))

def create_filters(options, paired, min_affected):
    filters = Filters(FilterFactory(paired, min_affected))
    
    if options.minimum_length is not None and options.minimum_length > 0:
        filters.add_filter(TooShortReadFilter, options.minimum_length)
    
    if options.maximum_length < sys.maxsize:
        filters.add_filter(TooLongReadFilter, options.maximum_length)
    
    if options.max_n >= 0:
        filters.add_filter(NContentFilter, options.max_n)
    
    if options.discard_trimmed:
        filters.add_filter(TrimmedFilter)
    
    if options.discard_untrimmed or options.untrimmed_output:
        filters.add_filter(UntrimmedFilter)
    
    return filters

def create_formatters(options, qualities, default_outfile):
    seq_formatter_args = dict(
        qualities=qualities,
        colorspace=options.colorspace,
        interleaved=options.interleaved
    )
    formatters = Formatters(options.output, seq_formatter_args)
    force_create = []
        
    if (options.merge_overlapping and options.merged_output):
        formatters.add_seq_formatter(MergedReadFilter, options.merged_output)
        
    if (options.minimum_length is not None
            and options.minimum_length > 0
            and options.too_short_output):
        formatters.add_seq_formatter(TooShortReadFilter,
            options.too_short_output, options.too_short_paired_output)

    if options.maximum_length < sys.maxsize and options.too_long_output is not None:
        formatters.add_seq_formatter(TooLongReadFilter,
            options.too_long_output, options.too_long_paired_output)
    
    if not formatters.multiplexed:
        if options.output is not None:
            formatters.add_seq_formatter(NoFilter,
                options.output, options.paired_output)
            if options.output != "-":
                force_create.append(options.output)
                if options.paired_output is not None:
                    force_create.append(options.paired_output)
        elif not (options.discard_trimmed and options.untrimmed_output):
            formatters.add_seq_formatter(NoFilter, default_outfile)
            if default_outfile != "-":
                force_create.append(default_outfile)
    
    if not options.discard_untrimmed:
        if formatters.multiplexed:
            untrimmed = options.untrimmed_output or options.output.format(name='unknown')
            formatters.add_seq_formatter(UntrimmedFilter, untrimmed)
            formatters.add_seq_formatter(NoFilter, untrimmed)
        elif options.untrimmed_output:
            formatters.add_seq_formatter(UntrimmedFilter,
                options.untrimmed_output, options.untrimmed_paired_output)

    if options.rest_file:
        formatters.add_info_formatter(RestFormatter(options.rest_file))
    if options.info_file:
        formatters.add_info_formatter(InfoFormatter(options.info_file))
    if options.wildcard_file:
        formatters.add_info_formatter(WildcardFormatter(options.wildcard_file))
    
    return (formatters, force_create)

def run_serial(reader, modifiers, filters, formatters, writers):
    try:
        n = 0
        total_bp1 = 0
        total_bp2 = 0
        
        for batch_size, batch in reader:
            n += batch_size
            result = defaultdict(lambda: [])
            for record in batch:
                reads, bp = modifiers.modify(record)
                total_bp1 += bp[0]
                total_bp2 += bp[1]
                dest = filters.filter(*reads)
                formatters.format(result, dest, *reads)
            result = dict((path, "".join(strings))
                for path, strings in result.items())
            writers.write_result(result)
        
        return (0, Summary(
            collect_process_statistics(n, total_bp1, total_bp2, modifiers, filters, formatters),
            summarize_adapters(modifiers.get_modifiers(AdapterCutter)),
            modifiers.get_trimmer_classes()
        ).finish())
    
    except KeyboardInterrupt as e:
        logging.getLogger().error("Interrupted")
        return (130, None)
    except IOError as e:
        if e.errno == errno.EPIPE:
            return (1, None)
        raise
    except (FormatError, EOFError) as e:
        logging.getLogger().error("Atropos error", exc_info=True)
        return(1, None)
    finally:
        reader.close()
        writers.close()

if __name__ == '__main__':
    main()
