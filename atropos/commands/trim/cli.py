"""Trim command line interface.
"""
import logging
import sys
from atropos.commands.cli import (
    BaseCommandParser, configure_threads, parse_stat_args, readable_file,
    readwriteable_file, writeable_file, positive, probability, CharList,
    Delimited, int_or_str)
from atropos.io import STDOUT, STDERR

class CommandParser(BaseCommandParser):
    name = 'trim'
    usage = """
atropos trim -a ADAPTER [options] [-o output.fastq] -se input.fastq
atropos trim -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq -pe1 in1.fastq -pe2 in2.fastq
"""
    description = """
Trim adapters and low-quality bases, and perform other NGS preprocessing. This
command provides most of Atropos' functionality.
"""
    details = """
Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. The reverse complement is *not* automatically
searched. All reads from input.fastq will be written to output.fastq with the
adapter sequence removed. Adapter matching is error-tolerant. Multiple adapter
sequences can be given (use further -a options), but only the best-matching
adapter will be removed.

Input may also be in FASTA, SAM, or BAM format. Compressed input and output is
supported and auto-detected from the file name (.gz, .xz, .bz2). Use the file
name '-' for standard input/output. Without the -o option, output is sent to
standard output.
"""

    def add_command_options(self):
        self.parser.set_defaults(
            zero_cap=None,
            action='trim',
            batch_size=None,
            known_adapter=None)

        group = self.add_group(
            "Adapters",
            title="Finding adapters",
            description=\
                "Parameters -a, -g, -b specify adapters to be removed from "
                "each read (or from the first read in a pair if data is "
                "paired). If specified multiple times, only the best matching "
                "adapter is trimmed (but see the --times option). When the "
                "special notation 'file:FILE' is used, adapter sequences are "
                "read from the given FASTA file. When the --adapter-file "
                "option is used, adapters can be specified by name rather "
                "than sequence.")
        group.add_argument(
            "-a",
            "--adapter",
            action="append", default=[], metavar="ADAPTER", dest="adapters",
            help="Sequence of an adapter ligated to the 3' end (paired data: "
                 "of the first read). The adapter and subsequent bases are "
                 "trimmed. If a '$' character is appended ('anchoring'), the "
                 "adapter is only found if it is a suffix of the read. (none)")
        group.add_argument(
            "-g",
            "--front",
            action="append", default=[], metavar="ADAPTER",
            help="Sequence of an adapter ligated to the 5' end (paired data: "
                 "of the first read). The adapter and any preceding bases are "
                 "trimmed. Partial matches at the 5' end are allowed. If a "
                 "'^' character is prepended ('anchoring'), the adapter is "
                 "only found if it is a prefix of the read. (none)")
        group.add_argument(
            "-b",
            "--anywhere",
            action="append", default=[], metavar="ADAPTER",
            help="Sequence of an adapter that may be ligated to the 5' or 3' "
                 "end (paired data: of the first read). Both types of matches "
                 "as described under -a und -g are allowed. If the first base "
                 "of the read is part of the match, the behavior is as with "
                 "-g, otherwise as with -a. This option is mostly for rescuing "
                 "failed library preparations - do not use if you know which "
                 "end your adapter was ligated to! (none)")
        group.add_argument(
            "-F",
            "--known-adapters-file",
            type=readable_file, action="append", default=None,
            help="Path or URL of a FASTA file containing adapter sequences.")
        group.add_argument(
            "--no-default-adapters",
            action="store_false", dest="default_adapters", default=True,
            help="Don't fetch the default adapter list (which is currently "
                 "stored in GitHub).")
        group.add_argument(
            "--adapter-cache-file",
            type=readwriteable_file, default='.adapters',
            help="File where adapter sequences will be cached, unless "
                 "--no-cache-adapters is set.")
        group.add_argument(
            "--no-cache-adapters",
            action="store_false", dest="cache_adapters", default=True,
            help="Don't cache adapters list as '.adapters' in the working "
                 "directory.")
        group.add_argument(
            "--no-trim",
            action='store_const', dest='action', const=None,
            help="Match and redirect reads to output/untrimmed-output as "
                 "usual, but do not remove adapters. (no)")
        group.add_argument(
            "--mask-adapter",
            action='store_const', dest='action', const='mask',
            help="Mask adapters with 'N' characters instead of trimming "
                 "them. (no)")
        group.add_argument(
            "--gc-content",
            type=probability, default=0.5,
            help="Expected GC content of sequences.")

        ## Arguments specific to the choice of aligner

        group.add_argument(
            "--aligner",
            choices=('adapter', 'insert'), default='adapter',
            help="Which alignment algorithm to use for identifying adapters. "
                 "Currently, you can choose between the semi-global alignment "
                 "strategy used in Cutdapt ('adapter') or the more accurate "
                 "insert-based alignment algorithm ('insert'). Note that "
                 "insert-based alignment can only be used with paired-end "
                 "reads containing 3' adapters. New algorithms are being "
                 "implemented and the default is likely to change. (adapter)")

        # TODO: all the different matching options are pretty confusing. Either
        # explain their usage better in the docs or find a way to simplify the
        # choices.

        # Arguments for adapter match
        group.add_argument(
            "-e",
            "--error-rate",
            type=probability, default=None,
            help="Maximum allowed error rate for adapter match (no. of errors "
                 "divided by the length of the matching region). (0.1)")
        group.add_argument(
            "--indel-cost",
            type=positive(int, True), default=None, metavar="COST",
            help="Integer cost of insertions and deletions during adapter "
                 "match. Substitutions always have a cost of 1. (1)")
        group.add_argument(
            "--no-indels",
            action='store_false', dest='indels', default=True,
            help="Allow only mismatches in alignments. (allow both mismatches "
                 "and indels)")
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
            help="If the overlap between the read and the adapter is shorter "
                 "than MINLENGTH, the read is not modified. Reduces the no. "
                 "of bases trimmed due to random adapter matches. (3)")
        group.add_argument(
            "--adapter-max-rmp",
            type=probability, default=None, metavar="PROB",
            help="If no minimum overlap (-O) is specified, then adapters are "
                 "only matched when the probabilty of observing k out of n "
                 "matching bases is <= PROB. (1E-6)")

        # Arguments for insert match
        group.add_argument(
            "--insert-max-rmp",
            type=probability, default=1E-6, metavar="PROB",
            help="Overlapping inserts only match when the probablity of "
                 "observing k of n matching bases is <= PROB. (1E-6)")
        group.add_argument(
            "--insert-match-error-rate",
            type=probability, default=None,
            help="Maximum allowed error rate for insert match (no. of errors "
                 "divided by the length of the matching region). (0.2)")
        group.add_argument(
            "--insert-match-adapter-error-rate",
            type=probability, default=None,
            help="Maximum allowed error rate for matching adapters after "
                 "successful insert match (no. of errors divided by the length "
                 "of the matching region). (0.2)")

        # Arguments for merging and error correction
        # TODO: add RMP parameter for MergeOverlap
        group.add_argument(
            "-R",
            "--merge-overlapping",
            action="store_true", default=False,
            help="Merge read pairs that overlap into a single sequence. This "
                 "is experimental. (no)")
        group.add_argument(
            "--merge-min-overlap",
            type=positive(float, True), default=0.9,
            help="The minimum overlap between reads required for merging. If "
                 "this number is (0,1.0], it specifies the minimum length as "
                 "the fraction of the length of the *shorter* read in the "
                 "pair; otherwise it specifies the minimum number of "
                 "overlapping base pairs (with an absolute minimum of 2 bp). "
                 "(0.9)")
        group.add_argument(
            "--merge-error-rate",
            type=probability, default=None,
            help="The maximum error rate for merging. (0.2)")
        group.add_argument(
            "--correct-mismatches",
            choices=("liberal", "conservative", "N"), default=None,
            help="How to handle mismatches while aligning/merging. 'Liberal' "
                 "and 'conservative' error correction both involve setting the "
                 "base to the one with the best quality. They differ only when "
                 "the qualities are equal -- liberal means set it to the base "
                 "from the read with the overall best median base quality, "
                 "while conservative means to leave it unchanged. 'N' means to "
                 "set the base to N. If exactly one base is ambiguous, the "
                 "non-ambiguous base is always used. (no error correction)")

        group = self.add_group(
            "Modifications", title="Additional read modifications")
        group.add_argument(
            "--op-order",
            type=CharList(choices=('A','C','G','Q','W')), default="CGQAW",
            help="The order in which trimming operations are be applied. This "
                 "is a string of 1-5 of the following characters: A = adapter "
                 "trimming; C = cutting (unconditional); G = NextSeq trimming; "
                 "Q = quality trimming; W = overwrite poor quality reads. The "
                 "default is 'WCGQA' to maintain compatibility with "
                 "Cutadapt; however, this is likely to change to 'GAWCQ' in "
                 "the near future.")
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
            type=Delimited(data_type=positive(int, True), min_len=1, max_len=2),
            default=None, metavar="[5'CUTOFF,]3'CUTOFF",
            help="Trim low-quality bases from 5' and/or 3' ends of each read "
                 "before adapter removal. Applied to both reads if data is "
                 "paired. If one value is given, only the 3' end is trimmed. "
                 "If two comma-separated cutoffs are given, the 5' end is "
                 "trimmed with the first cutoff, the 3' end with the second. "
                 "(no)")
        group.add_argument(
            "-i",
            "--cut-min",
            type=int, action='append', default=[], metavar="LENGTH",
            help="Similar to -u, except that cutting is done AFTER adapter "
                 "trimming, and only if a minimum of LENGTH bases was not "
                 "already removed. (no)")
        group.add_argument(
            "--nextseq-trim",
            type=positive(), default=None, metavar="3'CUTOFF",
            help="NextSeq-specific quality trimming (each read). Trims also "
                 "dark cycles appearing as high-quality G bases "
                 "(EXPERIMENTAL). (no)")
        group.add_argument(
            "--trim-n",
            action='store_true', default=False,
            help="Trim N's on ends of reads. (no)")
        group.add_argument(
            "-x",
            "--prefix",
            default='',
            help="Add this prefix to read names. Use {name} to insert the name "
                 "of the matching adapter. (no)")
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
            help="Search for TAG followed by a decimal number in the "
                 "description field of the read. Replace the decimal number "
                 "with the correct length of the trimmed read. For example, "
                 "use --length-tag 'length=' to correct fields like "
                 "'length=123'. (no)")

        group = self.add_group(
            "Filtering", title="Filtering of processed reads")
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
            help="Discard trimmed reads that are shorter than LENGTH. "
                 "Reads that are too short even before adapter removal are "
                 "also discarded. In colorspace, an initial primer is not "
                 "counted. (0)")
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
            help="Discard reads with too many N bases. If COUNT is an integer, "
                 "it is treated as the absolute number of N bases. If it is "
                 "between 0 and 1, it is treated as the proportion of N's "
                 "allowed in a read. (no)")

        group = self.add_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writeable_file, metavar="FILE",
            help="Write trimmed reads to FILE. FASTQ or FASTA format is chosen "
                 "depending on input. The summary report is sent to standard "
                 "output. Use '{name}' in FILE to demultiplex reads into "
                 "multiple files. (write to standard output)")
        group.add_argument(
            "--info-file",
            type=writeable_file, metavar="FILE",
            help="Write information about each read and its adapter matches "
                 "into FILE. See the documentation for the file format. (no)")
        group.add_argument(
            "-r",
            "--rest-file",
            type=writeable_file, metavar="FILE",
            help="When the adapter matches in the middle of a read, write the "
                 "rest (after the adapter) into FILE. (no)")
        group.add_argument(
            "--wildcard-file",
            type=writeable_file, metavar="FILE",
            help="When the adapter has N bases (wildcards), write adapter "
                 "bases matching wildcard positions to FILE. When there are "
                 "indels in the alignment, this will often not be accurate. "
                 "(no)")
        group.add_argument(
            "--too-short-output",
            type=writeable_file, metavar="FILE",
            help="Write reads that are too short (according to length "
                 "specified by -m) to FILE. (no - too short reads are "
                 "discarded)")
        group.add_argument(
            "--too-long-output",
            type=writeable_file, metavar="FILE",
            help="Write reads that are too long (according to length specified "
                 "by -M) to FILE. (no - too long reads are discarded)")
        group.add_argument(
            "--untrimmed-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write reads that do not contain the adapter to FILE. (no - "
                 "untrimmed reads are written to default output)")
        group.add_argument(
            "--merged-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write reads that have been merged to this file. (merged "
                 "reads are discarded)")
        group.add_argument(
            "--report-file",
            type=writeable_file, default="-", metavar="FILE",
            help="Write report to file rather than stdout/stderr. (no)")
        group.add_argument(
            "--report-formats",
            nargs="*", choices=("txt", "json", "yaml", "pickle"), default=None,
            metavar="FORMAT",
            help="Report type(s) to generate. If multiple, '--report-file' "
                 "is treated as a prefix and the appropriate extensions are "
                 "appended. If unspecified, the format is guessed from the "
                 "file extension. Supported formats are: txt (legacy text "
                 "format), json, yaml, pickle. See the documentation for a "
                 "full description of the structured output (json/yaml/pickle "
                 "formats).")
        group.add_argument(
            "--stats",
            nargs="*", default=None,
            help="Which read-level statistics to compute. Can be 'none' "
                 "(default), 'pre': only compute pre-trimming stats; 'post': "
                 "only compute post-trimming stats; or 'both'. The keyword can "
                 "be followed by ':' and then additional configuration "
                 "parameters. E.g. 'pre:tiles' means to also collect "
                 "tile-level statistics (Illumina data only), and "
                 "'pre:tiles=<regexp>' means to use the specified regular "
                 "expression to extract key portions of read names to "
                 "collect the tile statistics.")

        group = self.add_group("Colorspace options")
        group.add_argument(
            "-d",
            "--double-encode",
            action='store_true', default=False,
            help="Double-encode colors (map 0,1,2,3,4 to A,C,G,T,N). (no)")
        group.add_argument(
            "-t",
            "--trim-primer",
            action='store_true', default=False,
            help="Trim primer base and the first color (which is the "
                 "transition to the first nucleotide). (no)")
        group.add_argument(
            "--strip-f3",
            action='store_true', default=False,
            help="Strip the _F3 suffix of read names. (no)")
        group.add_argument(
            "--maq", "--bwa",
            action='store_true', default=False,
            help="MAQ- and BWA-compatible colorspace output. This enables -c, "
                 "-d, -t, --strip-f3 and -y '/1'. (no)")
        group.add_argument(
            "--no-zero-cap",
            dest='zero_cap', action='store_false',
            help="Do not change negative quality values to zero in colorspace "
                 "data. By default, they are since many tools have problems "
                 "with negative qualities. (no)")
        group.add_argument(
            "-z",
            "--zero-cap",
            action='store_true',
            help="Change negative quality values to zero. This is enabled "
                 "by default when -c/--colorspace is also enabled. Use the "
                 "above option to disable it. (no)")

        group = self.add_group(
            "Paired",
            title="Paired-end options",
            description=\
                "The -A/-G/-B/-U/-I options work like their -a/-b/-g/-u/-i "
                "counterparts, but are applied to the second read in each "
                "pair.")
        group.add_argument(
            "-A",
            "--adapter2",
            action='append', dest='adapters2', default=[], metavar='ADAPTER',
            help="3' adapter to be removed from second read in a pair. (no)")
        group.add_argument(
            "-G",
            "--front2",
            action='append', dest='front2', default=[], metavar='ADAPTER',
            help="5' adapter to be removed from second read in a pair. (no)")
        group.add_argument(
            "-B",
            "--anywhere2",
            action='append', dest='anywhere2',  default=[], metavar='ADAPTER',
            help="5'/3 adapter to be removed from second read in a pair. (no)")
        group.add_argument(
            "-U",
            "--cut2",
            type=int, action='append', dest='cut2', default=[],
            metavar="LENGTH",
            help="Remove LENGTH bases from second read in a pair (see --cut). "
                 "(no)")
        group.add_argument(
            "-I",
            "--cut-min2",
            type=int, action='append', default=[], metavar="LENGTH",
            help="Similar to -U, except that cutting is done AFTER adapter "
                 "trimming, and only if a minimum of LENGTH bases was not "
                 "already removed (see --cut-min). (no)")
        group.add_argument(
            "-w",
            "--overwrite-low-quality",
            type=Delimited(data_type=positive(int, True), min_len=3, max_len=3),
            default=None, metavar="LOWQ,HIGHQ,WINDOW",
            help="When one read has mean quality < LOWQ and the other read has "
                 "mean quality >= HIGHQ over the first WINDOW bases, overwrite "
                 "the worse read with the better read.")
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
        # Setting the default for pair_filter to None allows us to find out
        # whether the option was used at all.
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
                 "--untrimmed-output when trimming paired-end reads. (no - "
                 "output to same file as trimmed reads)")
        group.add_argument(
            "--too-short-paired-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write second read in a pair to this file if pair is too "
                 "short. Use together with --too-short-output. (no - too short "
                 "reads are discarded)")
        group.add_argument(
            "--too-long-paired-output",
            type=writeable_file, default=None, metavar="FILE",
            help="Write second read in a pair to this file if pair is too "
                 "long. Use together with --too-long-output. (no - too long "
                 "reads are discarded)")

        group = self.add_group("Method-specific options")
        group = group.add_mutually_exclusive_group()
        group.add_argument(
            "--bisulfite",
            default=False, metavar="METHOD",
            help="Set default option values for bisulfite-treated data. The "
                 "argument specifies the type of bisulfite library (rrbs, "
                 "non-directional, non-directional-rrbs, truseq, epignome, or "
                 "swift) or custom parameters for trimming: "
                 "'<read1>[;<read2>]' where trimming parameters for each read "
                 "are: '<5' cut>,<3' cut>,<include trimmed>,<only trimmed>' "
                 "where 'include trimmed' is 1 or 0 for whether or not the "
                 "bases already trimmed during/prior to adapter trimming "
                 "should be counted towards the total bases to be cut and "
                 "'only trimmed' is 1 or 0 for whether or not only trimmed "
                 "reads should be further cut. (no)")
        group.add_argument(
            "--mirna",
            action="store_true", default=False,
            help="Set default option values for miRNA data. (no)")

        group = self.add_group(
            "Parallel", title="Parallel (multi-core) options")
        group.add_argument(
            "-T",
            "--threads",
            type=positive(int, True), default=None, metavar="THREADS",
            help="Number of threads to use for read trimming. Set to 0 to use "
                 "max available threads. (Do not use multithreading)")
        group.add_argument(
            "--no-writer-process",
            action="store_false", dest="writer_process", default=True,
            help="Do not use a writer process; instead, each worker thread "
                 "writes its own output to a file with a '.N' suffix. (no)")
        group.add_argument(
            "--preserve-order",
            action="store_true", default=False,
            help="Preserve order of reads in input files (ignored if "
                 "--no-writer-process is set). (no)")
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
        group.add_argument(
            "--result-queue-size",
            type=int_or_str, default=None, metavar="SIZE",
            help="Size of queue for batches of results to be written. "
                 "(THREADS * 100)")
        group.add_argument(
            "--compression",
            choices=("worker", "writer"), default=None,
            help="Where data compression should be performed. Defaults to "
                 "'writer' if system-level compression can be used and "
                 "(1 < threads < 8), otherwise defaults to 'worker'.")

    def validate_command_options(self, options):
        parser = self.parser
        paired = options.paired

        if not paired:
            if not options.output:
                parser.error("An output file is required")
            if options.untrimmed_paired_output:
                parser.error(
                    "Option --untrimmed-paired-output can only be used when "
                    "trimming paired-end reads (with option -p).")
        else:
            if not options.interleaved_output:
                if not options.output:
                    parser.error(
                        "When you use -p or --paired-output, you must also "
                        "use the -o option.")
                if not options.paired_output:
                    parser.error(
                        "When paired-end trimming is enabled via -A/-G/-B/-U, "
                        "a second output file needs to be specified via -p "
                        "(--paired-output).")
                if (bool(options.untrimmed_output) !=
                        bool(options.untrimmed_paired_output)):
                    parser.error(
                        "When trimming paired-end reads, you must use either "
                        "none or both of the --untrimmed-output/"
                        "--untrimmed-paired-output options.")
                if (options.too_short_output and
                        not options.too_short_paired_output):
                    parser.error(
                        "When using --too-short-output with ""paired-end "
                        "reads, you also need to use --too-short-paired-output")
                if (options.too_long_output and
                        not options.too_long_paired_output):
                    parser.error(
                        "When using --too-long-output with paired-end "
                        "reads, you also need to use --too-long-paired-output")

            # Any of these options switch off legacy mode
            if (options.adapters2 or options.front2 or options.anywhere2 or
                    options.cut2 or options.cut_min2 or
                    options.quality_cutoff or options.trim_n or
                    options.interleaved_input or options.pair_filter or
                    options.too_short_paired_output or
                    options.too_long_paired_output or
                    options.overwrite_low_quality):
                # Full paired-end trimming when both -p and -A/-G/-B/-U given
                # Read modifications (such as quality trimming) are applied also
                # to second read.
                paired = 'both'
            else:
                # Modify first read only, keep second in sync (-p given, but not
                # -A/-G/-B/-U). This exists for backwards compatibility
                # ('legacy mode').
                paired = 'first'

            options.paired = paired

        # Send report to stderr if main output will be going to stdout
        if options.output is None and options.report_file == STDOUT:
            options.report_file = STDERR

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
                # TODO: should also be checking that there is exactly one 3'
                # adapter for each read
                # TODO: have the aligner tell us whether it can be used based on
                # options?
            if options.indels and options.indel_cost is None:
                options.indel_cost = 3
            if options.overlap is None:
                options.overlap = 1
                if options.adapter_max_rmp is None:
                    options.adapter_max_rmp = 1E-6
            if options.insert_match_error_rate is None:
                options.insert_match_error_rate = options.error_rate or 0.2
            if options.insert_match_adapter_error_rate is None:
                options.insert_match_adapter_error_rate = \
                    options.insert_match_error_rate

        if options.merge_overlapping:
            if options.merged_output is None:
                logging.getLogger().warning(
                    "--merge-output is not set; merged reads will be discarded")
            if options.merge_error_rate is None:
                options.merge_error_rate = options.error_rate or 0.2

        if options.mirna:
            if not (options.adapters or options.front or options.anywhere):
                options.adapters = ['TGGAATTCTCGG'] # illumina small RNA adapter
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
                parser.error(
                    "Swift trimming is only compatible with paired-end reads")
            if options.bisulfite not in (
                    "rrbs", "non-directional", "truseq", "epignome", "swift",
                    "non-directional-rrbs"):
                def parse_bisulfite_params(arg):
                    """Parse bisulfite trimming command line arguments into
                    arguments to the Modifiers.
                    """
                    try:
                        parts = [int(part) for part in arg.split(",")]
                        assert len(parts) == 4
                        if parts[0] <= 0 and parts[1] <= 0:
                            return None
                        return dict(zip(
                            ("lengths", "count_trimmed", "only_trimmed"),
                            ((parts[0], -1 * parts[1]), (False,True)[parts[2]],
                             (False,True)[parts[3]])))
                    except:
                        parser.error(
                            "Invalidate format for bisulfite parameters")

                temp = [
                    parse_bisulfite_params(arg)
                    for arg in options.bisulfite.split(";")]
                if paired == "both" and len(temp) == 1:
                    temp = [temp[0], temp[0]]
                elif paired != "both" and len(temp) > 1:
                    parser.error(
                        "Too many bisulfite parameters for single-end reads")
                options.bisulfite = temp

        if options.overwrite_low_quality:
            if not paired:
                parser.error(
                    "--overwrite-low-quality is not valid for single-end reads")
            if (options.overwrite_low_quality[0] >
                    options.overwrite_low_quality[1]):
                parser.error(
                    "For --overwrite-low-quality, LOWQ must be <= HIGHQ")

        if options.quality_cutoff:
            if all(c <= 0 for c in options.quality_cutoff):
                options.quality_cutoff = None
            elif len(options.quality_cutoff) == 1:
                options.quality_cutoff = [0] + options.quality_cutoff

        if options.pair_filter is None:
            options.pair_filter = 'any'

        if ((options.discard_trimmed or options.discard_untrimmed) and
                (options.untrimmed_output is not None)):
            parser.error(
                "Only one of the --discard-trimmed, --discard-untrimmed "
                "and --untrimmed-output options can be used at the same time.")

        if options.output is not None and '{name}' in options.output:
            if options.discard_trimmed:
                parser.error(
                    "Do not use --discard-trimmed when demultiplexing.")
            if paired:
                parser.error(
                    "Demultiplexing not supported for paired-end files, yet.")

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
                parser.error(
                    "Using --anywhere with colorspace reads is currently not "
                    "supported (if you think this may be useful, contact the "
                    "author).")
            if options.match_read_wildcards:
                parser.error('IUPAC wildcards not supported in colorspace')
            options.match_adapter_wildcards = False
        else:
            if options.trim_primer:
                parser.error(
                    "Trimming the primer makes only sense in colorspace.")
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
            if (len(options.cut_min) == 2 and
                    options.cut_min[0] * options.cut_min[1] > 0):
                parser.error("You cannot remove bases from the same end twice.")

        if paired == 'both' and options.cut2:
            if len(options.cut2) > 2:
                parser.error("You cannot remove bases from more than two ends.")
            if len(options.cut2) == 2 and options.cut2[0] * options.cut2[1] > 0:
                parser.error("You cannot remove bases from the same end twice.")

        if paired == 'both' and options.cut_min2:
            if len(options.cut_min2) > 2:
                parser.error("You cannot remove bases from more than two ends.")
            if (len(options.cut_min2) == 2 and
                    options.cut_min2[0] * options.cut_min2[1] > 0):
                parser.error("You cannot remove bases from the same end twice.")

        if not options.stats or options.stats == 'none':
            options.stats = None
        else:
            stats = {}
            for stat_spec in options.stats:
                parts = stat_spec.split(':')
                name = parts[0]
                args = {} if len(parts) == 1 else parse_stat_args(parts[1])
                if name == 'both':
                    stats['pre'] = stats['post'] = args
                else:
                    stats[name] = args
            options.stats = stats

        if options.threads is not None:
            threads = configure_threads(options, parser)

            if options.compression is None:
                # Our tests show that with 8 or more threads, worker compression
                # is more efficient.
                if options.writer_process and 2 < threads < 8:
                    from atropos.io import compression
                    if compression.can_use_system_compression():
                        options.compression = "writer"
                    else:
                        options.compression = "worker"
                else:
                    options.compression = "worker"
            elif options.compression == "writer":
                if not options.writer_process:
                    parser.error(
                        "Writer compression and --no-writer-process are "
                        "mutually exclusive")
                elif threads == 2:
                    logging.getLogger().warning(
                        "Writer compression requires > 2 threads; using "
                        "worker compression instead")
                    options.compression = "worker"

            # Set queue sizes if necessary.
            # If we are using writer compression, the back-up will be in the
            # result queue, otherwise it will be in the read queue.
            if options.read_queue_size is None:
                options.read_queue_size = (
                    threads * (100 if options.compression == "writer" else 500))
            elif (
                    options.read_queue_size > 0 and
                    options.read_queue_size < threads):
                parser.error("Read queue size must be >= 'threads'")

            if options.result_queue_size is None:
                options.result_queue_size = (
                    threads * (100 if options.compression == "worker" else 500))
            elif (
                    options.result_queue_size > 0 and
                    options.result_queue_size < threads):
                parser.error("Result queue size must be >= 'threads'")

            max_queue_size = options.read_queue_size + options.result_queue_size
            if options.batch_size is None:
                options.batch_size = max(1000, max_queue_size / 10e6)
            elif options.batch_size * max_queue_size > 10e6:
                logging.getLogger().warning(
                    "Combination of batch size %d and total queue size %d "
                    "may lead to excessive memory usage",
                    options.batch_size, max_queue_size)

        if options.batch_size is None:
            options.batch_size = 1000
