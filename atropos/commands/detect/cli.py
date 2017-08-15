"""Command-line interface for the detect command.
"""
from atropos.io import STDOUT, STDERR
from atropos.commands.cli import (
    BaseCommandParser, positive, readable_url, writeable_file,
    readwriteable_file, probability)

class CommandParser(BaseCommandParser):
    name = 'detect'
    usage = """
atropos detect -se input.fastq
atropos detect -pe1 in1.fq -pe2 in2.fq
"""
    description = """
Detect adapter sequences directly from read sequences. Use this command if you
are unsure if your data has adapters, or if you know that it has adapters but
you don't know what are the adapter sequences.
"""
    
    def add_command_options(self):
        parser = self.parser
        parser.set_defaults(
            max_reads=10000,
            counter_magnitude="K")
        
        group = self.add_group("Adapter Detection")
        group.add_argument(
            "-d",
            "--detector",
            choices=('known', 'heuristic', 'khmer'), default=None,
            help="Which detector to use. (automatically choose based on other "
                 "options)")
        group.add_argument(
            "-k",
            "--kmer-size",
            type=positive(), default=12,
            help="Size of k-mer used to scan reads for adapter sequences. (12)")
        group.add_argument(
            "-e",
            "--past-end-bases", nargs="*", default=('A',),
            help="On Illumina, long runs of A (and sometimes other bases) "
                 "can signify that the sequencer has read past the end of a "
                 "fragment that is shorter than the read length + adapter "
                 "length. Those bases will be removed from any sequencers "
                 "before looking for matching contaminants. Can also be a "
                 "regular expression.")
        group.add_argument(
            "-i",
            "--include-contaminants",
            choices=('all','known','unknown'), default='all',
            help="What conaminants to search for: all, only known "
                 "adapters/contaminants ('known'), or only unknown "
                 "contaminants ('unknown'). (all)")
        group.add_argument(
            "-x",
            "--known-contaminant",
            action="append", dest='known_adapter', default=None,
            help="Pass known contaminants in on the commandline as "
                 "'name=sequence'. Can be specified multiple times.")
        group.add_argument(
            "-F",
            "--known-contaminants-file",
            type=readable_url, action="append", dest='known_adapters_file',
            default=None,
            help="Points to FASTA File or URL with known contaminants.")
        group.add_argument(
            "--no-default-contaminants",
            action="store_false", dest="default_adapters", default=True,
            help="Don't fetch the default contaminant list (which is currently "
                 "stored as a GitHub gist).")
        group.add_argument(
            "--contaminant-cache-file",
            type=readwriteable_file, dest='adapter_cache_file',
            default='.adapters',
            help="File where known contaminant sequences will be cached, "
                 "unless --no-cache-contaminants is set.")
        group.add_argument(
            "--no-cache-contaminants",
            action="store_false", dest="cache_adapters", default=True,
            help="Don't cache contaminant list as '.contaminants' in the "
                 "working directory.")
        
        group = self.add_group("Known Detector Options")
        group.add_argument(
            "--min-kmer-match-frac",
            type=probability, default=0.5,
            help="Minimum fraction of contaminant kmers that must be found "
                 "in a read sequence to be considered a match.")
        
        group = self.add_group("Heuristic Detector Options")
        group.add_argument(
            "--min-frequency",
            type=probability, default=0.001,
            help="Minimum frequency required to retain a k-mer.")
        group.add_argument(
            "--min-contaminant-match-frac",
            type=probability, default=0.9,
            help="Minimum fraction of nucleotides that must align for a "
                 "detected contaminant to match a known adapter sequence.")
        
        group = self.add_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writeable_file, default=STDOUT, metavar="FILE",
            help="File in which to write the summary of detected adapters. "
                 "(stdout)")
        group.add_argument(
            "-O",
            "--output-formats",
            nargs="*", choices=("txt", "fasta", "json", "yaml", "pickle"), 
            default=None, metavar="FORMAT", dest="report_formats",
            help="Report type(s) to generate. If multiple, '--output' "
                 "is treated as a prefix and the appropriate extensions are "
                 "appended. If unspecified, the format is guessed from the "
                 "file extension. Supported formats are: txt, json, yaml, "
                 "pickle, fasta. Additional arguments for fasta output are "
                 "provided via the '--fasta' option. See the documentation "
                 "for a full description of the structured output "
                 "(json/yaml/pickle formats).")
        group.add_argument(
            "--fasta",
            nargs="*", choices=("union", "perinput"), default=None, 
            metavar="OPTIONS",
            help="Additional arguments for fasta output. Adds 'fasta' to the "
                 "list of output formats if not already specified. Options "
                 "are: perinput=generate one output file per input file, "
                 "union=generate a single output file with all sequences "
                 "merged.")
        group.add_argument(
            "-m",
            "--max-adapters",
            type=positive(), default=None,
            help="The maximum number of candidate adapters to report. "
                 "(report all)")
    
    def validate_command_options(self, options):
        options.report_file = options.output
        is_std = options.report_file in (STDOUT, STDERR)
        if options.fasta:
            if is_std and 'perinput' in options.fasta:
                self.parser.error("Per-input fasta cannot be written to stdout")
            if not options.report_formats:
                options.report_formats = ['fasta']
            elif 'fasta' not in options.report_formats:
                options.report_formats = list(options.report_formats) + ['fasta']
        elif is_std and options.report_formats and 'fasta' in options.report_formats:
            options.fasta = ['union']
