"""Command-line interface for the detect command.
"""
from atropos.io import STDOUT
from atropos.commands.cli import (
    BaseCommandParser, positive, readable_url, writeable_file,
    readwriteable_file)

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
            counter_magnitude="K",
            report_formats=['txt'])
        
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
        
        group = self.add_group("Output")
        group.add_argument(
            "-o",
            "--output",
            type=writeable_file, default=STDOUT, metavar="FILE",
            help="File in which to write the summary of detected adapters. "
                 "(stdout)")
        group.add_argument(
            "-m",
            "--max-adapters",
            type=positive(), default=None,
            help="The maximum number of candidate adapters to report. "
                 "(report all)")
    
    def validate_command_options(self, options):
        options.report_file = options.output
