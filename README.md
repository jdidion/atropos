[![Travis CI](https://travis-ci.org/jdidion/atropos.svg?branch=1.1)](https://travis-ci.org/jdidion/atropos)
[![PyPi](https://img.shields.io/pypi/v/atropos.svg)](https://pypi.python.org/pypi/atropos)
[![DOI](https://zenodo.org/badge/61393086.svg)](https://zenodo.org/badge/latestdoi/61393086)

# Atropos

Atropos is tool for specific, sensitive, and speedy trimming of NGS reads. It is a fork of the venerable Cutadapt read trimmer (https://github.com/marcelm/cutadapt, [DOI:10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)), with the primary improvements being:

1. Multi-threading support, including an extremely fast "parallel write" mode.
2. Implementation of a new insert alignment-based trimming algorithm for paired-end reads that is substantially more sensitive and specific than the original Cutadapt adapter alignment-based algorithm. This algorithm can also correct mismatches between the overlapping portions of the reads.
3. Options for trimming specific types of data (miRNA, bisulfite-seq).
4. A new command ('detect') that will detect adapter sequences and other potential contaminants.
5. A new command ('error') that will estimate the sequencing error rate, which helps to select the appropriate adapter- and quality- trimming parameter values.
6. A new command ('qc') that generates read statistics similar to FastQC. The trim command can also compute read statistics both before and after trimming (using the '--stats' option).
7. Improved summary reports, including support for serialization formats (JSON, YAML, pickle), support for user-defined templates (via the optional Jinja2 dependency), and integration with [MultiQC](http://multiqc.info).
8. The ability to merge overlapping reads (this is experimental and the functionality is limited).
9. The ability to write the summary report and log messages to separate files.
10. The ability to read SAM/BAM files and read/write interleaved FASTQ files.
11. Direct trimming of reads from an SRA accession.
12. A progress bar, and other minor usability enhancements.

## Manual installation

Atropos is available from [pypi](https://pypi.python.org/pypi/atropos) and can be installed using `pip`.

First install dependencies:

* Required
    * Python 3.4.5+
        * Python 2.x is NOT supported
        * Python 3.5.2 has a [bug](https://bugs.python.org/issue28387) that [has been reported to cause random SEGFAULTs](https://github.com/jdidion/atropos/issues/88)
    * Cython 0.25.2+ (`pip install Cython`)
* Maybe python libraries
    * pytest (for running unit tests)
    * progressbar2 or tqdm (progressbar support)
    * pysam (SAM/BAM input)
    * khmer 2.0+ (for detecting low-frequency adapter contamination)
    * jinja2 (for user-defined report formats)
    * ngstream (for SRA streaming), which requires [ngs](https://github.com/ncbi/ngs)

Pip can be used to install atropos and optional dependencies, e.g.:

pip install atropos[tqdm,pysam,ngstream]

## Conda

There is an Atropos recipe in [Bioconda](https://anaconda.org/bioconda/atropos).

`conda install -c bioconda atropos`

## Docker

A [Docker image](https://hub.docker.com/r/jdidion/atropos/) is available for Atropos in Docker Hub.

`docker run jdidion/atropos <arguments>`

## Usage

Atropos is almost fully backward-compatible with cutadapt. If you currently use cutadapt, you can simply install Atropos and then substitute the executable name in your command line, with one key difference: you need to use options to specify input file names. For example:

```
atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz -se reads.fq.gz
```

To take advantage of multi-threading, set the `--threads` option:

```
atropos --threads 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz -se reads.fq.gz
```

To take advantage of the new aligner (if you have paired-end reads with 3' adapters), set the `--aligner` option to 'insert':

```
atropos --aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed.1.fq.gz -p trimmed.2.fq.gz \
  -pe1 reads.1.fq.gz -pe2 reads.2.fq.gz
```

See the [Documentation](https://atropos.readthedocs.org/) for more complete usage information.

## Publications

Atropos is [published](https://peerj.com/articles/3720/) in PeerJ.

Please cite as:

> Didion JP, Martin M, Collins FS. (2017) Atropos: specific, sensitive, and speedy trimming of sequencing reads. PeerJ 5:e3720 https://doi.org/10.7717/peerj.3720

The results in the paper can be fully reproduced using the workflow defined in the [paper](paper/README.md) directory.

The citation for the original Cutadapt paper is:

> Marcel Martin. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200


## Links

* [Documentation](https://atropos.readthedocs.org/)
* [Source code](https://github.com/jdidion/atropos/)
* [Report an issue](https://github.com/jdidion/atropos/issues)
* [Code of conduct](https://github.com/jdidion/atropos/CODE_OF_CONDUCT.md)
* [Contributing](https://github.com/jdidion/atropos/CONTRIUBTING.md)

## Roadmap

### 1.2

* Migrate to [xphyle](https://github.com/jdidion/xphyle) for file management.
* Migrate to [pokrok](https://github.com/jdidion/pokrok) for progress bar management.
* Accept multiple input files.
* Support SAM output (including #33).
* Direct streaming and trimming of reads from SRA and htsget using [ngstream](https://github.com/jdidion/ngstream).
* Read "cropping" (#50)
* Support for ThruPlex-style adapters (in which barcode is part of query sequence; #55)
* Accessibility:
    * Create recipe for homebrew.
    * Automatically update conda and homebrew recipes for each release.
    * Create Galaxy tool description using [argparse2tool](https://github.com/erasche/argparse2tool#cwl-specific-functionality).
* Improve documentation (#24)
* Port over improvements in latest versions of Cutadapt https://cutadapt.readthedocs.io/en/stable/
* Switch to using entry point instead of Atropos executable.

### 1.3

* Add auto-trimming mode for paired-end reads.
* Support for UMIs.
* Provide PacBio- and nanopore-specific options (https://github.com/marcelm/cutadapt/issues/120).
* Provide option for RNA-seq data that will trim polyA sequence.
* Add formal config file support (#53)
* Automate crash reporting using [sentry](https://github.com/getsentry/raven-python).
* Look at [NGMerge] for improving read merging: https://github.com/harvardinformatics/NGmerge
* Look at replacing pysam with [pybam](https://github.com/JohnLonginotto/pybam)

### 1.4

* Currently, InsertAligner requires a single 3' adapter for each end. Adapter trimming will later be generalized so that A) the InsertAligner can handle multiple matched pairs of adapters and/or B) multiple different aligners can be used for different adapters.
* Integrate with [AdapterBase](https://github.com/NCBI-Hackathons/OnlineAdapterDatabase) for improved matching of detected contaminants to known adapters, automated trimming of datasets with known adapters, and (opt-in) submission of adapter information for novel datasets.
* Migrate to seqio (https://github.com/jdidion/seqio) for reading/writing sequence files.
    * Also look at using HTSeq instead: https://github.com/simon-anders/htseq
* General-purpose read filtering based on read ID: https://github.com/marcelm/cutadapt/issues/107.
* Currently, SAM/BAM input files must be name sorted; add an option to 1) pre-sort reads inline using samtools or sambamba, or 2) cache each read in memory until its mate is found.

### 1.5

* Provide more user control over anchoring of adapters: https://github.com/marcelm/cutadapt/issues/53.
* Enable user to define custom read structure: https://github.com/nh13/read-structure-examples
* Support for paired-end demultiplexing
* Demultiplexing based on barcodes: https://github.com/marcelm/cutadapt/issues/118.
* Consider supporting different error rates for read1 vs read2.
* Add a ClipOverlapping modifier that will remove read overlaps (as opposed to merging).
* Look more closely at providing solutions to the Illumina two-color chemistry issue:
    * Provide and option to exempt G calls from the assessment of quality
    * Trim 3′ Gs from reads
* Also look at addressing any issues with one-color chemistry (iSeq).
* Consider whether to support trimming/QC of raw IonTorrent data.

### 1.6

* Switch to using Click for CLI.
* Implement a public plugin API.
* Add more logging and convert log messages from old-style to new-style format strings.
* Add option to estimate bisulfite conversion rate from filled-in cytosine methylation status in reads that were MspI-digested.
* CPU and memory profiling. Try out:
    * https://github.com/nvdv/vprof
    * https://github.com/what-studio/profiling
    * https://github.com/fabianp/memory_profiler
    * https://github.com/rkern/line_profiler#line-profiler
* Look at some new trimming/qc programs; decide whether to add to benchmarks and/or incorporate any of their features
    * https://github.com/OpenGene/fastp/issues
    * http://tagcleaner.sourceforge.net/
    * https://github.com/mdshw5/fastqp/blob/master/README.md

### 2.0

* Simplification of command line options, perhaps by further splitting functionality up into different sub-commands, but also by more intelligent choices for default option values based on context.

* Consider adding additional report formats
* Add fuzzy matching: https://github.com/Daniel-Liu-c0deb0t/Java-Fuzzy-Search
* https://github.com/marcelm/cutadapt/issues/112
* Performance enhancements using
    * http://numba.pydata.org/
    * https://github.com/undefx/vecpy
    * https://github.com/serge-sans-paille/pythran
    * https://github.com/IntelLabs/hpat
    * https://github.com/cupy/cupy
* >90% test coverage
* Fuzz testing using AFL
    * http://lcamtuf.coredump.cx/afl/
    * https://github.com/jwilk/python-afl

### Beyond 2.0

* Implement additional alternate alignment algorithms.
* Implement the error detection algorithm in ADEPT: https://github.com/LANL-Bioinformatics/ADEPT
* Explore new quality trimming algorithms
  * UrQt: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450468/
  * InfoTrim: github.com/JacobPorter/InfoTrim
* Scythe is an interesting new trimmer. Depending on how the benchmarks look in the forthcoming paper, we will add it to the list of tools we compare against Atropos, and perhaps implement their Bayesian approach for adapter match.
* Experiment with replacing the multicore implementation with an asyncio-based implementation (using ProcessPoolExecutor and uvloop).
* Automatic adaptive tuning of queue sizes to maximize the balance between memory usage and latency.
* FastProNGS has some nice visualizations that could be included, rather than relying on MultiQC: https://github.com/Megagenomics/FastProNGS

While we consider the command-line interface to be stable, the internal code organization of Atropos is likely to change. At this time, we recommend to not directly interface with Atropos as a library (or to be prepared for your code to break). The internal code organization will be stabilized as of version 2.0, which is planned for sometime in 2017.

If you would like to suggest additional enhancements, you can submit issues and/or pull requests at our GitHub page.
