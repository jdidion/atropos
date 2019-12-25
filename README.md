[![Travis CI](https://travis-ci.org/jdidion/atropos.svg?branch=develop)](https://travis-ci.org/jdidion/atropos)
[![PyPi](https://img.shields.io/pypi/v/atropos.svg)](https://pypi.python.org/pypi/atropos)
[![DOI](https://zenodo.org/badge/61393086.svg)](https://zenodo.org/badge/latestdoi/61393086)
[![Coverage Status](https://img.shields.io/coveralls/jdidion/atropos/master.svg)](https://coveralls.io/github/jdidion/atropos?branch=develop)

# Atropos

Atropos is tool for specific, sensitive, and speedy trimming of NGS reads. It is a fork of the venerable Cutadapt read trimmer (https://github.com/marcelm/cutadapt, [DOI:10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)), with the primary improvements being:

1. Multi-threading support, including an extremely fast "parallel write" mode.
2. Implementation of a new insert alignment-based trimming algorithm for paired-end reads that is substantially more sensitive and specific than the original Cutadapt adapter alignment-based algorithm. This algorithm can also correct mismatches between the overlapping portions of the reads.
3. Options for trimming specific types of data (miRNA, bisulfite-seq).
4. A new command ('detect') that will detect adapter sequences and other potential contaminants.
5. A new command ('error') that will estimate the sequencing error rate, which helps to select the appropriate adapter- and quality- trimming parameter values.
6. A new command ('qc') that`                                           ``` generates read statistics similar to FastQC. The trim command can also compute read metrics both before and after trimming (using the '--metrics' option).
7. Improved summary reports, including support for serialization formats (JSON, YAML, pickle), support for user-defined templates (via the optional Jinja2 dependency), and integration with [MultiQC](http://multiqc.info).
8. The ability to merge overlapping reads (this is experimental and the functionality is limited).
9. The ability to write the summary report and log messages to separate files.
10. The ability to read/write SAM, BAM, and interleaved FASTQ files.
11. Direct trimming of reads from an SRA accession.
12. A progress bar, and other minor usability enhancements.

## Manual installation

Atropos is available from [pypi](https://pypi.python.org/pypi/atropos) and can be installed using `pip`.

First install dependencies:

* Required
    * Python 3.6+ (python 2.x is NOT supported)
    * Cython 0.25.2+/0.29+/0.29.14+, depending on whether you're using python 3.6/3.7/3.8 (`pip install Cython`)
    * [loguru]()
    * [pokrok]() 0.2.0+
    * [xphyle]() 4.2.1+
* Optional
    * pytest (for running unit tests)
    * bamnostic or pysam (SAM/BAM support)
    * khmer 2.0+ (for detecting low-frequency adapter contamination)
    * jinja2 (for user-defined report formats)
    * ngstream (for SRA streaming), which requires [ngs](https://github.com/ncbi/ngs)

Pip can be used to install atropos and optional dependencies, e.g.:

`pip install atropos[tqdm,bamnostic,ngstream]`

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

## Using Atropos as a library

While we consider the command-line interface to be stable, the internal code organization of Atropos is likely to change. At this time, we recommend to not directly interface with Atropos as a library (or to be prepared for your code to break).

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
