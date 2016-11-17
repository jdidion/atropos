![https://travis-ci.org/jdidion/atropos](https://travis-ci.org/jdidion/atropos.svg?branch=master)
![https://pypi.python.org/pypi/atropos](https://img.shields.io/pypi/v/atropos.svg?branch=master)

# Atropos

Atropos is tool for specific, sensitive, and speedy trimming of NGS reads. It is a fork of the venerable Cutadapt read trimmer (https://github.com/marcelm/cutadapt, [DOI:10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)), with the primary improvements being:

1. Multi-threading support, including an extremely fast "parallel write" mode.
2. Implementation of a new insert alignment-based trimming algorithm for paired-end reads that is substantially more sensitive and specific than the original Cutadapt adapter alignment-based algorithm. This algorithm can also correct mismatches between the overlapping portions of the reads.
3. Options for trimming specific types of data (miRNA, bisulfite-seq).
4. A new command ('detect') that will detect adapter sequences and other potential contaminants (this is experimental).
5. A new command ('error') that will estimate the sequencing error rate, which helps to select the appropriate adapter- and quality- trimming parameter values.
5. The ability to merge overlapping reads (this is experimental and the functionality is limited).
6. The ability to write the summary report and log messages to separate files.
7. The ability to read and write interleaved FASTQ files.
8. A progress bar, and other minor usability enhancements.

## Dependencies

* Python 3.3+ (python 2.x is NOT supported)
* Cython 0.24+ (`pip install Cython`)
* progressbar2 or tqdm (optional, if you want progressbar support)
* pysam (optional, if you want support for SAM/BAM input)
* khmer 2.0+ (`pip install khmer`) (optional, for detecting low-frequency adapter contamination)

## Installation

`pip install atropos`

## Usage

Atropos is almost fully backward-compatible with cutadapt. If you currently use cutadapt, you can simply install Atropos and then substitute the executable name in your command line, with one key difference: you need to use options to specify input file names. For example:

```
atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz -se reads.fq.gz
```

To take advantage of multi-threading, set the `--threads` option:

```
atropos --threads 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz -se reads.fq.gz
```

To take advantage of the new aligner (if you have paired-end reads with 3' adatpers), set the `--aligner` option to 'insert':

```
atropos --aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed.1.fq.gz -p trimmed.2.fq.gz \
  -pe1 reads.1.fq.gz -pe2 reads.2.fq.gz
```

See the [Documentation](https://atropos.readthedocs.org/) for more complete usage information.

## Links

* [Documentation](https://atropos.readthedocs.org/)
* [Source code](https://github.com/jdidion/atropos/)
* [Report an issue](https://github.com/jdidion/atropos/issues)

## Roadmap

### 1.1

* Migrate to xphyle (https://github.com/jdidion/xphyle) for file management.
* Integrate userstats (opt-in, of course) to gather user statistics and crash reports.
* Provide option for RNA-seq data that will trim polyA sequence.
* Accept multiple input files.
* Expand the list of contaminants that are detected by default.

### 1.2

* Provide PacBio-specific options (https://github.com/marcelm/cutadapt/issues/120).
* Currently, InsertAligner requires a single 3' adapter for each end. Adapter trimming will later be generalized so that A) the InsertAligner can handle multiple matched pairs of adapters and/or B) multiple different aligners can be used for different adapters.

### 1.3

* Migrate to seqio (https://github.com/jdidion/seqio) for reading/writing sequence files.
* General-purpose read filtering based on read ID: https://github.com/marcelm/cutadapt/issues/107.
* Currently, SAM/BAM input files must be name sorted; add an option to 1) pre-sort reads inline using samtools or sambamba, or 2) cache each read in memory until its mate is found.

### 1.4

* Provide more user control over anchoring of adapters: https://github.com/marcelm/cutadapt/issues/53.
* Support for paired-end demultiplexing (i.e. when barcodes are used in both paired-end adapters): https://github.com/marcelm/cutadapt/issues/118.
* Add option to estimate bisulfite conversion rate from filled-in cytosine methylation status in reads that were MspI-digested.
* Consider supporting different error rates for read1 vs read2.
* Add a ClipOverlapping modifier that will remove read overlaps (as opposed to merging).
* Add option to InsertAdapter to trim overhangs without adapter matching.

### 1.5

* Implement a public plugin API.
* Improvements to the summary report, and addition of a computer-parsable report for use in QC pipelines
    * https://github.com/marcelm/cutadapt/issues/112
    * Also look at the QCML used in ngs-bits

### 2.0

* Simplification of command line options, perhaps by splitting functionality up into different sub-commands, but also by more intelligent choices for default option values based on context.

### Beyond 2.0

* Implement additional alternate alignment algorithms.
* Implement the error detection algorithm in ADEPT: https://github.com/LANL-Bioinformatics/ADEPT
* Implement the quality trimming algorithm used in UrQt: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450468/
* Scythe is an interesting new trimmer. Depending on how the benchmarks look in the forthcomming paper, we will add it to the list of tools we compare against Atropos, and perhaps implement their Bayesian approach for adapter match.

While we consider the command-line interface to be stable, the internal code organization of Atropos is likely to change substantially. At this time, we recommend to not directly interface with Atropos as a library (or to be prepared for your code to break). The internal code organization will be stablized as of version 2.0, which is planned for early 2017.

If you would like to suggest additional enhancements, you can submit issues and/or pull requests at our GitHub page.

## Citations

The citation for the original Cutadapt paper is:
 
> Marcel Martin. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

Atropos is currently published as a pre-print on PeerJ, and will be submitted for peer review shortly. For now, you can cite it as:

> John P Didion, Marcel Martin, and Francis S Collins. "Atropos: specific, sensitive, and speedy trimming of sequencing reads." https://peerj.com/preprints/2452/
