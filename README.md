![https://travis-ci.org/jdidion/atropos](https://travis-ci.org/jdidion/atropos.svg?branch=master)
![https://pypi.python.org/pypi/atropos](https://img.shields.io/pypi/v/atropos.svg?branch=master)

# Atropos

Atropos is tool for specific, sensitive, and speedy trimming of NGS reads. It is a fork of the venerable Cutadapt read trimmer (https://github.com/marcelm/cutadapt, [DOI:10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)), with the primary improvements being:

1. Multi-threading support, including an extremely fast "parallel write" mode.
2. Implementation of a new insert alignment-based trimming algorith for paired-end reads that is substantially more sensitive and specific than the original Cutadapt adapter alignment-based algorithm.
3. Options for trimming specific types of data (miRNA, bisulfite-seq).
4. The ability (currently limited) to merge overlapping reads.
5. The ability to write the summary report and log messages to separate files.
6. The ability to write interleaved FASTQ output.
7. A progress bar, and other minor usability enhancements.

## Dependencies

* Python 3.3+ (python 2.x is NOT supported)
* Cython 0.24+ (`pip install Cython`)
* progressbar or tqdm (optional, if you want progressbar support)

## Installation

`pip install atropos`

## Usage

Atropos is fully backward-compatible with cutadapt. If you currently use cutadapt, you can simply install Atropos and then substitute the executable name in your command line. For example:

```{python}
atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz reads.fq.gz
```

To take advantage of multi-threading, set the `--threads` option:

```{python}
atropos --threads 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz reads.fq.gz
```

To take advantage of the new aligner (if you have paired-end reads with 3' adatpers), set the `--aligner` option to 'insert':

```{python}
atropos --aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed.1.fq.gz -p trimmed.2.fq.gz \
  reads.1.fq.gz reads.2.fq.gz
```

See the [Documentation](https://atropos.readthedocs.org/) for more complete usage information.

## Links

* [Documentation](https://atropos.readthedocs.org/)
* [Source code](https://github.com/jdidion/atropos/)
* [Report an issue](https://github.com/jdidion/atropos/issues)

## Planned enhancements and experiments

* Implement an auto-detect option for adapters similar to TrimGalore: read the first 1M reads, search for common adapters, and pick the one that appears most frequently.
* Optional error correction for overlapping pairs (if there are mismatches, change the lower quality base to the higher quality base)
* Currently, InsertAligner requires a single 3' adapter for each end. Adapter trimming will later be generalized so that A) the InsertAligner can handle multiple matched pairs of adapters and/or B) multiple different aligners can be used for different adapters.

## Citations

The citation for the original Cutadapt paper is:
 
> Marcel Martin. "Cutadapt removes adapter sequences from high-throughput sequencing reads." EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

A manuscript for Atropos is currently in preparation. For now, you can cite it as:

> John P Didion. "Atropos." 2016. https://github.com/jdidion/atropos
