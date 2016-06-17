![https://travis-ci.org/jdidion/atropos](https://travis-ci.org/jdidion/atropos.svg?branch=master)
![https://pypi.python.org/pypi/atropos](https://img.shields.io/pypi/v/atropos.svg?branch=master)

# Atropos

Atropos is tool for specific, sensitive, and speedy trimming of NGS reads. It is a fork of the venerable cutadapt read trimmer (https://github.com/marcelm/cutadapt, [DOI:10.14806/ej.17.1.200](http://dx.doi.org/10.14806/ej.17.1.200)), with the primary improvements being:

1. Multi-threading support that has been optimized for different parallel computing environments.
2. Options for trimming specific types of data (miRNA, bisulfite).
3. The ability (currently limited) to merge overlapping reads.
4. The ability to write the summary report and log messages to separate files.
5. The ability to write interleaved FASTQ output.
6. A progress bar, and other minor usability enhancements.

## Dependencies

* Python 2.7 or 3.3+
* Cython 0.24+ (`pip install Cython`)

## Installation

`pip install atropos`

## Usage

Atropos is fully backward-compatible with cutadapt. If you currently use cutadapt, you can simply install Atropos and then substitute the executable name in your command line. For example:

```{r}
atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz reads.fq.gz
```

To take advantage of multi-threading, set the `--threads` option:

```
{r}
atropos --threads 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz reads.fq.gz
```

See the [Documentation](https://atropos.readthedocs.org/) for more complete usage information.

## License

Atropos incorporates a substantial amount of code from cutadapt. Those portions of code are clearly marked, and are [governed by the MIT license](https://github.com/marcelm/cutadapt/blob/master/LICENSE)). The portions of code developed by John Didion (also clearly marked) are government work and thus in the public domain. Note that by contributing pull requests to this project you give up the right to copyright that code, as it will (if accepted) also belong to the public domain.

## Technical details

When the `--threads` option is set, the main process that loads batches of reads from the input file(s) and posts them to a Queue. One or more worker processes take batches from the Queue and process them in much the same manner as cutadapt. Writing the results to disk can be a bottleneck, and so there are multiple options for how to handle this:

1. Local mode: the program is being run on a system with a local disk (e.g. a laptop or desktop machine) with perhaps limited memory. In this mode, there is a worker thread that collects the results from the writer threads, compresses the data (if necessary), and writes it to disk. A relatively small maximum Queue size is used to keep from taking up too much memory.
2. Cluster mode: the program is being run on a system with ample memory, but writes to remote file storage (NAS), which typically has much higher latency than local storage. In this mode, larger maximum Queue sizes are used, and data compression is handled by the worker threads, so that the writer thread's only task is writing bytes to disk.
3. Parallel writing mode: in many cases, it is not actually necessary to write all results to the same file. For example, if the next processing step after trimming is alignment, and your aligner supports reading from multiple FASTQ files (or you are on a linux-based system and can use [process substitution](http://www.tldp.org/LDP/abs/html/process-sub.html)) to concatenate multiple files to a single input stream) then it can be much faster to have worker thread write results directly to separate files.

If you are going to be processing lots of data, we recommend taking some time to optimize Atropos for your particular environment. You can do this by using the `--max-reads` parameter to limit the number of input reads (we recommend 1-10 M reads to get a good idea of average processing time per read) and then experiment with multiple parameter combinations. Turning on DEBUG log messages (`--log-level DEBUG`) can also be helpful for this task. Note that increasing the number of available threads has diminishing returns, and should certainly not exceed the number of cores available on your system. We generally find 8 threads to offer the best trade-off between speed and resource usage, though this may differ for your own environment.

## The Parallel Pipeline Paridigm

We are currently improving Atropos to support the [Parallel Pipeline Paridigm]() (PPP). The general idea of PPP is that software pipelines waste too much time on splitting and merging data, and on disk I/O. The solution is to connect all components of the pipeline so that all data is transmitted in memory, and, furthermore, so that when adjacent tools in the pipeline are both multi-threaded, data can be transmitted directly from a parallel process in the first tool to a corresponding process in the second tool without the need for multiplexing and then demultiplexing the data. PPP advocates sockets for inter-process communication (IPC), but also allows for named pipes and other file-like methods of IPC, especially for integrating tools that do not yet support sockets. We are developing [Lachesis]() as a light-weight framework for creating and executing parallel pipleines.

## Links

* [Documentation](https://atropos.readthedocs.org/)
* [Source code](https://github.com/jdidion/atropos/)
* [Report an issue](https://github.com/jdidion/atropos/issues)
* [Project page on PyPI (Python package index)](https://pypi.python.org/pypi/atropos/)
