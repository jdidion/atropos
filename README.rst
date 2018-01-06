|Travis CI| |PyPi| |DOI|

Atropos
=======

Atropos is tool for specific, sensitive, and speedy trimming of NGS
reads. It is a fork of the venerable Cutadapt read trimmer
(https://github.com/marcelm/cutadapt,
`DOI:10.14806/ej.17.1.200 <http://dx.doi.org/10.14806/ej.17.1.200>`__),
with the primary improvements being:

1.  Multi-threading support, including an extremely fast "parallel
    write" mode.
2.  Implementation of a new insert alignment-based trimming algorithm
    for paired-end reads that is substantially more sensitive and
    specific than the original Cutadapt adapter alignment-based
    algorithm. This algorithm can also correct mismatches between the
    overlapping portions of the reads.
3.  Options for trimming specific types of data (miRNA, bisulfite-seq).
4.  A new command ('detect') that will detect adapter sequences and
    other potential contaminants.
5.  A new command ('error') that will estimate the sequencing error
    rate, which helps to select the appropriate adapter- and quality-
    trimming parameter values.
6.  A new command ('qc') that generates read statistics similar to
    FastQC. The trim command can also compute read statistics both
    before and after trimming (using the '--stats' option).
7.  Improved summary reports, including support for serialization
    formats (JSON, YAML, pickle), support for user-defined templates
    (via the optional Jinja2 dependency), and integration with
    `MultiQC <http://multiqc.info>`__.
8.  The ability to merge overlapping reads (this is experimental and the
    functionality is limited).
9.  The ability to write the summary report and log messages to separate
    files.
10. The ability to read SAM/BAM files and read/write interleaved FASTQ
    files.
11. A progress bar, and other minor usability enhancements.

Manual installation
-------------------

Atropos is available from
`pypi <https://pypi.python.org/pypi/atropos>`__ and can be installed
using ``pip``.

First install dependencies:

-  Required

   -  Python 3.3+ (python 2.x is NOT supported) - note: we have
      identified a possible bug in python 3.4.2 that causes random
      segmentation faults. We think this mainly affects unit testing
      (and thus specifically test on 3.4.3). If you encounter this bug,
      we recommend upgrading to a newer python version.
   -  Cython 0.25.2+ (``pip install Cython``)

-  Maybe

   -  pytest (for running unit tests)
   -  progressbar2 or tqdm (progressbar support)
   -  pysam (SAM/BAM input)
   -  khmer 2.0+ (``pip install khmer``) (for detecting low-frequency
      adapter contamination)
   -  jinja2 (for user-defined report formats)

Then run:

``pip install atropos``

Conda
-----

There is an Atropos recipe in
`Bioconda <https://anaconda.org/bioconda/atropos>`__.

``conda install -c bioconda atropos``

Docker
------

A `Docker image <https://hub.docker.com/r/jdidion/atropos/>`__ is
available for Atropos in Docker Hub.

``docker run jdidion/atropos <arguments>``

Usage
-----

Atropos is almost fully backward-compatible with cutadapt. If you
currently use cutadapt, you can simply install Atropos and then
substitute the executable name in your command line, with one key
difference: you need to use options to specify input file names. For
example:

::

    atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz -se reads.fq.gz

To take advantage of multi-threading, set the ``--threads`` option:

::

    atropos --threads 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA -o trimmed.fq.gz -se reads.fq.gz

To take advantage of the new aligner (if you have paired-end reads with
3' adatpers), set the ``--aligner`` option to 'insert':

::

    atropos --aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o trimmed.1.fq.gz -p trimmed.2.fq.gz \
      -pe1 reads.1.fq.gz -pe2 reads.2.fq.gz

See the `Documentation <https://atropos.readthedocs.org/>`__ for more
complete usage information.

Publication
-----------

A `preprint <https://peerj.com/preprints/2452/>`__ is available and the
submitted paper is currently under review. The results in the paper can
be fully reproduced using the workflow defined in the
`paper <paper/README.md>`__ directory.

Developers
----------

We welcome any contributions via GitHub issues and pull requests. See
the `documentation <https://atropos.readthedocs.org/>`__ for style
guidelines and best practices. We enforce the `Contributor
Covenant <http://contributor-covenant.org/>`__ code of conduct.

Links
-----

-  `Documentation <https://atropos.readthedocs.org/>`__
-  `Source code <https://github.com/jdidion/atropos/>`__
-  `Report an issue <https://github.com/jdidion/atropos/issues>`__

Roadmap
-------

1.2
~~~

-  Migrate to xphyle (https://github.com/jdidion/xphyle) for file
   management.
-  Provide option for RNA-seq data that will trim polyA sequence.
-  Accept multiple input files.
-  Support SAM output.
-  Expand the list of contaminants that are detected by default.
-  Accessibility:

   -  Create recipe for homebrew.
   -  Automatically update conda and homebrew recipes for each release.
   -  Create Galaxy tool description using
      `argparse2tool <https://github.com/erasche/argparse2tool#cwl-specific-functionality>`__.

1.3
~~~

-  Provide PacBio- and nanopore-specific options
   (https://github.com/marcelm/cutadapt/issues/120).
-  Currently, InsertAligner requires a single 3' adapter for each end.
   Adapter trimming will later be generalized so that A) the
   InsertAligner can handle multiple matched pairs of adapters and/or B)
   multiple different aligners can be used for different adapters.
-  Automate creation and sending of user statistics and crash reports
   using `pytattle <https://github.com/biologyguy/PyTattle>`__.

1.4
~~~

-  Migrate to seqio (https://github.com/jdidion/seqio) for
   reading/writing sequence files.
-  General-purpose read filtering based on read ID:
   https://github.com/marcelm/cutadapt/issues/107.
-  Currently, SAM/BAM input files must be name sorted; add an option to
   1) pre-sort reads inline using samtools or sambamba, or 2) cache each
   read in memory until its mate is found.

1.5
~~~

-  Provide more user control over anchoring of adapters:
   https://github.com/marcelm/cutadapt/issues/53.
-  Support for paired-end demultiplexing (i.e. when barcodes are used in
   both paired-end adapters):
   https://github.com/marcelm/cutadapt/issues/118.
-  Add option to estimate bisulfite conversion rate from filled-in
   cytosine methylation status in reads that were MspI-digested.
-  Consider supporting different error rates for read1 vs read2.
-  Add a ClipOverlapping modifier that will remove read overlaps (as
   opposed to merging).
-  Add option to InsertAdapter to trim overhangs without adapter
   matching.

1.6
~~~

-  Implement a public plugin API.
-  Add more logging and convert log messages from old-style to new-style
   format strings.

2.0
~~~

-  Simplification of command line options, perhaps by further splitting
   functionality up into different sub-commands, but also by more
   intelligent choices for default option values based on context.
-  Consider adding additional report formats

   -  https://github.com/marcelm/cutadapt/issues/112

-  Performance enhancements using

   -  http://numba.pydata.org/
   -  https://github.com/undefx/vecpy
   -  https://github.com/serge-sans-paille/pythran

-  90% test coverage

-  Fuzz testing using AFL

   -  http://lcamtuf.coredump.cx/afl/
   -  https://github.com/jwilk/python-afl

Beyond 2.0
~~~~~~~~~~

-  Implement additional alternate alignment algorithms.
-  Implement the error detection algorithm in ADEPT:
   https://github.com/LANL-Bioinformatics/ADEPT
-  Implement the quality trimming algorithm used in UrQt:
   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450468/
-  Scythe is an interesting new trimmer. Depending on how the benchmarks
   look in the forthcoming paper, we will add it to the list of tools we
   compare against Atropos, and perhaps implement their Bayesian
   approach for adapter match.
-  Experiment with replacing the multicore implementation with an
   asyncio-based implementation (using ProcessPoolExecutor and uvloop).

While we consider the command-line interface to be stable, the internal
code organization of Atropos is likely to change. At this time, we
recommend to not directly interface with Atropos as a library (or to be
prepared for your code to break). The internal code organization will be
stabilized as of version 2.0, which is planned for sometime in 2017.

If you would like to suggest additional enhancements, you can submit
issues and/or pull requests at our GitHub page.

Citations
---------

The citation for the original Cutadapt paper is:

    Marcel Martin. "Cutadapt removes adapter sequences from
    high-throughput sequencing reads." EMBnet.Journal, 17(1):10-12, May
    2011. http://dx.doi.org/10.14806/ej.17.1.200

Atropos is currently published as a pre-print on PeerJ, and will be
submitted for peer review shortly. For now, you can cite it as:

    John P Didion, Marcel Martin, and Francis S Collins. "Atropos:
    specific, sensitive, and speedy trimming of sequencing reads."
    https://peerj.com/preprints/2452/

.. |Travis CI| image:: https://travis-ci.org/jdidion/atropos.svg?branch=master
   :target: https://travis-ci.org/jdidion/atropos]
.. |PyPi| image:: https://img.shields.io/pypi/v/atropos.svg?branch=master
   :target: https://pypi.python.org/pypi/atropos
.. |DOI| image:: https://zenodo.org/badge/61393086.svg
   :target: https://zenodo.org/badge/latestdoi/61393086
