.. image:: https://travis-ci.org/jdidion/atropos.svg?branch=master
    :target: https://travis-ci.org/marcelm/cutadapt

.. image:: https://img.shields.io/pypi/v/atropos.svg?branch=master
    :target: https://pypi.python.org/pypi/atropos

=======
atropos
=======

Atropos is tool for specific, sensitive, and speedy trimming of NGS reads. It is a fork of the
venerable `Cutadapt <https://github.com/marcelm/cutadapt>`_ read trimmer, with the primary
improvements being:

1. Multi-threading support that has been optimized for different parallel computing environments.
2. Options for trimming specific types of data (miRNA, bisulfite).
3. The ability (currently limited) to merge overlapping reads.
4. The ability to write the summary report and log messages to separate files.
5. The ability to write interleaved FASTQ output.
6. A progress bar, and other minor usability enhancements.

Also note that one potential downside of Atropos is that it requires python 3.3+. If you require
python 2.7 support, you can still use `Cutadapt <https://github.com/marcelm/cutadapt>`_.

Atropos finds and removes adapter sequences, primers, poly-A tails and other
types of unwanted sequence from your high-throughput sequencing reads.

Cleaning your data in this way is often required: Reads from small-RNA
sequencing contain the 3’ sequencing adapter because the read is longer than
the molecule that is sequenced. Amplicon reads start with a primer sequence.
Poly-A tails are useful for pulling out RNA from your sample, but often you
don’t want them to be in your reads.

Atropos helps with these trimming tasks by finding the adapter or primer
sequences in an error-tolerant way. It can also modify and filter reads in
various ways. Adapter sequences can contain IUPAC wildcard characters. Also,
paired-end reads and even colorspace data is supported. If you want, you can
also just demultiplex your input data, without removing adapter sequences at all.

The parts of Atropos that were carried over from Cutadapt are under an MIT license,
while improvements and additions are in the public domain.

The citation for the original Cutadapt paper is:
 
  Marcel Martin. "Cutadapt removes adapter sequences from high-throughput sequencing reads."
  EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

A manuscript for Atropos is currently in preparation. For now, you can cite it as:

  John P Didion. "Atropos." 2016. https://github.com/jdidion/atropos

Links
-----

* `Documentation <https://atropos.readthedocs.org/>`_
* `Source code <https://github.com/jdidion/atropos/>`_
* `Report an issue <https://github.com/jdidion/atropos/issues>`_
* `Project page on PyPI (Python package index) <https://pypi.python.org/pypi/atropos/>`_
