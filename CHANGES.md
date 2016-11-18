# Changes

v1.0.16 (dev)
-------------

* Migrate to Versioneer for version management.
* Enable stderr as a valid output using the '_' shortcut.
* Add ability to specify SAM/BAM as input format.
* Add option to select which read to use when treating a paired-end interleaved or SAM/BAM file as single-end.
* Remove restrictions on the use of --merge-overlapping, and enable error correction during merging.
* We are beginning to move towards the use of commands for all operations, and read-trimming now falls under the 'trim' command. Currently, calling atropos without any command will default to the 'trim' command.
* When InsertAdapterCutter.symmetric is True and mismatch_action is not None, insert match fails, at least one adapter match succeeds, and the adapter matches (if there are two) are complementary, then the reads are treated as overlapping and error correction is performed. This leads to substantial improvements when one read is of good quality while the other is other is of poor quality.

v1.0.15 (2016.09.14)
--------------------

* Fixed missing import bug in 'detect' command.
* Added estimate of fraction of contaminated reads to output of 'detect' command.
* Optionally cache list of known contaminants rather than re-download it every time.

v1.0.14 (2016.09.13)
--------------------

* Implemented \_align.MultiAligner, which returns all matches that satisfy the overlap and error thresholds. align.InsertAligner now uses MultiAligner for insert matching, and tests all matches in decreasing size order until it finds one with adapter matches (if any).
* Major improvements to the accuracy of the 'detect' command.
* Added options for how to correct mismatched bases for which qualities are equal.
* Added option to select a single pair of adapters from multiple sequences in a fasta file.
* Fixed report when insert-match is used.
* Fixed several bugs when using the "message" progress bar (thanks to Thomas Cokelaer!).
* Fixed a segmentation fault that occurs when trying to trim zero-length reads with the insert aligner.
* Sevaral other bugfixes.

v1.0.13 (2016.08.31)
--------------------

* Add options to specify max error rates for insert and adapter matching within insert aligner.
* Add new command to estimate empirical error rate in data set from base qualities.

v1.0.12 (2016.08.30)
--------------------

* Add ability to correct errors during insert-match adapter trimming.
* Implement additional adapter-detection algorithms.
* Fix bug where default output file is force-created in parallel-write mode

v1.0.11 (2016.08.24)
--------------------

* Clarify and fix issues with bisulfite trimming. Notably, rrbs and non-directional are now allowed independently or in combination.

v1.0.10 (2016.08.23)
--------------------

* Introduced new 'detect' command for automatically detecting adapter sequences.
* Options are now required to specify input files.
* Major updates to documentation.

v1.0.9 (2016.08.22)
-------------------

* Bugfix release

v1.0.8 (2016.08.19)
-------------------

* Reverted previously introduced (and no longer necessary) dependency on bitarray).
* Switched the insert aligner back to the default implementation, as the one that ignores indels is not any faster.

v1.0.7 (2016.08.18)
-------------------

* Re-engineered modifiers.py (and all dependent code) to enable use of modifiers that simultaneously edit both reads in a pair.
* Add --op-order option to enable use to specify order of first four trimming operations.
* Implemented insert-based alignment for paired-end adapter trimming. This is currently experimental. Benchmarking against SeqPurge and Skewer using simulated reads showed that the method Cutadapt uses to align adapters, while optimal for single-end reads, is much less sensitive and specific than the insert match algorithms used by SeqPurge and Skewer. Our algorithm is similar to the one used by SeqPurge but leverages the dynamic programming model of Cutadapt.

v1.0.6 (2016.08.08)
-------------------

* Based on tests, worker compression is faster than writer compression when more than 8 threads are available, so set this to be the default.

v1.0.5 (2016.08.06)
-------------------

* Interanal code reorganization - compression code moved to separate module
* Eliminated the --worker-compression option in favor of --compression (whose value is either 'worker' or 'writer')
* More documentation improvements

v1.0.3 (2016.08.05)
-------------------

* Significant performance improvements:
    * Start an extra worker once the main process is finished loading reads
    * Use system-level gzip for writer compression
    * Use writer compression by default
* More documentation fixes
* Disable quality trimming if all cutoffs are set to 0
* Eliminated the --parallel-environment option

v1.0.1 (2016.08.04)
-------------------

* Fix documentation bugs associated with migration from optparse to argparse

v1.0 (2016.07.29)
-----------------
* Initial release (forked from cutadapt 1.10)
* Re-wrote much of filters.py and modifiers.py to separate modifying/filtering from file writing.
    * File writing is now managed by a separate class (seqio.Writers)
    * There are container classes for managing filters (filters.Filters) and modifiers (modifiers.Modifiers)
* Re-wrote all of the output-oriented code in seqio.py
    * Formatting Sequence objects is now separate from writing data
    * There is a container class (seqio.Formatters) that manages the formatters for output files
    * Added support for interleaved output
* Implemented multiprocessing
    * Added several new options in scripts.atropos to control parallelization
    * Wrote all of the parallel processing code in atropos.multicore
    * Renamed scripts.atropos.process_single_reads() to scripts.atropos.run_serial() and rewrote to work similarly to atropos.multicore.run_parallel()
    * Added ability to merge report statistics from multiple worker threads
* Added miRNA and bisulfite sequencing options to scripts.atropos
* Added progress bar support
* Switched argument parsing to argparse
* Reorganized the monolithic scripts.atropos.main() into multiple functions
* Dropped all support for python 2.x
