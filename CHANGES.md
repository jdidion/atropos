# Changes

v1.0.7 (8/18/16)
----------------

* Re-engineered modifiers.py (and all dependent code) to enable use of modifiers that simultaneously edit both reads in a pair.
* Add --op-order option to enable use to specify order of first four trimming operations.
* Implemented insert-based alignment for paired-end adapter trimming. This is currently experimental. Benchmarking against SeqPurge and Skewer using simulated reads showed that the method Cutadapt uses to align adapters, while optimal for single-end reads, is much less sensitive and specific than the insert match algorithms used by SeqPurge and Skewer. Our algorithm is similar to the one used by SeqPurge but leverages the dynamic programming model of Cutadapt.

v1.0.6 (8/8/16)
---------------

* Based on tests, worker compression is faster than writer compression when more than 8 threads are available, so set this to be the default.

v1.0.5 (8/6/16)
---------------

* Interanal code reorganization - compression code moved to separate module
* Eliminated the --worker-compression option in favor of --compression (whose value is either 'worker' or 'writer')
* More documentation improvements

v1.0.3 (8/5/16)
---------------

* Significant performance improvements:
    * Start an extra worker once the main process is finished loading reads
    * Use system-level gzip for writer compression
    * Use writer compression by default
* More documentation fixes
* Disable quality trimming if all cutoffs are set to 0
* Eliminated the --parallel-environment option

v1.0.1 (8/4/16)
---------------

* Fix documentation bugs associated with migration from optparse to argparse

v1.0 (7/29/16)
--------------
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
