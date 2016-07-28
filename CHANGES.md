# Changes

v1.0
----
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
