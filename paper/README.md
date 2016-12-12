The scripts in this directory should enable you to re-run the analyses in the Atropos paper.

# 1. Install software

The scripts/install_software.sh script installs all the software in the 'software' directory into the software/bin directory. It also installs
Atropos from pypi. Finally, it builds the STAR and bwa-meth indexes if they don't already exist. You may need to edit this script to set paths on your local environment.

For the script to work, you must first have the following pre-requisites. We recommend using conda to install all of these packages. Conda is available [stand-alone](http://conda.pydata.org/miniconda.html) or as part of the [Anaconda python distribution](https://www.continuum.io/downloads).

1. modern C++ compiler (we use gcc 5.1)
2. python 3.3+
    * cython 0.25+
    * mako
3. To install ART you need:
    * automake
4. To install SeqPurge, you need
    * Qt 5.3+ with xmlpatterns and mysql packages
    * git
    * cmake
5. To run the analyses on real data, you need the following, which are all
available in the bioconda repository: https://anaconda.org/bioconda):
    * BWA (the MEM command is required; we used version 0.7.15)
    * bwa-meth
    * STAR aligner
    * samtools
    * bedops
    * SRA tools
6. Download the reference sequence and annotations and place them in $root/data
    * Reference genome NCBI GRCh37:
     ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3
        * If you download individual chromsome files, concatenate them all together into a single fasta
    * ENCODE v19 Annotations:
     ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

 You can create and use a new virtual environment for testing purposes - and this is required if your normal environment has python < 3.3:

     conda create --name atropos_test python=3.5

To run the script:

    cd scripts
    ./install_software.sh

# 2. Simulate reads

There are already 3 simulated datasets in the 'data' directory. If you'd like to re-create these, run:
  
    ./simulate_reads.sh

# 3. Prepare the command scripts

    for threads in 4 8 16
    do
         ./prepare_analyses.sh -t $threads -o <mode>
    done

Where <mode> is either 'local' or 'cluster'. For each number of threads, this will generate several command scripts, and a main script, named run_t{threads}_{mode}, that will execute all the necessary steps when run. This script will also download the real datasets if they don't exist already.

To submit jobs on the cluster, we use the 'swarm' script written by Peter Chines (NHGRI), which is specific to Sun Grid Engine. Similar tools should be available for other queue management systems.

# TODO

* These papers do a nice job of benchmarking trimmers. Consider adding some more benchmarks.
    * http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-16-S1-S2
    * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0454-y
