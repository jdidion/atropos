#!/bin/bash
# This script will run the analyses for the Atropos paper, as follows:
#
# For each error profile, we run Atropos, SeqPurge, and Skewer using
# equivalent (or as similar as possible) arguments. The tests are run
# both locally (Late 2013 Mac Pro, 3.7 GHz quad-core Xeon E5, 32 GB
# memory) and on a cluster (SL6, ???, connected to Isilon NAS over
# InfiniBand??). Atropos is run both with and without a separate Writer
# process. The tests are run on 1) simulated data (see simulate_reads.sh)
# and 2) real data (from the SeqPurge paper, Sturm et al. 2016,
# ftp://ftp.sra.ebi.ac.uk/vol1/ERA494/ERA494451/).

# Call: run_analyses.sh -m <local|cluster> -t <threads>

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
mode="local"
threads=8

while getopts "m:t:" opt; do
    case "$opt" in
    m)
        mode=$OPTARG
        ;;
    t)  threads=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

echo "Running with arguments mode: $mode, threads: $threads, unused args: $@"

if [ $mode == local ]
then
    ../bin/atropos
else
    
fi
