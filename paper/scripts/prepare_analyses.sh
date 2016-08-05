#!/bin/bash

# This script will generate the commands to run the analyses for the
# Atropos paper.
#
# For each error profile, we run Atropos, SeqPurge, and Skewer using
# equivalent (or as similar as possible) arguments. The tests are run
# both locally (Late 2013 Mac Pro, 3.7 GHz quad-core Xeon E5, 32 GB
# memory) and on a cluster (SL6, ???, connected to Isilon NAS over
# InfiniBand??). Atropos is run both with and without a separate Writer
# process. The tests are run on 1) simulated data (see simulate_reads.sh)
# and 2) real data (from the SeqPurge paper, Sturm et al. 2016,
# ftp://ftp.sra.ebi.ac.uk/vol1/ERA494/ERA494451/).
#
# Call: run_analyses.sh -m <local|cluster> -t <threads> -r <root dir>

# A POSIX variable; reset in case getopts has been used
# previously in the shell.
OPTIND=1

# Set default values
mode="local"
threads=8
root="."

while getopts "m:t:o:r:" opt; do
    case "$opt" in
    m)
        mode=$OPTARG
        ;;
    t)
        threads=$OPTARG
        ;;
    r)
        root=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

commands="${mode}_commands"

echo "Running with arguments " \
     "mode: $mode, threads: $threads, command file: $commands, root: $root" \
     "unused args: $@"

## Constants

# binaries
ATROPOS=atropos
SEQPURGE=$root/software/bin/SeqPurge
SKEWER=$root/software/bin/skewer
# adapter sequences
ADAPTER1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
# minimum read length after trimming
MIN_LEN=25
# minimum number of adapter bases that must overlap
MIN_OVERLAP=7
# number of reads to process in a batch
# (also used as prefetch size for SeqPurge)
BATCH_SIZE=5000

for prog in atropos seqpurge skewer
do
  mkdir $prog
done

## simulated reads

for err in 001 005 01
do
  fq1=$root/data/simulated/sim_${err}.1.fq
  fq2=$root/data/simulated/sim_${err}.2.fq
  
  atropos_command="$ATROPOS " \
  "-a $ADAPTER1 -A $ADAPTER2"
  
  if [ $mode == local ]
  then
    echo "$atropos_command --parallel-environment local " \
    "$fq1 $fq2" >> $commands
  else
    echo "$atropos_command --parallel-environment cluster " \
    "$fq1 $fq2" >> $commands
    echo "$atropos_command --parallel-environment cluster " \
    "--no-writer-process $fq1 $fq2" >> $commands
  fi

  # NOTE: SeqPurge provides two parameters, -match_perc and
  # -mep, which may be important for tuning performance;
  # however, the use of these parameters is unclear from the
  # command-line help, and no other guidance is provided.
  # Furthermore, it is not stated in the SeqPurge paper whether
  # these parameters were modified from their defaults (e.g. if
  # -mep is supposed to match the exected sequence error rate).
  # Thus, we leave these parameters set to their default values.
  for qcut in 0 20
  do
    if [ "$qcut" -eq "0" ]
    then
      n=''
      ncut=0
    else
      n='-n'
      ncut=7
    fi
    profile="${mode}_${err}_q${qcut}"

    echo "$SEQPURGE -in1 $fq1 -in2 $fq2" \
    "-out1 seqpurge/sim_${profile}.1.fq.gz" \
    "-out2 seqpurge/sim_${profile}.2.fq.gz" \
    "-a1 $ADAPTER1 -a2 $ADAPTER2 -threads $threads" \
    "-qcut $qcut -ncut $ncut -min_len $MIN_LEN" \
    "-prefetch $BATCH_SIZE" \
    "-summary seqpurge/sim_${profile}.summary" >> $commands
  done

  echo "$SKEWER -m pe -l $MIN_LEN -k $MIN_OVERLAP" \
  "-o skewer/sim_${profile} -z" \
  "-x $ADAPTER1 -y $ADAPTER2 -t $threads" \
  "-q $qcut $n $fq1 $fq2" >> $commands
done
