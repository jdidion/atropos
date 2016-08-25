#!/bin/bash
# Generate analysis commands and run.
# Execution depends on the environment (local or cluster).

script_dir=`pwd`
ATROPOS_ROOT=`dirname $script_dir`
ATROPOS_RESULT='results'
GB_PER_PROCESS=4

env='local'

while getopts "e:" opt; do
    case "$opt" in
    e)
        env=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [ "$env" == "local" ]
then
    # summarize timing
    python summarize_timing_info.py -i ../results/timing_local_t4.txt --output-format latex \
      -o ../results/timing_local_table.latex --table-name "local-timing" \
      --table-caption "Execution time for programs running on desktop with 4 threads."
    # summarize simulated accuracy (we only do this locally since, theoretically,
    # the results should be the same as when run on the cluster).
    mkdir ../results/simulated_accuracy
    for err in 001 005 01
    do
      for qcut in 0 20
      do
        i=1
        for read in 1 2
        do
          for file in in \
            atropos_4_${err}_q${qcut}_adapter_writercomp.${read}.fq.gz \
            atropos_4_${err}_q${qcut}_insert_writercomp.${read}.fq.gz \
            seqpurge_4_${err}_q${qcut}.${read}.fq.gz \
            skewer_4_${err}_q${qcut}-trimmed-pair${read}.fastq.gz
          do
            files[$i] = $file
            i=$((i+1))
          done
        done
        python summarize_simulated_trimming_accuracy.py \
          -a1 ../data/simulated/sim_${err}.1.aln
          -a2 ../data/simulated/sim_${err}.2.aln \
          -r1 ${file[$i]} -r2 ${file[$((i+4))]} \
          -o ../results/
          -s ../results/
else
    # summarize timing
    cat commands_t*.e* | python summarize_timing_info.py --output-format latex \
      -o ../results/timing_cluster_table.latex --table-name "cluster-timing" \
      --table-caption "Execution time for programs running on a cluster."
    # summarize accuracy on real data (we only ran this on the cluster)
    
fi
