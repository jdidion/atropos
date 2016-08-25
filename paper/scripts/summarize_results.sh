#!/bin/bash
# Generate analysis commands and run.
# Execution depends on the environment (local or cluster).

script_dir=`pwd`
root=`dirname $script_dir`
outdir='results'
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
    mkdir -p $root/results/simulated_accuracy
    for err in 001 005 01
    do
        for profile in in \
          atropos_4_${err}_q0_adapter_writercomp \
          atropos_4_${err}_q0_insert_writercomp \
          seqpurge_4_${err}_q0 \
          skewer_4_${err}_q0-trimmed-pair
        do
            python summarize_simulated_trimming_accuracy.py \
              -a1 $root/data/simulated/sim_${err}.1.aln
              -a2 $root/data/simulated/sim_${err}.2.aln \
              -r1 $outdir/$profile.1.fq.gz -r2 $outdir/$profile.2.fq.gz \
              -o $root/results/simulated_accuracy/$profile.txt \
              -s $root/results/simulated_accuracy/$profile.summary.txt \
              -t $root/results/simulated_accuracy/table.txt \
              --name $profile
        done
    done
else
    # summarize timing
    cat commands_t*.e* | python summarize_timing_info.py --output-format latex \
      -o ../results/timing_cluster_table.latex --table-name "cluster-timing" \
      --table-caption "Execution time for programs running on a cluster."
    # summarize accuracy on real data (we only ran this on the cluster)
fi
