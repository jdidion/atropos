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
    python summarize_timing_info.py -i ../results/timing_local_t4.txt --output-format latex \
      -o ../results/timing_local_table.latex --table-name "local-timing" \
      --table-caption "Execution time for programs running on desktop with 4 threads."
else
    # summarize timing
    cat commands_t*.e* | python summarize_timing_info.py --output-format latex \
      -o ../results/timing_cluster_table.latex --table-name "cluster-timing" \
      --table-caption "Execution time for programs running on a cluster."
fi
