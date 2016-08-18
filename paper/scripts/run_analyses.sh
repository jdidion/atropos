#!/bin/bash
# Generate analysis commands and run.
# Execution depends on the environment (local or cluster).

script_dir=`pwd`
ATROPOS_ROOT=`dirname $script_dir`
ATROPOS_RESULT='results'
GB_PER_PROCESS=4

mkdir -p $ATROPOS_RESULT

while getopts "e:" opt; do
    case "$opt" in
    t)
        env=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

for threads in 4 8 16
do
    rm -f ./commands_t${threads}.sh
    ./prepare_analyses.sh -t $threads -r $ATROPOS_ROOT -o $ATROPOS_RESULT
    
    if [ "$env" == "local" ]
    then
        ./commands_t${threads}.sh 2>> timing_log.txt
    else
        swarm --threads-per-process ${threads} --gb-per-process $GB_PER_PROCESS --file commands_t${threads}.sh
    fi
done
