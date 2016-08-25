#!/bin/bash
# Generate analysis commands and run.
# Execution depends on the environment (local or cluster).

script_dir=`pwd`
ATROPOS_ROOT=`dirname $script_dir`
ATROPOS_RESULT='results'
GB_PER_PROCESS=4

export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ATROPOS_ROOT/software/bin

mkdir -p $ATROPOS_RESULT

env='local'
thread_list='4 8 16'

while getopts "e:t:" opt; do
    case "$opt" in
    t)
        thread_list=$OPTARG
        ;;
    e)
        env=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

for threads in $thread_list
do
    rm -f ./commands_t${threads}.sh ./commands_t${threads}_shuf.sh
    ./prepare_analyses.sh -t $threads -r $ATROPOS_ROOT -o $ATROPOS_RESULT
    # shuffle the commands just to make sure there's no bias associated with the ordering
    shuf -o ./commands_t${threads}_shuf.sh ./commands_t${threads}.sh
    
    if [ "$env" == "local" ]
    then
        rm -f timing_log_${threads}.txt
        ./commands_t${threads}.sh 2>> ../results/timing_log_${threads}.txt
        ./rename_outputs.sh
    else
        swarm --threads-per-process ${threads} --gb-per-process $GB_PER_PROCESS \
          --file commands_t${threads}.sh
        # rename skewer outputs
        qsub <dependency> rename_outputs.sh
        # map reads
        swarm <dependency> --threads-per-process ${threads} \
          --gb-per-process $BWAMETH_GB_PER_PROCESS \
          --file bwameth_commands.sh
    fi
done
