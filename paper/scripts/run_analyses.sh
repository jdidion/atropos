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
    rm -f ./commands_t${threads}.sh
    ./prepare_analyses.sh -t $threads -r $ATROPOS_ROOT -o $ATROPOS_RESULT
    
    if [ "$env" == "local" ]
    then
        ./commands_t${threads}.sh 2>> timing_log.txt
    else
        swarm --threads-per-process ${threads} --gb-per-process $GB_PER_PROCESS --file commands_t${threads}.sh
    fi
done
