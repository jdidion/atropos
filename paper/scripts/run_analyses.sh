#!/bin/bash
# Generate analysis commands and run.
# Execution depends on the environment (local or cluster).

script_dir=`pwd`
ATROPOS_ROOT=`dirname $script_dir`
ATROPOS_RESULT='results'
GB_PER_PROCESS=4
ALIGN_GB_PER_PROCESS=64
SORT_GB_PER_PROCESS=8
OVERLAP_GB_PER_PROCESS=16

export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ATROPOS_ROOT/software/bin

mkdir -p $ATROPOS_RESULT

env='local'
thread_list='4 8 16'
mode='dry'

while getopts "e:t:o:" opt; do
    case "$opt" in
    t)
        thread_list=$OPTARG
        ;;
    e)
        env=$OPTARG
        ;;
    o)
        mode=$OPTARG
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
    
    if [ "$mode" == "run"]
    then
        if [ "$env" == "local" ]
        then
            rm -f timing_log_${threads}.txt
            ./commands_t${threads}_shuf.sh 2>> ../results/timing_log_${threads}.txt
            ./rename_outputs.sh
            ./align_commands_t${threads}.sh
            ./sort_commands_t${threads}.sh
            ./bedops_commands_t${threads}.sh
            ./summarize.sh
        else
            trimJID=`swarm --jobid --threads-per-process ${threads} --gb-per-process $GB_PER_PROCESS --file commands_t${threads}.sh`
            # rename skewer outputs
            renameJID=`qsub -hold_jid $trimJID rename_outputs.sh`
            # map reads
            alignJID=`swarm --hold_jid $renameJID --threads-per-process ${threads} --gb-per-process $ALIGN_GB_PER_PROCESS --file align_commands_t${threads}.sh`
            # name-sort reads
            sortJID=`swarm --hold_jid $alignJID --threads-per-process ${threads} --gb-per-process $SORT_GB_PER_PROCESS --file sort_commands_t${threads}.sh`
            # overlap RNA-seq alignments with GENCODE annotations
            overlapJID=`swarm --hold_jid $sortJID --gb-per-process $OVERLAP_GB_PER_PROCESS --file bedops_commands_t${threads}.sh`
            # summarize trimming accuracy
            swarm --hold_jid $overlapJID --file summarize_commands.sh
        fi
    fi
done
