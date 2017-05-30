#!/bin/bash
# Run the nextflow workflows. Usage: 
# ./run-workflows.sh -p <local|cluster> [-w <workdir>] [-n <nextflow exe>]

OPTIND=1 # Reset in case getopts has been used previously in the shell.
NEXTFLOW=nextflow

# Initialize our own variables:
while getopts "p:w:n:" opt; do
    case "$opt" in
    p)  PROFILE=$OPTARG
        ;;
    w)  export NXF_HOME="$OPTARG/nf"
        export NXF_WORK="$OPTARG/work"
        export NXF_TEMP="$OPTARG/tmp"
        ;;
    n)  NEXTFLOW=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [ $PROFILE="cluster" ]
then
  # Load the Nextflow module; you may need to change this on your cluster
  module load nextflow/0.23.3
fi
$NEXTFLOW run -c nextflow.config -profile $PROFILE simulated.nf
$NEXTFLOW run -c nextflow.config -profile $PROFILE rnaseq.nf
$NEXTFLOW run -c nextflow.config -profile $PROFILE wgbs.nf
