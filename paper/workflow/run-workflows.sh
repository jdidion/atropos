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
  # These statements are required on our cluster; YMMV
  ## Load the Nextflow module
  module load nextflow/0.23.3
  ## Unset the Singularity bindpath, which interferes with paths
  ## in our data containers.
  export SINGULARITY_BINDPATH=
fi
$NEXTFLOW run -c nextflow.config -profile $PROFILE simulated.nf
$NEXTFLOW run -c nextflow.config -profile $PROFILE rnaseq.nf
$NEXTFLOW run -c nextflow.config -profile $PROFILE wgbs.nf
