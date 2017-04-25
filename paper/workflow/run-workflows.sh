#!/bin/bash
# Run the nextflow workflows. Usage: ./run-workflows.sh <local|cluster>
PROFILE=$1
NEXTFLOW=nextflow # change to point to nextflow executable if not in path
$NEXTFLOW run -c nextflow.config -profile $PROFILE simulated.nf
$NEXTFLOW run -c nextflow.config -profile $PROFILE rnaseq.nf
$NEXTFLOW run -c nextflow.config -profile $PROFILE wgbs.nf
