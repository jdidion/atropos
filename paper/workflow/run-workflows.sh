#!/bin/bash
# Run the nextflow workflows. Usage: ./run-workflows.sh <local|cluster>
PROFILE=$1
NEXTFLOW=nextflow # change to point to nextflow executable if not in path
$NEXTFLOW run -c nextflow.config -profile local simulated.nf
$NEXTFLOW run -c nextflow.config -profile local rnaseq.nf
$NEXTFLOW run -c nextflow.config -profile local wgbs.nf