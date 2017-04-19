#!/bin/bash
scripts=`pwd`
root=`dirname $scripts`
mkdir $root/software/build
# Set this to the location of the hg19 reference multifasta.
# If this file doesn't exist, the bwameth index won't be built.
genome_dir=$root/data/reference
genome=$genome_dir/ref.fa
annotations=$genome_dir/gencode.v19.annotation.gtf

# TODO: switch to using Docker instances for all these software

# Install version of Atropos on which manuscript is based
pip install atropos==1.1.1

# Install Skewer from conda
conda install skewer

# TODO: install bwa, star

# Build the aligner indexes
if [ -f $genome ]
then
    if [ ! -f $genome_dir/STAR ]
    then
    fi
    
    if [ ! -f $genome_dir/bwa-meth ]
    then
        # Build the bwa-meth index
        mkdir $genome_dir/bwa-meth
        ln -s $genome $genome_dir/bwa-meth/ref.fa
        bwameth.py index $genome_dir/bwa-meth/ref.fa
    fi
    
    if [ ! -f $genome_dir/STAR ]
    then
        # Build the STAR index
        # Set --runThreadN to the number of threads available on your machine
        STAR --runMode genomeGenerate --genomeDir $genome_dir/STAR --genomeFastaFiles \
        $genome --runThreadN 24 --sjdbGTFfile $annotations --sjdbOverhang 75
    fi
fi
