#!/bin/bash
scripts=`pwd`
root=`dirname $scripts`
mkdir $root/software/build
automake_dir=/usr/local/Cellar/automake/1.15/share/automake-1.15
# Set this to the location of the hg19 reference multifasta.
# If this file doesn't exist, the bwameth index won't be built.
genome_dir=$root/data/reference
genome=$genome_dir/ref.fa
annotations=$genome_dir/gencode.v19.annotation.gtf

# For the benchmarking script, there is the option to compute
# edit distance between the untrimmed and trimmed reads. We
# did not use those metrics in the paper. If you want to enable
# the edit distance calculation, you need to install the 'editdistance'
# pyton library.
# pip install editdistance
# Note: if this doesn't work for you, you'll need to checkout
# the editdistance repository and edit setup.py to enable the
# 'cythonize' command, which will recompile the cython code for
# your local environment. So you would run:
# python setup.py build_ext -i && python setup.py install

# Install modified ART
mkdir ../software/build/art &&
    cd ../software/build/art &&
    cp ../../art_illumina_src151.tar.gz . &&
    tar -xzf art_illumina_src151.tar.gz &&
    cp ../../art_illumina_src151-adapter-enabled.tar.gz . &&
    tar -xzf art_illumina_src151-adapter-enabled.tar.gz &&
    cd art_illumina_dir &&
    for f in config.sub config.guess install-sh depcomp missing INSTALL
    do
    rm $f
    ln -s $automake_dir/$f .
    done &&
    ./configure --prefix $root/software &&
    make &&
    make install &&
    cd ../../../scripts

# Install version of Atropos on which manuscript is based
pip install atropos==1.0.23

# Install Skewer from conda
conda install skewer

# Install SeqPurge
# We tested using the version from GitHub represented by commit
# 8713481a9a7404cb3e69f7660b94d9847dbe632b
cd ../software/build &&
    git clone --recursive https://github.com/imgag/ngs-bits.git &&
    cd ngs-bits &&
    git checkout 8713481a9a7404cb3e69f7660b94d9847dbe632b &&
    make build_3rdparty &&
    make build_tools_release &&
    ln -s bin/SeqPurge ../../../bin &&
    cd ../../../scripts &&
    echo "You may need to 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$root/software/bin'"

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
