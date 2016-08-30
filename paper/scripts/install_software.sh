#!/bin/bash
# This script installs all the software in the ../software
# directory into the ../bin directory. It also installs
# Atropos from pypi, and it installs bwameth from GitHub.
# Prerequisites:
# 1. modern C++ compiler (we use homebrew/linuxbrew to install gcc 5.1)
# 2. python 3.3+ (we recommend Anaconda)
#    * cython 0.24+ ('conda install cython' or 'pip install cython')
# 3. To install SeqPurge, you need
#    * Qt 5.3+ with xmlpatterns and mysql packages
#    * git
#    * cmake
# 4. To run the analyses on real data, you need
#    * BWA (the MEM command is required; we used version 0.7.15)
#    * samtools

scripts=`pwd`
root=`dirname $scripts`
mkdir ../software/build
automake_dir=/usr/local/Cellar/automake/1.15/share/automake-1.15
# Set this to the location of the hg19 reference multifasta.
# If this file doesn't exist, the bwameth index won't be built.
genome=../data/ref.fa

# We need some python libraries for the benchmarking scripts
pip install editdistance

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
pip install atropos==1.0.12

# Install Skewer
mkdir ../software/build/skewer &&
    cd ../software/build/skewer &&
    cp ../../skewer_2016.08.04.zip . &&
    unzip skewer_2016.08.04.zip &&
    cd skewer-master &&
    sed -i -e 's/\/usr\/local\/bin/..\/..\/..\/bin/' Makefile &&
    make &&
    make install &&
    cd ../../../scripts

# Install SeqPurge
mkdir ../software/build/seqpurge &&
    cd ../software/build/seqpurge &&
    cp ../../ngs-bits_2016.08.04.zip . &&
    unzip ngs-bits_2016.08.04.zip &&
    cd ngs-bits &&
    make build_3rdparty &&
    make build_tools_release &&
    ln -s bin/SeqPurge ../../../bin &&
    cd ../../../scripts &&
    echo "You may need to 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$root/software/bin'"

if [ -f $genome ]
then
    # Install bwameth
    pip install toolshed
    mkdir ../software/build/bwameth &&
        cd ../software/build/bwameth &&
        wget https://github.com/brentp/bwa-meth/archive/v0.10.tar.gz &&
        tar xzvf v0.10.tar.gz &&
        cd bwa-meth-0.10/ &&
        sudo python setup.py install
    # Build the index
    bwameth.py index $genome
fi
