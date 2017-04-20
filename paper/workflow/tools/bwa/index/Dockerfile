#################################################################
# Dockerfile
#
# Software:         bwa and bwa-meth, plus genome indexes
# Software Version: 0.7.15, 0.10
# Description:      Alpine bwa and bwa-meth Image with genome indexes
# Website:          https://github.com/lh3/bwa
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa, bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwaindex:latest .
# Pull Cmd:         docker pull jdidion/bwaindex
# Run Cmd:          docker run --rm jdidion/bwaindex bwa
#                   docker run --rm jdidion/bwaindex python bwa-meth.py
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
WORKDIR /tmp
RUN mkdir /genomes

# Fetch hg19/37
RUN scp 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/2*' . &&
    zcat *.gz > /genomes/hg37.fa &&
    samtools faidx /genomes/hg37.fa

# Build bwa index
RUN bwa index /genomes/hg37.fa

# Fetch hg38
#RUN scp 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' . &&
#    gunzip hg38.fa.gz &&
#    mv hg38.fa.gz /genomes &&
#    samtools faidx /genomes/hg38.fa
