#################################################################
# Dockerfile
#
# Software:         bwa and bwa-meth, plus genome indexes
# Software Version: 0.7.15|0.10
# Description:      Alpine bwa and bwa-meth Image with genome indexes
# Website:          https://github.com/lh3/bwa
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa|bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwaindex:latest .
# Pull Cmd:         docker pull jdidion/bwaindex
# Run Cmd:          docker run --rm jdidion/bwaindex bwa
#                   docker run --rm jdidion/bwaindex python bwa-meth.py
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
WORKDIR /tmp

RUN apk --update add wget

# Build bwa indexes
RUN bwa index /genomes/hg37.fa
# Build bwa-meth indexes
RUN bwameth.py index /genomes/hg37.fa

RUN bwa index /genomes/hg38.fa
RUN bwameth.py index /genomes/hg38.fa

RUN rm -rf /var/cache/apk/* /tmp/* && \
    apk del deps
