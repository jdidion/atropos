#################################################################
# Dockerfile
#
# Software:         bwa-meth plus hg37 genome indexes
# Software Version: 0.10
# Description:      Alpine bwa-meth Image with genome indexes
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwameth_hg37index:latest .
# Pull Cmd:         docker pull jdidion/bwameth_hg37index
# Run Cmd:          docker run --rm jdidion/bwameth_hg37index python bwa-meth.py
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD hg37.fa* /data/index/bwameth/hg37/
