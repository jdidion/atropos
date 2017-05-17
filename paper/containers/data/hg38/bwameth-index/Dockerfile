#################################################################
# Dockerfile
#
# Software:         bwa-meth plus hg38 genome indexes
# Software Version: 0.10
# Description:      Alpine bwa-meth image with genome indexes
# Website:          https://github.com/brentp/bwa-meth
# Provides:         bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwameth_hg38index:latest .
# Pull Cmd:         docker pull jdidion/bwameth_hg38index
# Run Cmd:          docker run --rm jdidion/bwameth_hg38index python bwa-meth.py
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD index/* /data/index/bwameth/hg38/
