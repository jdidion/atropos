#################################################################
# Dockerfile
#
# Software:         bwa-meth plus genome indexes
# Software Version: 0.10
# Description:      Alpine bwa-meth Image with genome indexes
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwamethindex:latest .
# Pull Cmd:         docker pull jdidion/bwamethindex
# Run Cmd:          docker run --rm jdidion/bwaindex python bwa-meth.py
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD bwameth-index/* /data/index/bwameth/hg37
