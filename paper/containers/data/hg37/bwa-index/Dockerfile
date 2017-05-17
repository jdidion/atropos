#################################################################
# Dockerfile
#
# Software:         bwa plus hg37 genome indexes
# Software Version: 0.7.15
# Description:      Alpine bwa image with genome indexes
# Website:          https://github.com/lh3/bwa
# Provides:         bwa|/data/index/bwa/hg37
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwa_hg37index:latest .
# Pull Cmd:         docker pull jdidion/bwa_hg37index
# Run Cmd:          docker run --rm jdidion/bwa_hg37index bwa
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD hg37.fa* /data/index/bwa/hg37/
