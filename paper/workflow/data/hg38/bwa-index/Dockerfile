#################################################################
# Dockerfile
#
# Software:         bwa plus hg38 genome indexes
# Software Version: 0.7.15
# Description:      Alpine bwa image with genome indexes
# Website:          https://github.com/lh3/bwa
# Provides:         bwa
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwaindex:latest .
# Pull Cmd:         docker pull jdidion/bwaindex
# Run Cmd:          docker run --rm jdidion/bwaindex bwa
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD hg38* /data/index/bwa/hg38/
