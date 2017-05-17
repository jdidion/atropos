#################################################################
# Dockerfile
#
# Software:         bwa plus hg38 genome indexes
# Software Version: 0.7.15
# Description:      Alpine bwa image with genome indexes
# Website:          https://github.com/lh3/bwa
# Provides:         bwa
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwa_hg38index:latest .
# Pull Cmd:         docker pull jdidion/bwa_hg38index
# Run Cmd:          docker run --rm jdidion/bwa_hg38index bwa
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD index/* /data/index/bwa/hg38/
