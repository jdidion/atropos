#################################################################
# Dockerfile
#
# Software:         bwa plus genome indexes
# Software Version: 0.7.15
# Description:      Alpine bwa image with genome indexes
# Website:          https://github.com/lh3/bwa
# Provides:         bwa|/data/index/bwa/hg37
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwaindex:latest .
# Pull Cmd:         docker pull jdidion/bwaindex
# Run Cmd:          docker run --rm jdidion/bwaindex bwa
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
ADD bwa-index/* /data/index/bwa/hg37/
