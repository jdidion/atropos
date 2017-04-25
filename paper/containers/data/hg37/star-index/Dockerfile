#################################################################
# Dockerfile
#
# Software:         STAR, plus genome indexes
# Software Version: 2.5.3a
# Description:      Phusion STAR RNA-Seq aligner Image with genome indexes
# Website:          https://github.com/alexdobin/STAR
# Provides:         STAR
# Base Image:       jedisct1/phusion-baseimage-latest
# Build Cmd:        docker build -t jdidion/starbase:latest .
# Pull Cmd:         docker pull jdidion/starbase
# Run Cmd:          docker run --rm jdidion/starbase STAR
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/starbase
ADD index/* /data/index/star/hg37/
