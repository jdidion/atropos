#################################################################
# Dockerfile
#
# Software:         STAR, plus genome indexes
# Software Version: 2.5.3a
# Description:      Phusion STAR RNA-Seq aligner Image with genome indexes
# Website:          https://github.com/alexdobin/STAR
# Provides:         STAR
# Base Image:       jedisct1/phusion-baseimage-latest
# Build Cmd:        docker build -t jdidion/star_hg38index:latest .
# Pull Cmd:         docker pull jdidion/star_hg38index
# Run Cmd:          docker run --rm jdidion/star_hg38index STAR
#################################################################
FROM jdidion/bwa
ADD index/* /data/index/star/hg38/
