#################################################################
# Data:             Human genome build 37 (from UCSC)
# Website:          ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes
# Provides:         /data/reference/hg37/hg37.fa|/data/reference/hg37/hg37.fa.fai
# Base Image:       busybox
#################################################################
FROM blang/busybox-bash
RUN mkdir -p /data/reference/hg37 \
 && mkdir -p /data/annotations/hg37
ADD hg37.fa /data/reference/hg37/hg37.fa
ADD hg37.fa.fai /data/reference/hg37/hg37.fa.fai
ADD gencode.v26lift37.annotation.gtf /data/annotations/hg37/gencode.v26lift37.annotation.gtf
