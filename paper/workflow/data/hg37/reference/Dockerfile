#################################################################
# Data:             Human genome build 37 (from UCSC)
# Website:          ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes
# Provides:         /data/hg37/reference/hg37.fa|/data/hg37/reference/hg37.fa.fai
# Base Image:       busybox
#################################################################
FROM busybox
RUN mkdir -p /data/reference/hg37
ADD hg37.fa /data/reference/hg37/hg37.fa
ADD hg37.fa.fai /data/reference/hg37/hg37.fa.fai