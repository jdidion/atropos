#################################################################
# Data:             Human genome build 38 (from UCSC)
# Website:          http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# Provides:         /data/hg38/reference/hg38.fa|/data/hg38/reference/hg38.fa.fai
# Base Image:       busybox
#################################################################
FROM busybox
RUN mkdir -p /data/reference/hg38
ADD hg38.fa /data/reference/hg38/hg38.fa
ADD hg38.fa.fai /data/reference/hg38/hg38.fa.fai