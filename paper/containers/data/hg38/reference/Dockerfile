#################################################################
# Data:             Human genome build 38 (from UCSC)
# Website:          http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# Provides:         /data/reference/ng38/hg38.fa|/data/reference/hg38/hg38.fa.fai
# Base Image:       busybox
#################################################################
FROM blang/busybox-bash
RUN mkdir -p /data/reference/hg38 \
 && mkdir -p /data/annotations/hg38
ADD hg38.fa /data/reference/hg38/hg38.fa
ADD hg38.fa.fai /data/reference/hg38/hg38.fa.fai
ADD gencode.v26.annotation.gtf /data/annotations/hg38/gencode.v26.annotation.gtf
