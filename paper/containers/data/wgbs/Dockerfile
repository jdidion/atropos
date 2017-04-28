#################################################################
# Data:             WGBS data from SRA ENCODE accession ENCSR890UQO
# Website:          https://www.encodeproject.org/experiments/ENCSR890UQO/
# Provides:         /data/wgbs/wgbs.1.fq.gz|/data/rna/wgbs.2.fq.gz
# Adapters:         AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG # TruSeq index 7
#                   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT # TruSeq universal
# Base Image:       busybox
#################################################################
FROM blang/busybox-bash
RUN mkdir -p /data/wgbs
ADD wgbs.1.fq.gz /data/wgbs/wgbs.1.fq.gz
ADD wgbs.2.fq.gz /data/wgbs/wgbs.2.fq.gz
