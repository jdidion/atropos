#!/bin/bash
# download first 1M reads
wget -qO- https://www.encodeproject.org/files/ENCFF798RSS/@@download/ENCFF798RSS.fastq.gz | gunzip | head -4000000 | gzip > wgbs.1.fq.gz \
&& wget -qO- https://www.encodeproject.org/files/ENCFF113KRQ/@@download/ENCFF113KRQ.fastq.gz | gunzip | head -4000000 | gzip > wgbs.2.fq.gz \
&& docker build -f Dockerfile -t jdidion/atropos_wgbs:latest . \
&& rm wgbs.*.fq.gz

