#################################################################
# Data:             RNA-Seq data from SRA accession SRR521458
# Website:          https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR521458
# Provides:         /data/rna/rna.1.fq.gz|/data/rna/rna.2.fq.gz
# Adapters:         AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATATCGTATGCCGTCTTCTGCTTG # Custom?
#                   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT # TruSeq Universal
# Base Image:       busybox
#################################################################
FROM blang/busybox-bash
RUN mkdir -p /data/rna
ADD SRR521458_1.fastq.gz /data/rna/rna.1.fq.gz
ADD SRR521458_2.fastq.gz /data/rna/rna.2.fq.gz
