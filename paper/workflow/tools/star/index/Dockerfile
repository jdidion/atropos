#################################################################
# Dockerfile
#
# Software:         STAR
# Software Version: 2.5.3a
# Description:      Alpine STAR RNA-Seq aligner Image with genome index
# Website:          https://github.com/alexdobin/STAR
# Provides:         STAR
# Base Image:       alpine
# Build Cmd:        docker build -t jdidion/starindex:latest .
# Pull Cmd:         docker pull jdidion/starindex
# Run Cmd:          docker run jdidion/starindex
#################################################################
# Using jdidion/bwabase as the base image
FROM jdidion/bwabase
WORKDIR /tmp
RUN mkdir /genomes
RUN mkdir /annotations

# Fetch hg19/37
RUN scp 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/2*' . &&
    zcat *.gz > /genomes/hg37.fa &&
    samtools faidx /genomes/hg37.fa

# Fetch GENCODE GTF
RUN wget -o /annotations 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'

# Build star index
RUN STAR --runMode genomeGenerate --genomeDir /genomes --genomeFastaFiles /genomes/hg37.fa --sjdbGTFfile /annotations/gencode.v19.annotation.gtf.gz

# Fetch hg38
#RUN scp 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' . &&
#    gunzip hg38.fa.gz &&
#    mv hg38.fa.gz /genomes &&
#    samtools faidx /genomes/hg38.fa

