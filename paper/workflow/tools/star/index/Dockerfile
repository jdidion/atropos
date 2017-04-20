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
FROM jdidion/starbase
WORKDIR /tmp
RUN mkdir /genomes
RUN mkdir /annotations
ENV GENCODE_VERSION 26
# Number of threads to use for building index
ENV THREADS 2

RUN apt-get update && apt-get install --yes wget gzip

# hg19/37

## download reference genome
RUN wget -r --no-parent --no-directories -A 'chr*' 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes' && \
    zcat *.gz > /genomes/hg37.fa && \
    samtools faidx /genomes/hg37.fa

## download GTF
ENV FILE "gencode.v${GENCODE_VERSION}lift37.annotation.gtf.gz"
ENV URL "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh37_mapping/${FILE}"
RUN wget --no-parent --no-directories -O /annotations/${FILE} ${URL} && \
    gunzip /annotations/${FILE}

## build index
RUN STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /genomes --genomeFastaFiles /genomes/hg37.fa \
        --sjdbGTFfile /annotations/gencode.v${GENCODE_VERSION}lift37.annotation.gtf \
        --sjdbOverhang 75

# hg38

## download reference genome
RUN wget -O hg38.fa.gz 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' && \
    gunzip hg38.fa.gz && \
    mv hg38.fa /genomes && \
    samtools faidx /genomes/hg38.fa

# Fetch GENCODE GTF
ENV FILE "gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
ENV URL "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/${FILE}"
RUN wget --no-parent --no-directories -O /annotations/${FILE} ${URL} && \
    gunzip /annotations/${FILE}

# Build star index
RUN STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /genomes --genomeFastaFiles /genomes/hg38.fa \
        --sjdbGTFfile /annotations/gencode.v${GENCODE_VERSION}.annotation.gtf \
        --sjdbOverhang 75

RUN apt-get remove --purge -y $BUILD_PACKAGES $(apt-mark showauto) && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/*
