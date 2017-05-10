#################################################################
# Dockerfile
#
# Software:         bedops, SRA tools
# Software Version: 2.4.26, 2.8.2
# Description:      Phusion image with tools to support the Atropos paper workflow
# Website:          https://bedops.readthedocs.io/en/latest/
#                   https://github.com/ncbi/sra-tools
# Provides:         bedops|/opt/fastq-dump/fastq-dump-wrapper.sh
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -t jdidion/atropos_paper_analysis:latest .
# Pull Cmd:         docker pull jdidion/atropos_paper_analysis
# Run Cmd:          docker run --rm jdidion/atropos_paper_analysis <cmd>
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
RUN mkdir /annotations
RUN mkdir /data

ENV GENCODE_VERSION '26'
ENV GENCODE_37_FILE "gencode.v${GENCODE_VERSION}lift37.annotation.gtf.gz"
ENV GENCODE_38_FILE "gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
ENV BEDOPS_VERSION 2.4.26
ENV BEDOPS_URL "https://github.com/bedops/bedops/archive/v${BEDOPS_VERSION}.tar.gz"
ENV SRA_VERSION 2.8.2-1
ENV SRA_URL "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz"
ENV BUILD_PACKAGES \
    libc6-dev \
    build-essential \
    tcsh \
    devscripts \
    debhelper \
    git \
    wget \
    ca-certificates \
    openssl
ENV RUNTIME_PACKAGES \
    time

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && update-ca-certificates \
    && wget -o /annotations/${GENCODE_37_FILE} "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh37_mapping/${GENCODE_37_FILE}" \
    && wget -o /annotations/${GENCODE_38_FILE} "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/${GENCODE_38_FILE}" \
    && wget -q -O - $BEDOPS_URL | tar -zxv \
    && cd bedops-${BEDOPS_VERSION} \
    && make \
    && make install \
    && cp bin/* /usr/local/bin \
    && cd .. \
    && wget -q -O - $SRA_URL | tar -zxv \
    && mkdir /opt/fastq-dump/ \
    && mv sratoolkit.${SRA_VERSION}-ubuntu64 /opt/fastq-dump/ \
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8

ADD fastq-dump-wrapper.sh /opt/fastq-dump/
WORKDIR /data
