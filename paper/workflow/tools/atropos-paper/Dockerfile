#################################################################
# Dockerfile
#
# Software:         bedops, SRA tools
# Software Version: 2.4.26, 2.8.2
# Description:      Alpine image with tools to support the Atropos paper workflow
# Website:          https://bedops.readthedocs.io/en/latest/
#                   https://github.com/ncbi/sra-tools
# Provides:         bedops|/opt/fastq-dump/fastq-dump-wrapper.sh
# Base Image:       alpine
# Build Cmd:        docker build -t jdidion/atropos_paper:latest .
# Pull Cmd:         docker pull jdidion/atropos_paper
# Run Cmd:          docker run --rm jdidion/atropos_paper <cmd>
#################################################################
FROM jedisct1/phusion-baseimage-latest
WORKDIR /tmp
RUN mkdir /annotations

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y libc6-dev build-essential tcsh devscripts debhelper git wget
RUN apt-get install ca-certificates openssl && update-ca-certificates
RUN apt-get clean

ENV GENCODE_VERSION '26'

# Fetch GENCODE GTFs
ENV FILE "gencode.v${GENCODE_VERSION}lift37.annotation.gtf.gz"
RUN wget -o /annotations/${FILE} "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh37_mapping/${FILE}"
ENV FILE "gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
RUN wget -o /annotations/${FILE} "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/${FILE}"

# bedops
ENV VERSION 2.4.26
ENV URL "https://github.com/bedops/bedops/archive/v${VERSION}.tar.gz"
RUN wget -q -O - $URL | tar -zxv && \
    cd bedops-${VERSION} && \
    make && \
    make install && \
    cp bin/* /usr/local/bin

# sra-tools
ENV VERSION 2.8.2-1
ENV URL "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz"
RUN wget 
RUN tar -xvf sratoolkit.2.8.2-1-ubuntu64.tar.gz
RUN rm sratoolkit.2.8.2-1-ubuntu64.tar.gz
RUN mkdir /opt/fastq-dump/
RUN sudo mv sratoolkit.2.8.2-1-ubuntu64 /opt/fastq-dump/

ADD fastq-dump-wrapper.sh /opt/fastq-dump/
RUN mkdir /data
WORKDIR /data

RUN apt-get remove --purge -y $BUILD_PACKAGES $(apt-mark showauto) && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/*
