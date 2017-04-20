#################################################################
# Dockerfile
#
# Software:         bedops, SRA tools
# Software Version: 2.4.26, 2.8.2
# Description:      Alpine image with tools to support the Atropos paper workflow
# Website:          https://bedops.readthedocs.io/en/latest/
#                   https://github.com/ncbi/sra-tools
# Provides:         bedops|fastq-dump
# Base Image:       alpine
# Build Cmd:        docker build -t jdidion/atropos_paper:latest .
# Pull Cmd:         docker pull jdidion/atropos_paper
# Run Cmd:          docker run --rm jdidion/atropos_paper <cmd>
#################################################################
FROM python:alpine
WORKDIR /tmp

RUN apk add --no-cache g++ make git wget

ENV GENCODE_VERSION '26'

# Fetch GENCODE GTFs
RUN wget -o /annotations 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh37_mapping/gencode.v${GENCODE_VERSION}lift37.annotation.gtf.gz'
RUN wget -o /annotations 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz'

# bedops
ENV VERSION 2.4.26
ENV URL "https://github.com/bedops/bedops/archive/v${VERISON}.tar.gz"
RUN wget -q -O - $URL | tar -zxv && \
    cd v${VERSION} && \
    make && \
    make install && \
    cp bin/* /usr/local/bin

# ngs toolkit
ENV VERSION 2.8.2
ENV URL "https://github.com/ncbi/ngs-tools/archive/${VERSION}.tar.gz"
RUN wget -q -O - $URL | tar -zxv && \
    cd ${VERSION} &&
    ./configure --build-prefix=/usr/local/ &&
    make &&
    make install

# add my scripts
ADD ../../scripts /scripts

RUN rm -rf /var/cache/apk/* /tmp/* && \
    apk del deps