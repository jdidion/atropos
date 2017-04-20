#################################################################
# Dockerfile
#
# Software:         bwa and bwa-meth
# Software Version: 0.7.15, 0.10
# Description:      Alpine bwa and bwa-meth Image
# Website:          https://github.com/lh3/bwa
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa, bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwabase:latest .
# Pull Cmd:         docker pull jdidion/bwabase
# Run Cmd:          docker run --rm jdidion/bwabase bwa
#                   docker run --rm jdidion/bwabase python bwa-meth.py
#################################################################
# python3 Alpine image with gcc
FROM python:alpine
ENV DEBIAN_FRONTEND noninteractive
WORKDIR /tmp

### Install bwa
ENV VERSION 0.7.15
ENV NAME bwa
ENV URL "https://github.com/lh3/bwa/archive/${VERSION}.tar.gz"
RUN wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION} && \
    make -j 4 && \
    cd .. && \
    cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
    strip /usr/local/bin/${NAME}; true && \
    rm -rf ./${NAME}-${VERSION}/

# Install bwa-meth
RUN git clone --recursive https://github.com/brentp/bwa-meth && \
    cd bwa-meth && \
    python setup.py install

# Install samtools
ENV VERSION 1.4
ENV URL "https://github.com/samtools/samtools/releases/download/1.4/samtools-${VERSION}.tar.bz2"
RUN wget -q -O - $URL | tar -zxv && \
    cd samtools-${VERSION} &&
    make &&
    make --prefix=/usr/local/bin/ install
