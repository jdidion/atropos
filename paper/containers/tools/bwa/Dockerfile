#################################################################
# Dockerfile
#
# Software:         bwa, bwa-meth, samtools
# Software Version: 0.7.17, 0.20, 1.7
# Description:      Alpine bwa and bwa-meth Image
# Website:          https://github.com/lh3/bwa
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa|bwa-meth.py|samtools
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwabase:latest .
# Pull Cmd:         docker pull jdidion/bwabase
# Run Cmd:          docker run --rm jdidion/bwabase bwa
#                   docker run --rm jdidion/bwabase python bwa-meth.py
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp

ENV BUILD_PACKAGES \
    build-essential \
    g++ \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    git \
    wget
ENV RUNTIME_PACKAGES \
    libncurses5-dev \
    libncursesw5-dev \
    python3-pip \
    time
ENV BWA_VERSION 0.7.17
ENV BWA_URL "https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2"
ENV SAMTOOLS_VERSION 1.7
ENV SAMTOOLS_URL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

RUN apt-get update && apt-get install --yes \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && ln -s /usr/bin/pip3 /usr/bin/pip \
    && wget -q -O - $BWA_URL | tar -jxv \
    && cd bwa-${BWA_VERSION} \
    && sed -i '1i#include <stdint.h>' kthread.c \
    && sed -i[.bak] "s/u_int32_t/uint32_t/g" *.c \
    && sed -i[.bak] "s/u_int32_t/uint32_t/g" *.h \
    && make \
    && mv /tmp/bwa-${BWA_VERSION}/bwa /usr/local/bin \
    && cd .. \
    && git clone --recursive https://github.com/brentp/bwa-meth \
    && cd bwa-meth \
    && python setup.py install \
    && wget -q -O - $SAMTOOLS_URL | tar -jxv \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local/ --without-curses \
    && make \
    && make install \
    && rm -Rf /tmp/samtools-${SAMTOOLS_VERSION}

RUN apt-get remove --purge -y $BUILD_PACKAGES \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*
