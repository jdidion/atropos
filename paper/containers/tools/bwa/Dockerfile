#################################################################
# Dockerfile
#
# Software:         bwa and bwa-meth
# Software Version: 0.7.15, 0.10
# Description:      Alpine bwa and bwa-meth Image
# Website:          https://github.com/lh3/bwa
#                   https://github.com/brentp/bwa-meth
# Provides:         bwa|bwa-meth.py
# Base Image:       alpine:python
# Build Cmd:        docker build -t jdidion/bwa:latest .
# Pull Cmd:         docker pull jdidion/bwa
# Run Cmd:          docker run --rm jdidion/bwa bwa
#                   docker run --rm jdidion/bwa python bwa-meth.py
#################################################################
# python3 Alpine image with gcc
FROM python:alpine
WORKDIR /tmp

RUN apk add --update --no-cache ncurses
RUN apk add --virtual=deps --update --no-cache musl-dev zlib-dev bzip2-dev xz-dev make gcc git
RUN apk add --update ca-certificates openssl && update-ca-certificates

### Install bwa
ENV VERSION 0.7.15
ENV URL "https://github.com/lh3/bwa/releases/download/v${VERSION}/bwa-${VERSION}.tar.bz2"
RUN wget -q -O - $URL | tar -jxv && \
    cd bwa-${VERSION} && \
    sed -i '1i#include <stdint.h>' kthread.c && \
    sed -i[.bak] "s/u_int32_t/uint32_t/g" *.c && \
    sed -i[.bak] "s/u_int32_t/uint32_t/g" *.h && \
    make && \
    mv /tmp/bwa-${VERSION}/bwa /usr/local/bin

# Install bwa-meth
RUN git clone --recursive https://github.com/brentp/bwa-meth && \
    cd bwa-meth && \
    python setup.py install

# Install samtools
ENV VERSION 1.4
ENV URL "https://github.com/samtools/samtools/releases/download/1.4/samtools-${VERSION}.tar.bz2"
RUN wget -q -O - $URL | tar -jxv && \
    cd samtools-${VERSION} && \
    ./configure --prefix=/usr/local/ --without-curses && \
    make && \
    make install && \
    rm -Rf /tmp/samtools-${VERSION}

RUN rm -rf /var/cache/apk/* /tmp/* && \
    apk del deps
