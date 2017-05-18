#################################################################
# Dockerfile
#
# Software:         STAR
# Software Version: 2.5.3a
# Description:      Phusion STAR RNA-Seq aligner Image
# Website:          https://github.com/alexdobin/STAR
# Provides:         STAR
# Base Image:       jedisct1/phusion-baseimage-latest
# Build Cmd:        docker build -t jdidion/starbase:latest .
# Pull Cmd:         docker pull jdidion/starbase
# Run Cmd:          docker run --rm jdidion/starbase STAR
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp

# Install compiler and perl stuff
RUN apt-get update && apt-get install --yes \
    build-essential \
    gcc-multilib \
    apt-utils \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    wget

# STAR
ENV VERSION 2.5.3a
ENV URL "https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz"
RUN wget -q -O - $URL | tar -zxv && \
    cd STAR-${VERSION}/source && \
    make STAR LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=native" && \
    cd .. && \
    cp ./bin/Linux_x86_64_static/STAR /usr/local/bin/

# Install samtools
ENV VERSION 1.4
ENV URL "https://github.com/samtools/samtools/releases/download/1.4/samtools-${VERSION}.tar.bz2"
RUN wget -q -O - $URL | tar -jxv && \
    cd samtools-${VERSION} && \
    ./configure --prefix=/usr/local/ --without-curses && \
    make && \
    make install && \
    rm -Rf /tmp/samtools-${VERSION}

RUN apt-get remove --purge -y $BUILD_PACKAGES $(apt-mark showauto) && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/*
