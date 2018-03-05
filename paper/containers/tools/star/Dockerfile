#################################################################
# Dockerfile
#
# Software:         STAR
# Software Version: 2.5.4b
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

ENV BUILD_PACKAGES \
    build-essential \
    gcc-multilib \
    apt-utils \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    wget
ENV RUNTIME_PACKAGES \
    time

# STAR
ENV STAR_VERSION 2.5.4b
ENV STAR_URL "https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz"

# Install samtools
ENV SAMTOOLS_VERSION 1.7
ENV SAMTOOLS_URL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

# Install compiler and perl stuff
RUN apt-get update && apt-get install --yes \
    $BUILD_PACKAGES \
    $RUNTIME_PACKAGES \
&& wget -q -O - $STAR_URL | tar -zxv \
&& cd STAR-${STAR_VERSION}/source \
&& make STAR LDFLAGSextra=-flto CXXFLAGSextra="-flto -march=native" \
&& cd .. \
&& cp ./bin/Linux_x86_64_static/STAR /usr/local/bin/ \
&& cd .. \
&& wget -q -O - $SAMTOOLS_URL | tar -jxv \
&& cd samtools-${SAMTOOLS_VERSION} \
&& ./configure --prefix=/usr/local/ --without-curses \
&& make \
&& make install \
&& rm -Rf /tmp/samtools-${SAMTOOLS_VERSION} \
&& apt-get remove --purge -y $BUILD_PACKAGES \
&& rm -rf /var/lib/apt/lists/* \
&& rm -rf /tmp/*
