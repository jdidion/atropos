#################################################################
# Dockerfile
#
# Software:         Atropos
# Software Version: 1.1.4
# Description:      Atropos image for paper
# Website:          https://github.com/jdidion/atropos
# Provides:         atropos
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -f Dockerfile -t jdidion/atropos_paper:latest .
# Pull Cmd:         docker pull jdidion/atropos_paper
# Run Cmd:          docker run --rm jdidion/atropos_paper atropos -h
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
ENV VERSION 1.1.5
ENV BUILD_PACKAGES \
    git
ENV RUNTIME_PACKAGES \
    python3-pip \
    time
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && ln -s /usr/bin/pip3 /usr/bin/pip \
    && pip install cython \
    && git clone --recursive --branch $VERSION https://github.com/jdidion/atropos \
    && cd atropos \
    && make install \
    && cd ..\
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8
