#################################################################
# Dockerfile
#
# Software:         python
# Software Version: 3.6
# Description:      Python 3.6 pandas, seaborn, and other libraries
# Provides:         python
# Base Image:       phusion/baseimage
# Build Cmd:        docker build -f Dockerfile -t jdidion/python_bash .
# Pull Cmd:         docker pull jdidion/python_bash
# Run Cmd:          docker run --rm jdidion/python_bash python -h
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp

ENV BUILD_PACKAGES \
    build-essential \
    g++ \
    gfortran \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libpng-dev \
    libfreetype6-dev \
    git
ENV RUNTIME_PACKAGES \
    libncurses5-dev \
    libncursesw5-dev \
    python3-pip \
    time

RUN apt-get update && apt-get install --yes \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && ln -s /usr/bin/pip3 /usr/bin/pip \
    && pip install --upgrade pip \
    && pip install cython \
    && pip install editdistance \
    && pip install numpy \
    && pip install pandas \
    && pip install matplotlib \
    && pip install seaborn \
    && pip install mako \
    && pip install pysam \
    && pip uninstall --yes cython \
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*
