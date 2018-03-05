#################################################################
# Dockerfile
#
# Software:         Cutadapt
# Software Version: 1.16
# Description:      Cutadapt trimming software
# Website:          https://github.com/marcelm/cutadapt
# Provides:         cutadapt
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -f Dockerfile -t jdidion/cutadapt:latest .
# Pull Cmd:         docker pull jdidion/cutadapt
# Run Cmd:          docker run --rm jdidion/cutadapt cutadapt -h
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
ENV VERSION v1.16
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
    && git clone --recursive --branch $VERSION https://github.com/marcelm/cutadapt \
    && cd cutadapt \
    && python setup.py install \
    && python setup.py build_ext -i \
    && cd ..\
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8
