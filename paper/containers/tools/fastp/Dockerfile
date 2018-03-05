#################################################################
# Dockerfile
#
# Software:         fastp
# Software Version: 0.12.3
# Description:      fastp image
# Website:          https://github.com/OpenGene/fastp
# Provides:         fastp
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -f Dockerfile -t jdidion/fastp:latest .
#                   docker tag jdidion/fastp:latest jdidion/fastp
# Pull Cmd:         docker pull jdidion/fastp
# Run Cmd:          docker run jdidion/fastp fastp --help
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
ENV VERSION v0.12.3
ENV BUILD_PACKAGES \
    build-essential \
    g++ \
    libbz2-dev \
    zlib1g-dev \
    git \ 
    perl
ENV RUNTIME_PACKAGES \
    time
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && git clone --recursive --branch $VERSION https://github.com/OpenGene/fastp \
    && cd fastp \
    && make \
    && mv /tmp/fastp/fastp /usr/local/bin \
    && cd .. \
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8
