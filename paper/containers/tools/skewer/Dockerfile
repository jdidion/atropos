#################################################################
# Dockerfile
#
# Software:         skewer
# Software Version: 0.2.2
# Description:      skewer image
# Website:          https://github.com/relipmoc/skewer
# Provides:         skewer
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -f Dockerfile -t relipmoc/skewer:latest .
#                   docker tag relipmoc/skewer:latest jdidion/seqpurge
# Pull Cmd:         docker pull jdidion/skewer
# Run Cmd:          docker run jdidion/skewer
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
ENV VERSION 0.2.2
ENV BUILD_PACKAGES \
    build-essential \
    g++ \
    git
ENV RUNTIME_PACKAGES \
    time    
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && git clone --recursive --branch $VERSION https://github.com/relipmoc/skewer \
    && cd skewer \
    && make \
    && mv /tmp/skewer/skewer /usr/local/bin \
    && cd .. \
    && apt-get remove --purge -y $BUILD_PACKAGES \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8
