#################################################################
# Dockerfile
#
# Software:         SeqPurge
# Software Version: 0.1 (9e8d99d)
# Description:      SeqPurge image
# Website:          https://github.com/imgag/ngs-bits
# Provides:         SeqPurge
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -f Dockerfile -t imgag/ngs-bits:latest .
#                   docker tag imgag/ngs-bits:latest jdidion/seqpurge
# Pull Cmd:         docker pull jdidion/seqpurge
# Run Cmd:          docker run
# Note: Apparently they do not use tags or releases for ngs-bits, so I
# specifically checkout the commit that contains the most recent update
# to SeqPurge (as of this writing).
# Note: For me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
ENV VERSION 9e8d99d
ENV BUILD_PACKAGES \
    software-properties-common \
    git \
    cmake \
    build-essential \
    g++ \
    libbz2-dev \
    zlib1g-dev \
    qt5-qmake \
    libqt5sql5-mysql
ENV RUNTIME_PACKAGES \
    qt5-default \
    libqt5xmlpatterns5-dev \
    time
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
        $RUNTIME_PACKAGES \
    && git clone --recursive https://github.com/imgag/ngs-bits.git \
    && git checkout $VERSION \
    && cd ngs-bits \
    && make build_3rdparty \
    && make build_tools_release \
    && mv /tmp/ngs-bits/bin/* /usr/local/bin \
    && cd .. \
    && rm -rf /var/lib/apt/lists/* /tmp/* \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8
