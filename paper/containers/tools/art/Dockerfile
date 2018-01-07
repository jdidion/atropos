#################################################################
# Dockerfile
#
# Software:         ART (modified)
# Software Version: 1.5.1
# Description:      ART read simulator, modified to introduce adapter sequences
# Website:          https://sourceforge.net/projects/skewer/files/Simulator/
# Provides:         ART
# Base Image:       phusion/baseimage:latest
# Build Cmd:        docker build -f Dockerfile -t jdidion/art_skewer:latest .
# Pull Cmd:         docker pull jdidion/art_skewer
# Run Cmd:          docker run jdidion/art_skewer
# Note: I worry that SourceForge is going to die any day now, so I chose to 
# mirror the installation files here install them into the container from a 
# local directory (hence the '-v $(pwd):/art' argument to docker build).
# Note: For me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM phusion/baseimage:latest
WORKDIR /tmp
ENV BUILD_PACKAGES \
    locales \
    build-essential \
    gcc-multilib \
    automake \
    intltool \
    pkg-config \
    libgsl-dev
ADD art_illumina_src151.tar.gz .
ADD art_illumina_src151-adapter-enabled.tar.gz .
RUN apt-get clean \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
        $BUILD_PACKAGES \
    && cd art_illumina_dir \
    && for f in config.sub config.guess install-sh depcomp missing INSTALL ; do rm -f $f ; ln -s /usr/share/automake-1.15/$f . ; done \
    && ./configure --prefix /usr/local \
    && make \
    && make install \
    && cd .. \
    && locale-gen en_US.UTF-8 \
    && update-locale LANG=en_US.UTF-8 \
    && apt-get remove --purge -y $BUILD_PACKAGES $(apt-mark showauto) \
    && rm -rf /var/lib/apt/lists/* /tmp/*
