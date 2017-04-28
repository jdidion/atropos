#################################################################
# Dockerfile
#
# Software:         ART (modified)
# Software Version: 1.5.1
# Description:      ART read simulator, modified to introduce adapter sequences
# Website:          https://sourceforge.net/projects/skewer/files/Simulator/
# Provides:         ART
# Base Image:       jedisct1/phusion-baseimage-latest
# Build Cmd:        docker build -f Dockerfile -t jdidion/art_skewer:latest .
# Pull Cmd:         docker pull jdidion/art_skewer
# Run Cmd:          docker run jdidion/art_skewer
# Note: I worry that SourceForge is going to die any day now, so I chose to 
# mirror the installation files here install them into the container from a 
# local directory (hence the '-v $(pwd):/art' argument to docker build).
# Note: For me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM jedisct1/phusion-baseimage-latest
WORKDIR /tmp
# Install compiler and perl stuff
RUN apt-get update && \
    apt-get install --yes \
        apt-utils \
        build-essential \
        gcc-multilib \
        automake \
        intltool \
        pkg-config \
        libgsl-dev \
        bash
# Install modified ART
ADD art_illumina_src151.tar.gz .
ADD art_illumina_src151-adapter-enabled.tar.gz .
RUN cd art_illumina_dir && \
    for f in config.sub config.guess install-sh depcomp missing INSTALL ; do rm -f $f ; ln -s /usr/share/automake-1.15/$f . ; done && \
    ./configure --prefix /usr/local && \
    make && \
    make install
# Cleanup
RUN rm -rf /tmp/*