#################################################################
# Dockerfile
#
# Software:         STAR
# Software Version: 2.5.0b
# Description:      STAR image for SciDAP
# Website:          https://github.com/alexdobin/STAR, http://scidap.com/
# Provides:         STAR
# Base Image:       scidap/scidap:v0.0.1
# Build Cmd:        docker build --rm -t scidap/star:v2.5.0b .
# Pull Cmd:         docker pull scidap/star:v2.5.0b
# Run Cmd:          docker run --rm scidap/star:v2.5.0b STAR
#################################################################

### Base Image
FROM scidap/scidap:v0.0.1
MAINTAINER Andrey V Kartashov "porter@porter.st"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp

### Install STAR

ENV VERSION 2.5.0b
ENV NAME STAR
ENV URL "https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION}/source && \
    make -j 4 && \
    cd .. && \
    cp ./bin/Linux_x86_64_static/STAR /usr/local/bin/ && \
    cd .. && \
    strip /usr/local/bin/${NAME}; true && \
    rm -rf ./${NAME}-${VERSION}/