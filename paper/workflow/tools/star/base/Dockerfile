#################################################################
# Dockerfile
#
# Software:         STAR
# Software Version: 2.5.3a
# Description:      Alpine STAR RNA-Seq aligner Image
# Website:          https://github.com/alexdobin/STAR
# Provides:         STAR
# Base Image:       alpine
# Build Cmd:        docker build -t jdidion/star:latest .
# Pull Cmd:         docker pull jdidion/star
# Run Cmd:          docker run --rm jdidion/star STAR
#################################################################
FROM alpine
ENV DEBIAN_FRONTEND noninteractive
WORKDIR /tmp

# GCC
RUN apk add --no-cache g++ make git

# STAR
ENV VERSION 2.5.3a
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
