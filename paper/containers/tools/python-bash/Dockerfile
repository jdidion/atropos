#################################################################
# Dockerfile
#
# Software:         python
# Software Version: 3.6
# Description:      Python 3.6 with bash
# Website:          https://hub.docker.com/_/python/
# Provides:         python
# Base Image:       alpine
# Build Cmd:        docker build -f Dockerfile -t jdidion/python_bash .
# Pull Cmd:         docker pull jdidion/python_bash
# Run Cmd:          docker run --rm jdidion/python_bash python -h
#################################################################
FROM python:3-alpine
RUN apk add --no-cache libstdc++ && \
    apk add --no-cache \
        --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community \
        lapack-dev && \
    apk add --no-cache --virtual=.build-dependencies \
        pkgconfig     \
        g++           \
        gfortran      \
        musl-dev      \
        python3-dev   \
        libpng-dev    \
        freetype-dev && \
    apk add --no-cache \
        bash        \
        libpng      \
        freetype && \
    ln -s locale.h /usr/include/xlocale.h && \
    pip install --upgrade pip && \
    pip install cython && \
    pip install numpy && \
    pip install pandas && \
    pip install matplotlib && \
    pip install seaborn && \
    pip install mako && \
    pip uninstall --yes cython && \
    rm /usr/include/xlocale.h && \
    rm -r /root/.cache && \
    apk del .build-dependencies
