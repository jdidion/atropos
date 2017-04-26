#################################################################
# Dockerfile
#
# Software:         skewer
# Software Version: latest
# Description:      skewer image
# Website:          https://github.com/relipmoc/skewer
# Provides:         skewer
# Base Image:       alpine
# Build Cmd:        docker build -f Dockerfile -t relipmoc/skewer:latest .
#                   docker tag relipmoc/skewer:latest jdidion/seqpurge
# Pull Cmd:         docker pull jdidion/skewer
# Run Cmd:          docker run
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM alpine
WORKDIR /tmp
# install deps
RUN apk add --no-cache \
    g++  \
    make \
    git  \
    bash
# install skewer
RUN git clone --recursive https://github.com/relipmoc/skewer && \
    cd skewer && \
    make && \
    mv /tmp/skewer/skewer /usr/local/bin && \
    rm -Rf /tmp/skewer
# cleanup
# TODO: remove deps that are not necessary at runtime
RUN rm -rf /var/cache/apk/* /tmp/*
# Entrypoints don't work well with NextFlow - it wants to execute commands as
# if they were in the path.
#ENTRYPOINT ["/usr/local/bin/skewer"]
#CMD ["--help"]
