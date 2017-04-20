#################################################################
# Dockerfile
#
# Software:         AdapterRemoval
# Software Version: 2.2.0
# Description:      AdapterRemoval2 image
# Website:          https://github.com/MikkelSchubert/adapterremoval
# Provides:         AdapterRemoval
# Base Image:       alpine
# Build Cmd:        docker build -f Dockerfile -t mikkelschubert/adapterremoval:latest .
#                   docker tag mikkelschubert/adapterremoval:latest jdidion/adapterremoval
# Pull Cmd:         docker pull jdidion/adapterremoval
# Run Cmd:          docker run
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM alpine
WORKDIR /tmp
RUN apk add --no-cache g++ git make zlib-dev bzip2-dev perl
RUN git clone --recursive https://github.com/MikkelSchubert/adapterremoval && \
    cd adapterremoval && \
    make && \
    mv /tmp/adapterremoval/build/AdapterRemoval /usr/local/bin && \
    rm -Rf /tmp/adapterremoval
RUN rm -rf /var/cache/apk/* /tmp/* && \
    apk del deps
ENTRYPOINT ["/usr/local/bin/AdapterRemoval"]
CMD ["--help"]
