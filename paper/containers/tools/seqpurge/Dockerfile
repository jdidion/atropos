#################################################################
# Dockerfile
#
# Software:         SeqPurge
# Software Version: latest
# Description:      SeqPurge image
# Website:          https://github.com/imgag/ngs-bits
# Provides:         SeqPurge
# Base Image:       ubuntu:14.04
# Build Cmd:        docker build -f Dockerfile -t imgag/ngs-bits:latest .
#                   docker tag imgag/ngs-bits:latest jdidion/seqpurge
# Pull Cmd:         docker pull jdidion/seqpurge
# Run Cmd:          docker run
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM ubuntu:16.04
WORKDIR /tmp
RUN apt-get update          &&  \
    DEBIAN_FRONTEND=noninteractive apt-get -y upgrade      &&  \
    DEBIAN_FRONTEND=noninteractive apt-get install -y -qq      \
    software-properties-common \
    git                     \
    cmake                   \
    build-essential         \
    g++                     \
    libbz2-dev              \
    zlib1g-dev              \
    qt5-default             \
    qt5-qmake               \
    libqt5xmlpatterns5-dev  \
    libqt5sql5-mysql
RUN git clone --recursive https://github.com/imgag/ngs-bits.git && \
    cd ngs-bits && \
    make build_3rdparty && \
    make build_tools_release && \
    mv /tmp/ngs-bits/bin/* /usr/local/bin
# cleanup build dependencies and then reinstall runtime dependencies
RUN apt-get remove --purge -y $BUILD_PACKAGES $(apt-mark showauto)
RUN apt-get install -y qt5-default libqt5xmlpatterns5-dev
RUN rm -rf /var/lib/apt/lists/* /tmp/*
# Entrypoints don't work well with NextFlow - it wants to execute commands as
# if they were in the path.
#ENTRYPOINT ["/usr/local/bin/SeqPurge"]
#CMD ["--help"]
