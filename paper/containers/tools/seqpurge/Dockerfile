#################################################################
# Dockerfile
#
# Software:         SeqPurge
# Software Version: latest
# Description:      SeqPurge image
# Website:          https://github.com/imgag/ngs-bits
# Provides:         SeqPurge
# Base Image:       haffmans/cmake-qt5
# Build Cmd:        docker build -f Dockerfile -t imgag/ngs-bits:latest .
#                   docker tag imgag/ngs-bits:latest jdidion/seqpurge
# Pull Cmd:         docker pull jdidion/seqpurge
# Run Cmd:          docker run
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
#################################################################
FROM haffmans/cmake-qt5
WORKDIR /tmp
# install deps
RUN apt-get -y install g++ libqt5xmlpatterns5-dev libqt5sql5-mysql
# install ngs-bits
RUN git clone --recursive https://github.com/imgag/ngs-bits.git && \
    cd ngs-bits && \
    make build_3rdparty && \
    make build_tools_release && \
    mv /ngs-bits/bin/SeqPurge /usr/local/bin && \
    rm -Rf /tmp/ngs-bits
# cleanup
RUN apt-get remove --purge -y $BUILD_PACKAGES $(apt-mark showauto) && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/*
# Entrypoints don't work well with NextFlow - it wants to execute commands as
# if they were in the path.
#ENTRYPOINT ["/usr/local/bin/SeqPurge"]
#CMD ["--help"]
