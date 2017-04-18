FROM haffmans/cmake-qt5 # Ubuntu 14.04 image with qt5 and cmake

# GCC, git
RUN apt-get install g++ git libqt5xmlpatterns5-dev libqt5sql5-mysql
# Clone ngs-bits
RUN git clone --recursive https://github.com/imgag/ngs-bits.git && \
    cd ngs-bits &&
    make build_3rdparty &&
    make build_tools_release

ENTRYPOINT ["/ngs-bits/bin/SeqPurge"]
CMD ["--help"]
