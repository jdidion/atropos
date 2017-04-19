# Here's how to build and deploy.
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
# 1. docker build -f Dockerfile -t imgag/ngs-bits:latest .
# 2. docker tag imgag/ngs-bits:latest jdidion/ngs-bits
# 3. docker login -e <email> -u <user> -p <password>
# 4. docker push jdidion/ngs-bits

# Ubuntu 14.04 image with qt5 and cmake
FROM haffmans/cmake-qt5

# GCC, git
RUN apt-get -y install g++ libqt5xmlpatterns5-dev libqt5sql5-mysql

# Clone ngs-bits
RUN git clone --recursive https://github.com/imgag/ngs-bits.git && \
    cd ngs-bits && \
    make build_3rdparty && \
    make build_tools_release

ENTRYPOINT ["/ngs-bits/bin/SeqPurge"]
CMD ["--help"]
