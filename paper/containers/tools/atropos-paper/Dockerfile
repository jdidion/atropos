#################################################################
# Dockerfile
#
# Software:         Atropos
# Software Version: 1.1.0
# Description:      Atropos image for paper
# Website:          https://github.com/jdidion/atropos
# Provides:         atropos
# Base Image:       snakeego/cython
# Build Cmd:        docker build -f Dockerfile -t jdidion/atropos_paper:latest .
# Pull Cmd:         docker pull jdidion/atropos_paper
# Run Cmd:          docker run --rm jdidion/atropos_paper atropos -h
#################################################################
FROM snakeego/cython
WORKDIR /tmp

# GCC
RUN export NPROC=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || 1) \
    && apk --no-cache add --virtual build-deps \
        musl-dev \
        linux-headers \
        g++ \
        git \
        zlib-dev bzip2-dev xz-dev

# Additional Atropos dependencies
RUN pip install tqdm pysam jinja2

# Attach project directory and install
RUN git clone --recursive https://github.com/jdidion/atropos && \
    cd atropos && \
    make install

# Cleanup
RUN rm -rf /var/cache/apk/* /tmp/*

# Entrypoints don't work well with NextFlow - it wants to execute commands as
# if they were in the path.
#ENTRYPOINT ["/atropos/bin/atropos"]
#CMD ["--help"]
