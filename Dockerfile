#################################################################
# Dockerfile
#
# Software:         Atropos
# Software Version: 1.1.0
# Description:      Atropos image
# Website:          https://github.com/jdidion/atropos
# Provides:         atropos
# Base Image:       snakeego/cython
# Build Cmd:        docker build -f Dockerfile -t jdidion/atropos:latest .
# Pull Cmd:         docker pull jdidion/atropos
# Run Cmd:          docker run --rm jdidion/atropos
# Note: This Dockerfile is for building a container from a copy of the
# git repository. If you want to build a container without cloning the
# repo, use the Dockerfile in paper/containers/tools/atropos-paper instead.
#################################################################
FROM snakeego/cython

# bash support to enable container usage from cromwell and other workflow engines
RUN apk add --update bash && rm -rf /var/cache/apk/*

# GCC
RUN export NPROC=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || 1) \
    && apk --no-cache add --virtual build-deps \
        musl-dev \
        linux-headers \
        g++ \
        zlib-dev bzip2-dev xz-dev

# Additional Atropos dependencies
RUN pip install tqdm pysam jinja2

# Attach project directory and install
ADD . /atropos/
RUN cd /atropos/ && make install

# Cleanup
RUN rm -rf /var/cache/apk/* /atropos

# Entrypoints don't work well with NextFlow - it wants to execute commands as
# if they were in the path.
ENTRYPOINT ["/usr/local/bin/atropos"]
CMD ["--help"]
