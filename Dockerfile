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
# Run Cmd:          docker run
#################################################################
FROM snakeego/cython

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
RUN rm -rf /var/cache/apk/* /atropos && apk del deps

ENTRYPOINT ["/atropos/bin/atropos"]
CMD ["--help"]
