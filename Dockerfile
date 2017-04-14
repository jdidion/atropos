FROM snakeego/cython # python3 Alpine image with cython

# GCC
RUN export NPROC=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || 1) \
    && apk --no-cache add --virtual build-deps \
        musl-dev \
        linux-headers \
        g++ \
        zlib-dev bzip2-dev xz-dev

RUN pip install tqdm pysam

ADD . /atropos/

RUN cd /atropos/ && make install

ENTRYPOINT ["/atropos/bin/atropos"]
CMD ["--help"]
