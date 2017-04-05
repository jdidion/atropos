FROM ubuntu-core:alpha-3

RUN snappy install python3 cython python3-pip
RUN pip install tqdm pysam

ADD . /atropos/

RUN cd /atropos/ && make install

ENTRYPOINT ["/atropos/bin/atropos"]
CMD ["--help"]

# git clone https://github.com/jdidion/atropos.git
# cd atropos
# docker build -t jdidion/atropos:latest .
