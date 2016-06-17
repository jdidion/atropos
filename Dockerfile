FROM ubuntu-core:alpha-3

RUN snappy install python3 cython

ADD . /atropos/

RUN cd /atropos/ && python setup.py install && python setup.py build_ext -i

ENTRYPOINT ["/atropos/bin/atropos"]
CMD ["--help"]

# git clone https://github.com/jdidion/atropos.git
# cd atropos
# docker build -t jdidion/atropos:latest .
