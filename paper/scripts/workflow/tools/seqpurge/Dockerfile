FROM haffmans/cmake-qt5 # Ubuntu 14.04 image with qt5 and cmake

# GCC, git
RUN apt-get install g++ git
# Clone ngs-bits
RUN git clone --recursive https://github.com/imgag/ngs-bits.git
RUN pip install tqdm pysam

ADD . /atropos/

RUN cd /atropos/ && make install

ENTRYPOINT ["/atropos/bin/atropos"]
CMD ["--help"]
