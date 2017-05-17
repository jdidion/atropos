# Note: The index commands require ~32G free memory -- more than we have on a
# desktop -- and Singularity can't mount Docker data volumes. Thus, we don't
# currently use this script. Instead, we copy the reference fasta to the 
# cluster rather than getting it from the container. The Singularity command is:
# singularity run -H $(pwd) docker://jdidion/bwa /usr/local/bin bwameth.py index hg37.fa

# create a data volume from the reference genome image
docker create -v /data/reference/hg37 --name hg37 jdidion/hg37_reference && \
# build the bwa index
docker run \
    # create a local volume to store the output
    -v $(pwd):/data/index/bwameth/hg37 \
    # bind reference data volume
    --volumes-from hg37 \
    # run bwa mem from bwa image
    --rm jdidion/bwabase bash -c \
    "cp /data/reference/hg37/hg37.fa /data/index/bwa/hg37 && \
     bwameth.py index /data/index/bwameth/hg37/hg37.fa" && \
# create a new image that includes the index
docker build -f Dockerfile -t jdidion/bwameth_hg37index . && \
# cleanup
rm -Rf bwa*index

