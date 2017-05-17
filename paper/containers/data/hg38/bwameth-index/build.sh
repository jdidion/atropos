# Note: The index commands require ~32G free memory -- more than we have on a
# desktop -- and Singularity can't mount Docker data volumes. Thus, we don't
# currently use this script. Instead, we copy the reference fasta to the 
# cluster rather than getting it from the container. The Singularity command is:
# singularity run -H $(pwd) docker://jdidion/bwa /usr/local/bin bwameth.py index hg38.fa

# create a data volume from the reference genome image
#docker create -v /data/reference/hg38 --name hg38 jdidion/hg38_reference && \
# build the bwa index
# create a local volume to store the output bind reference data volume
# and run bwa mem from bwa image
# create a new image that includes the index
mkdir index \
&& docker run \
    -v $(pwd)/index:/data/index/bwameth/hg38 \
    --volumes-from hg38 \
    --rm jdidion/bwabase bash -c \
    "cp /data/reference/hg38/hg38.fa /data/index/bwameth/hg38 && \
     bwameth.py index /data/index/bwameth/hg38/hg38.fa" \
&& docker build -f Dockerfile -t jdidion/bwameth_hg38index . \
&& rm -Rf index
