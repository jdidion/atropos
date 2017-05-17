# Note: The index commands require ~32G free memory -- more than we have on a
# desktop -- and Singularity can't mount Docker data volumes. Thus, we don't
# currently use this script. Instead, we copy the reference fasta to the 
# cluster rather than getting it from the container. The Singularity command is:
# singularity run -H $(pwd) docker://jdidion/bwa /usr/local/bin bwa index hg38.fa

# create a data volume from the reference genome image
docker create -v /data/reference/hg38 --name hg38 jdidion/hg38_reference && \
# build the bwa index
mkdir index && ]
docker run \
    # create a local volume to store the output
    -v $(pwd)/index:/data/index/bwa/hg38 \
    # bind reference data volume
    --volumes-from hg38 \
    # run bwa mem from bwa image
    --rm jdidion/bwabase bash -c \
    "cp /data/reference/hg38/hg38.fa /data/index/bwa/hg38 && \
     bwa mem index /data/index/bwa/hg38/hg38.fa" && \
# create a new image that includes the index
docker build -f Dockerfile -t jdidion/bwa_hg38index . && \
# cleanup
rm -Rf bwa*index

