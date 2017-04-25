# Note: The index commands require ~32G free memory -- more than we have on a
# desktop -- and Singularity can't mount Docker data volumes. Thus, we don't
# currently use this script. Instead, we copy the reference fasta and annotation
# file to the cluster rather than getting it from the container. The Singularity 
# command is:
# singularity run -H $(pwd) docker://jdidion/star /usr/local/bin $INDEX_CMD 8

THREADS=$1
GENCODE_VERSION=26
INDEX_CMD="STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /data/index/hg37 --genomeFastaFiles /data/index/hg37/hg37.fa \
        --sjdbGTFfile /data/annotations/hg37/gencode.v${GENCODE_VERSION}lift37.annotation.gtf \
        --sjdbOverhang 75 --limitGenomeGenerateRAM 16000000000"

# create a data volume from the reference genome image
docker create -v /data/reference/hg37 --name hg37 jdidion/hg37_reference && \
# build the STAR index
mkdir index && \
docker run \
    # create a local volume to store the output
    -v $(pwd)/index:/star-index:/data/index/star/hg37 \
    # bind reference data volume
    --volumes-from hg37 \
    # run bwa mem from bwa image
    --rm jdidion/star bash -c \
    "cp /data/reference/hg37/hg37.fa /data/index/star/hg37 && ${INDEX_CMD}" && \
# create a new image that includes the index
docker build -f Dockerfile -t jdidion/starindex . && \
# cleanup
rm -Rf star-index
