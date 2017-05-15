#!/bin/bash
# Create a data volume from the reference genome image and run star index from 
# starbase image.
#
# Note: If using Docker, don't forget to increase the memory allocation to
# more than 16 GB for this command.

THREADS=$1
GENCODE_VERSION=26
INDEX_CMD="STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /data/index/star/hg38 \
        --genomeFastaFiles /data/index/star/hg38/hg38.fa \
        --sjdbGTFfile /data/annotations/hg38/gencode.v${GENCODE_VERSION}.annotation.gtf \
        --sjdbOverhang 75 --limitGenomeGenerateRAM 16000000000"

docker create \
  -v /data/reference/hg38 \
  -v /data/annotations/hg38 \
  --name hg38 jdidion/hg38_reference \
&& docker run \
    --volumes-from hg38 \
    --rm jdidion/starbase bash -c \
    "mkdir -p /data/index/star/hg38 && \
     cp /data/reference/hg38/hg38.fa /data/index/star/hg38 && \
     ${INDEX_CMD}" \
&& docker commit jdidion/starbase jdidion/star_hg38index \
&& docker rm -v hg38