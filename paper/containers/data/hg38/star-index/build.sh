#!/bin/bash
# Create a data volume from the reference genome image and run star index from 
# starbase image.
#
# Note: If using Docker, don't forget to increase the memory allocation to
# more than 16 GB for this command.

THREADS=$1
GENCODE_VERSION=26
INDEX_CMD="STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /index-build/hg38 \
        --genomeFastaFiles /index-build/hg38.fa \
        --sjdbGTFfile /data/annotations/hg38/gencode.v${GENCODE_VERSION}.annotation.gtf \
        --sjdbOverhang 75 --limitGenomeGenerateRAM 16000000000"

docker create \
  -v /data/reference/hg38 \
  -v /data/annotations/hg38 \
  --name hg38 jdidion/hg38_reference \
&& mkdir -p index \
&& docker run --rm \
  -v $(pwd)/index:/index-build --volumes-from hg38 \
  jdidion/starbase bash -c \
    "mkdir -p /data/index/star/hg38 && \
     cp /data/reference/hg38/hg38.fa /index-build/hg38 && \
     ${INDEX_CMD}" \
&& docker build -f Dockerfile -t jdidion/star_hg38index . \
&& docker rm -v hg38 \
&& rm -Rf index