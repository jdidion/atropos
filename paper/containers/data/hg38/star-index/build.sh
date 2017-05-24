#!/bin/bash
# Create a data volume from the reference genome image and run star index from 
# starbase image.
#
# Note: If using Docker, don't forget to increase the memory allocation to
# more than 16 GB for this command.
#
# Note: I find it most convenient to save the final image to a tar file that
# I can then copy to the HPC environment, rather than having to build the
# index twice. It's also nice to have a backup in case Docker or your computer
# crashes.

THREADS=$1
GENCODE_VERSION=26
INDEX_CMD="STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /index-build \
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
    "cp /data/reference/hg38/hg38.fa /index-build && \
     ${INDEX_CMD}" \
&& docker build -f Dockerfile -t jdidion/star_hg38index . \
&& docker rm -v hg38 \
&& rm -Rf index
