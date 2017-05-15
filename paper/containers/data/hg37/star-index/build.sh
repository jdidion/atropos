# Note: The index commands require ~32G free memory -- more than we have on a
# desktop -- and Singularity can't mount Docker data volumes. Thus, we don't
# currently use this script. Instead, we copy the reference fasta and annotation
# file to the cluster rather than getting it from the container. The Singularity 
# command is:
# singularity run -H $(pwd) docker://jdidion/star /usr/local/bin $INDEX_CMD 8

THREADS=$1
GENCODE_VERSION=26
INDEX_CMD="STAR --runMode genomeGenerate --runThreadN ${THREADS} \
        --genomeDir /index-build \
        --genomeFastaFiles /index-build/hg37.fa \
        --sjdbGTFfile /data/annotations/hg37/gencode.v${GENCODE_VERSION}lift37.annotation.gtf \
        --sjdbOverhang 75 --limitGenomeGenerateRAM 16000000000"

docker create -v /data/reference/hg37 --name hg37 jdidion/hg37_reference \
&& mkdir -p index \
&& docker run --rm \
    -v $(pwd)/index:/index-build --volumes-from hg37 \
    jdidion/starbase bash -c \
      "cp /data/reference/hg37/hg37.fa /index-build && \
       ${INDEX_CMD}" \
&& docker build -f Dockerfile -t jdidion/star_hg38index . \
&& docker rm -v hg37 \
&& rm -Rf index
