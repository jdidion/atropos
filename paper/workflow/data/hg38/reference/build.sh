# Download the reference genome
URL='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
wget -O hg38.fa.gz $URL && \
 gunzip hg38.fa.gz && \
 docker run --rm jdidion/bwa samtools faidx hg38.fa && \
# Build data container
docker build -f Dockerfile -t jdidion/hg38_reference:latest . &&  \
rm hg38.fa*