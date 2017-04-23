# Download the reference genome
GENCODE_VERSION=26
FILE="gencode.v${GENCODE_VERSION}.annotation.gtf.gz"
ANNOURL="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/${FILE}"
URL='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
wget -O hg38.fa.gz $URL && \
 gunzip hg38.fa.gz && \
 docker run --rm jdidion/bwa samtools faidx hg38.fa && \
wget --no-parent --no-directories -O ${FILE} ${ANNOURL} && \
gunzip ${FILE} && \
# Build data container
docker build -f Dockerfile -t jdidion/hg38_reference:latest . &&  \
rm hg38.fa*
