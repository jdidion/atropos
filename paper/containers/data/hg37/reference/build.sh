GENCODE_VERSION=26
FILE="gencode.v${GENCODE_VERSION}lift37.annotation.gtf.gz"
ANNOURL="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh37_mapping/${FILE}"
# Download the reference chromosomes and concat into single file
if [ ! -f hg37.fa ]
do
  URL='ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes'
  wget -r --no-parent --no-directories -A 'chr*' $URL && \
   zcat *.gz > hg37.fa && rm *.gz && \
   docker run -v $(pwd):/data --rm jdidion/bwa samtools faidx /tmp/hg37.fa \
done && \
if [ ! -f ${FILE} ] \
do \
  wget --no-parent --no-directories -O ${FILE} ${ANNOURL} && \
  gunzip ${FILE} \
done && \
# Build data container
docker build -f Dockerfile -t jdidion/hg37_reference:latest . && \
rm -Rf hg37.fa*

