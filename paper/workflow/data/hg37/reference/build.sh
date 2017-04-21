# Download the reference chromosomes and concat into single file
if [ ! -f hg37.fa ]
do
  URL='ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes'
  wget -r --no-parent --no-directories -A 'chr*' $URL && \
   zcat *.gz > hg37.fa && rm *.gz && \
   docker run -v $(pwd):/data --rm jdidion/bwa samtools faidx /data/hg37.fa
done && \
# Build data container
docker build -f Dockerfile -t jdidion/hg37_reference:latest . && \
rm -Rf hg37.fa*
