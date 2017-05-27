# Use fasatq-dump in atropos-tools to fetch some RNA-seq.
ACCN=SRR521458
docker run -v $(pwd):/dump --rm jdidion/atropos_paper_analysis bash -c \
  "/opt/fastq-dump/fastq-dump-wrapper.sh --split-files --gzip -O /dump ${ACCN}" && \
docker build -f Dockerfile -t jdidion/atropos_rnaseq:latest . && \
rm ${ACCN}_*.fastq.gz
