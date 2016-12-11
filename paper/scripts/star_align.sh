name=$1
base=$2
threads=$3
genome=$4
outdir=`dirname $base`
fq1=${base}.1.fq.gz
fq2=${base}.2.fq.gz
mkdir -p $name
cd $name
mkfifo Read1 Read2
export TMPDIR='.'
echo "Aligning $fq1 $fq2"
/usr/bin/gzip -cd $fq1 > Read1 &
/usr/bin/gzip -cd $fq2 > Read2 &
# Don't exclude multi-mappers, just randomly select one
STAR --runThreadN $threads --genomeDir $genome --readFilesIn Read1 Read2 \
  --outMultimapperOrder Random --outFilterMultimapNmax 100000 --outSAMmultNmax 1 \
  --outFileNamePrefix ${outdir}/${name}_rna \
  --outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs
rm Read1 Read2
cd ..
rm -Rf $name
