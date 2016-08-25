# Rename all fastq files to have a .[12].fq.gz extension
for file in results/skewer*1.fastq.gz
do
    base="${file%%1.*}"
    mv $file $base.1.fq.gz
done
for file in results/skewer*2.fastq.gz
do
    base="${file%%2.*}"
    mv $file $base.2.fq.gz
done
