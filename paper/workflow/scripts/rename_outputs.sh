# Rename all fastq files to have a .[12].fq.gz extension
for file in results/**/skewer*trimmed-pair1.fastq.gz
do
    base="${file%%-trimmed-pair1.*}"
    mv $file $base.1.fq.gz
done
for file in results/**/skewer*trimmed-pair2.fastq.gz
do
    base="${file%%-trimmed-pair2.*}"
    mv $file $base.2.fq.gz
done
