#!/bin/bash
# This script simulates reads using ART (Huang et al. 2012) that was modified to add adapter sequences
# (Jiang et al. 2014). We simulated ~1M 125 bp PE reads (0.1x genome coverage) using the empirically
# determined quality score profile (based on data derived from real HiSeq 2500 sequencing reads) -
# which has an overall error rate of 0.11% - and transformations of the empircal profile for which the
# quality score distribution was artificially left-shifited to result in higher overall error rates of
# 0.5% and 1%.
#
# References:
# Huang et al. 2012, http://bioinformatics.oxfordjournals.org/content/28/4/593.full
# Jiang et al. 2014, http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182

# Generate error-shifted read quality profiles
R_CMD="docker run --rm -v "$PWD":/home/docker -w /home/docker -u docker r-base Rscript --vanilla"
for err in 005 01
do
  for mate in 1 2
  do
    $R_CMD adjust_error_profiles.R -g -p art_profiles/HiSeq2500L125R${mate}_001.txt -o art_profiles/data/HiSeq2500L125R${mate}_${err}.txt -e 0.${err} -t art_profiles.txt
  done
done

# Variables used in read simulation
seed=42
read_len=125
mean_frag=200
sd_frag=100
cov=0.1

# prepare the reference genome data volume
docker create -v /data/reference/hg37 --name hg37 jdidion/hg37_reference
REF_GENOME="/data/reference/hg37/hg37.fa"
# Simulate reads
ART_CMD="docker run --rm -v "$PWD":/art --volumes-from hg37 -w /art jdidion/art_skewer art_illumina"
for prof in 001 005 01
do
  $ART_CMD -i $REF_GENOME -p -l $read_len -f $cov -m $mean_frag -s $sd_frag -rs $seed -o sim_${prof} \
    -1 art_profiles/HiSeq2500L125R1_${prof}.txt \
    -2 art_profiles/HiSeq2500L125R2_${prof}.txt
  # rename and fix files
  for mate in 1 2
  do
    for ext in aln fq
    do
      mv sim_${prof}${mate}.${ext} sim_${prof}.${mate}.${ext}
      # There is a bug in the modified ART in which reads that are shorter
      # than the read length even with the adapter apended have null characters
      # tacked on to the end. We replace them with 'A's.
      perl -pe 's/\x00/A/g' sim_${prof}.${mate}.${ext}
    done
  done
done

docker build -f Dockerfile -t jdidion/atropos_simulated .
