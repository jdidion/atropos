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
for err in 005 01
do
  for mate in 1 2
  do
    Rscript --vanilla adjust_error_profiles.R -g ../data/HiSeq2500L125R${mate}_001.txt ../data/HiSeq2500L125R${mate}_${err}.txt 0.${err}
  done
done

# Variables used in read simulation
seed=42
read_len=125
mean_frag=200
sd_frag=100
cov=0.1

# Simulate reads
for prof in 001 005 01
do
  # NOTE: REF_GENOME is externally defined to point to the fasta file for hg19/GRCh37.
  # This assumes that the modified version of ART has been downloaded and compiled in the bin folder.
  ../software/bin/art_illumina -i $REF_GENOME -p -l $read_len -f $cov -m $mean_frag -s $sd_frag -rs $seed -o sim_${prof} \
    -1 ../data/HiSeq2500L125R1_${prof}.txt \
    -2 ../data/HiSeq2500L125R2_${prof}.txt
  # move and rename files
  for mate in 1 2
  do
    for ext in aln fq
    do
      mv sim_${prof}${mate}.${ext} ../data/simulated/sim_${prof}.${mate}.${ext}
    done
  done
done
