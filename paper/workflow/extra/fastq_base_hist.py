# Generate base histograms of the first and last 20 bp of read1 and read2

import argparse
from common import *
import csv
from xphyle import xopen
import tqdm

nuc = ('A','C','G','T','N')

def make_hists(fq1, fq2, n_bases=20):
    read1_hist = [
        dict((base, [0] * n_bases) for base in nuc)
        for i in range(2)
    ]
    read2_hist = [
        dict((base, [0] * n_bases) for base in nuc)
        for i in range(2)
    ]
    
    def add_to_hist(r, hist):
        for i, b in enumerate_range(r[2].upper(), 0, n_bases):
            hist[0][b][i] += 1
        for i, b in enumerate_range(reversed(r[2].upper()), 0, n_bases):
            hist[1][b][i] += 1
    
    for r1, r2 in tqdm.tqdm(zip(fq_iterator(fq1), fq_iterator(fq2))):
        add_to_hist(r1, read1_hist)
        add_to_hist(r2, read2_hist)
    
    return (read1_hist, read2_hist)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", "--fastq1")
    parser.add_argument("-2", "--fastq2")
    parser.add_argument("-o", "--output", default="-")
    args = parser.parse_args()
    
    with xopen(args.fastq1) as fq1, xopen(args.fastq2) as fq2:
        hists = make_hists(fq1, fq2)
    
    with xopen(args.output, 'w') as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(('read', 'side', 'pos', 'base', 'count'))
        for i, h in enumerate(hists, 1):
            for j in range(2):
                for b in nuc:
                    for k, count in enumerate(h[j][b], 1):
                        w.writerow((i, j, k, b, count))

if __name__ == "__main__":
    main()
