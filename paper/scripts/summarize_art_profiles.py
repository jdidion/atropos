#!/usr/bin/env
from argparse import ArgumentParser
import csv
import sys
import editdistance
from common import *

def summarize_art_profile(p, w, i, adapter):
    while True:
        line = next(p)
        if line.startswith('#') or line.startswith('@'): continue
        break

    total_reads = 0
    num_reads_mismatch = 0
    total_edit_dist = 0
    num_adapters = 0
    total_ref_bp = 0
    total_adapter_bp = 0
    total_adapter_edit_dist = 0
    
    while True:
        if total_reads > 0:
            try:
                next(p)
            except:
                break
        r1 = next(p).rstrip()
        r2 = next(p).rstrip()
        total_reads += 1
        total_ref_bp += len(r1)
        
        if r1 == r2:
            continue
        
        if len(r1) != len(r2):
            num_adapters += 1
            l1 = len(r1)
            n = min(len(r2) - l1, len(adapter))
            total_adapter_bp += n
            total_adapter_edit_dist += editdistance.eval(r2[l1:(l1+n)], adapter[0:n])
            r2 = r2[0:l1]
        
        if r1 != r2:
            num_reads_mismatch += 1
            total_edit_dist += editdistance.eval(r1,r2)
    
    w.writerow((i,
        num_adapters, num_reads_mismatch, total_reads,
        total_edit_dist, total_ref_bp,
        total_adapter_edit_dist, total_adapter_bp
    ))

def main():
    parser = ArgumentParser()
    parser.add_argument('-a', '--adapters', nargs=2, defaults=DEFAULT_ADAPTERS)
    parser.add_argument('-1', '--profile1')
    parser.add_argument('-2', '--profile2')
    parser.add_argument('-o', '--output', default='-')
    args = parser.parse_args()
    
    with open(args.profile1, 'rt') as p1, open(args.profile2, 'rt') as p2, fileoutput(args.output) as o:
        w = csv.writer(o, delimiter="\t")
        w.writeline(('adapter',
            'num_adapters', 'num_reads_mismatch', 'total_reads',
            'total_edit_dist', 'total_ref_bp',
            'total_adapter_edit_dist', 'total_adapter_bp'
        ))
        summarize_art_profile(p1, w, 1, args.adapters[0])
        summarize_art_profile(p2, w, 2, args.adapters[1])

if __name__ == "__main__":
    main()
