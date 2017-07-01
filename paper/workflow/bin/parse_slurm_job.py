#!/usr/bin/env python
"""Parse job info from sacct. The output is two lines - the header row and the
job info row with the largest MaxRSS value.
"""
import argparse
from common import fileopen, parse_profile, parse_size
import csv
from subprocess import check_output

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("-p", "--profile", nargs='+')
    args = parser.parse_args()
    
    profile = list(parse_profile(args.profile[0]))
    if len(args.profile) > 1:
        profile.append(args.profile[1])
                
    with fileopen(args.input, 'rt') as inp:
        lines = list(csv.reader(inp, delimiter="|"))
    
    assert len(lines) >= 2
    header = lines[0]
    maxrss_col = header.index('MaxRSS')
    maxrss = 0
    maxrow = None
    for row in lines[1:]:
        row_maxrss = parse_size(row[maxrss_col])
        if row_maxrss > maxrss:
            maxrss = row_maxrss
            maxrow = row
    
    with fileopen(args.output, 'wt') as out:
        print(
            *('prog', 'prog2', 'threads', 'dataset', 'qcut'), *header,
            sep="\t", file=out)
        print(*profile, *maxrow, sep="\t", file=out)

if __name__ == '__main__':
    main()
