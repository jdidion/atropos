#!/usr/bin/env python
"""Parse the output of GNU time (gtime), along with a profile name ($item)
and write a tab-delimited row to stdout with the following fields:

1. Program name
2. Number of threads
3. The dataset - either an error rate (for simulated data) or a data type (WGBS or RNA-Seq)
4. Quality cutoff
5. Duration of program execution
6. Max CPU usage
7. Max memory usage
"""
import argparse
from common import fileopen
import os
import re
import sys
import fileinput

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default="-")
parser.add_argument("-o", "--output", default="-")
parser.add_argument("-p", "--profile")
args = parser.parse_args()

profile = args.profile.split("_")

if profile[0] == 'atropos':
    if len(profile) == 7:
        prog, threads, dataset1, dataset2, qcut, aligner, writer = profile
        dataset = '{}_{}'.format(dataset1, dataset2)
    else:
        prog, threads, dataset, qcut, aligner, writer = profile
    prog = "{} ({})".format(profile[0], aligner)
    prog2 = "{} ({} + {})".format(profile[0], aligner, writer)
else:
    if len(profile) == 5:
        prog, threads, dataset1, dataset2, qcut = profile
        dataset = '{}_{}'.format(dataset1, dataset2)
    else:
        prog, threads, dataset, qcut = profile
    prog2 = prog
# If this is a simulated dataset, format the error rate
try:
    float(dataset)
    dataset = '0.' + dataset
except:
    pass
qcut = qcut[1:]

with fileopen(args.input, 'rt') as i:
    lines = [line.strip() for line in i.readlines()]

cpu = lines[3]
cpu_match = re.match("Percent of CPU this job got: (\\d+)%", cpu)
assert cpu_match is not None
cpu_pct = float(cpu_match.group(1))

wc_time = lines[4]
wc_match_prefix = "Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): "
wc_match = re.match(wc_match_prefix + "(?:(\\d+)h )?(\\d+)m ([\\d\\.]+)s", wc_time)
if wc_match is None:
    wc_match = re.match(wc_match_prefix + "(?:(\\d+):)?(\\d+):([\\d\\.]+)", wc_time)
assert wc_match is not None
hrs = int(wc_match.group(1) or 0)
mins = int(wc_match.group(2))
secs = float(wc_match.group(3))
duration = ':'.join((
    '{:02d}'.format(hrs),
    '{:02d}'.format(mins),
    '{:0.2f}'.format(secs)))
duration_secs = (hrs * 3600) + (mins * 60) + secs

memory = lines[9]
memory_match = re.match("Maximum resident set size \\(kbytes\\): (\\d+)", memory)
assert memory_match is not None
memory_mbytes = int(memory_match.group(1)) / 1000

with fileopen(args.output, "wt") as out:
    print(prog, prog2, threads, dataset, qcut, duration_secs,
    cpu_pct, memory_mbytes, sep="\t", file=out)
