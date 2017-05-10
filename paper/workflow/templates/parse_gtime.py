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
import os
import re
import sys
import fileinput

profile = os.environ.get("item", "$item").split("_")
timing_file = os.environ.get("timing", "$timing")

if profile[0] == 'atropos':
    if len(profile) == 7:
        prog, threads, dataset1, dataset2, qcut, aligner, writer = profile
        dataset = '{}_{}'.format(dataset1, dataset2)
    else:
        prog, threads, dataset, qcut, aligner, writer = profile
    prog = "{} ({} + {})".format(profile[0], aligner, writer)
elif len(profile) == 5:
    prog, threads, dataset1, dataset2, qcut = profile
    dataset = '{}_{}'.format(dataset1, dataset2)
else:
    prog, threads, dataset, qcut = profile
# If this is a simulated dataset, format the error rate
try:
    float(dataset)
    dataset = '0.' + dataset
except:
    pass
qcut = qcut[1:]

with open(timing_file, 'rt') as i:
    lines = [line.strip() for line in i.readlines()]

cpu = lines[3]
cpu_match = re.match("Percent of CPU this job got: (\\d+)%", cpu)
assert cpu_match is not None
cpu_frac = float(cpu_match.group(1))

wc_time = lines[4]
wc_match_prefix = "Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): "
wc_match = re.match(wc_match_prefix + "(?:(\\d+)h )?(\\d+)m ([\\d\\.]+)s", wc_time)
if wc_match is None:
    wc_match = re.match(wc_match_prefix + "(?:(\\d+):)?(\\d+):([\\d\\.]+)", wc_time)
assert wc_match is not None
duration = ':'.join((
    '{:02d}'.format(int(wc_match.group(1) or 0)),
    '{:02d}'.format(int(wc_match.group(2))),
    '{:0.2f}'.format(float(wc_match.group(3)))))

memory = lines[9]
memory_match = re.match("Maximum resident set size \\(kbytes\\): (\\d+)", memory)
assert memory_match is not None
memory_bytes = int(memory_match.group(1)) * 1000

print(prog, threads, dataset, qcut, duration, cpu_frac, memory_bytes, sep="\t")
