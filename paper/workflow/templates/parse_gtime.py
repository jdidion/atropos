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
import sys
import fileinput

profile = $item.split("_")
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

lines = list(fileinput())

cpu = lines[3]
cpu_match = re.match("Percent of CPU this job got: (\d+)%", cpu)
assert cpu_match is not None
cpu_frac = float(cpu_match.group(1)) / 100

wc_time = lines[4]
wc_match = re.match(
    "Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (?:(\d+)h )?(\d+)m ([\d\.]+)s",
    wc_time)
assert wc_match is not None
duration = ':'.join(
    '{:2d}'.format(wc_match.group(1) or 0),
    '{:2d}'.format(wc_match.group(2)),
    '{:0.2f}'.format(wc_match.group(3)))

memory = lines[9]
memory_match = re.match("Maximum resident set size \(kbytes\): (\d+)", memory)
assert memory_match is not None
memory_bytes = int(memory_match.group(1)) * 1000

print(prog, threads, dataset, qcut, duration, cpu, memory, sep="\t")
