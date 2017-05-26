#!/usr/bin/env python
"""Parse machine info (cat /proc/cpuinfo /proc/meminfo) and write a 
tab-delimited row to stdout with the following fields:

1. Program name
2. Total memory
3. Number of processors
4. Speed of processors
"""
import argparse
from collections import defaultdict
from common import fileopen, parse_profile
import re

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("-p", "--profile", nargs='+')
    args = parser.parse_args()
    
    profile = list(parse_profile(args.profile[0]))
    if len(args.profile) > 1:
        profile.append(args.profile[1])
    
    with fileopen(args.input, 'rt') as i:
        lines = [line.strip() for line in i.readlines()]
    
    mem_matcher = re.compile("MemTotal:\s+(\d+ .*)")
    mem = None
    
    cpu_matcher = re.compile("model name\s*:\s*(.*)")
    cpus = defaultdict(int)
    
    for line in lines:
        if mem is None:
            mem_match = mem_matcher.match(line)
            if mem_match:
                mem = mem_match.group(1)
                continue
        
        cpu_match = cpu_matcher.match(line)
        if cpu_match:
            cpus[cpu_match.group(1)] += 1
    
    total_cpus = sum(cpus.values())
    cpu_str = "; ".join("{} {}".format(count, cpu) for cpu, count in cpus.items())
    
    with fileopen(args.output, "wt") as out:
        print(*profile, total_cpus, mem, cpu_str, sep="\t", file=out)

if __name__ == "__main__":
    main()