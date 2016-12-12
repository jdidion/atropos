"""Summarize timing and memory usage from SGE result emails.
"""
from argparse import ArgumentParser, FileType
from collections import defaultdict
import csv
import sys
import re

name_re = re.compile('Job-array task \d+\.(\d+) \((.+)\) Complete')
time_re = re.compile('Wallclock Time.*= (\d+):(\d+):(\d+)')
mem_re = re.compile('Max vmem.*= (.*)')

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-i', '--infile', type=FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    
    results = defaultdict(lambda: {})
    
    lines = args.infile.readlines()
    for i in range(0, len(lines), 18):
        name = name_re.match(lines[i+5])
        if not name:
            raise Exception("failed job: {}".format(lines[i+5]))
        job_index, job_array = name.groups()
        hr,mn,sec = (int(x) for x in time_re.match(lines[i+13]).groups())
        mem = mem_re.match(lines[i+15]).group(1)
        
        results[job_array][int(job_index)] = (hr,mn,sec,mem)
    
    writer = csv.writer(args.outfile, delimiter="\t")
    writer.writerow(('job_array','job_index','hrs','mins','secs','mem'))
    for job_array, jobs in results.items():
        job_indexes = sorted(jobs.keys())
        for i in job_indexes:
            writer.writerow((job_array, i) + jobs[i])
