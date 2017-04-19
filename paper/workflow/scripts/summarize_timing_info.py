#!/usr/bin/env python
"""Summarize the timing information output by the command scripts into a nicely
formatted tsv.

For analyses run locally, there will be a input file:

    python summarize_timing_info.py -i timing_log.txt -o timing_summary.txt

For analyses run on the cluster, there will be one stderr log for each process:

    cat commands_t16.e* | python summarize_timing_info.py -i - -o timing_summary.txt
"""

from argparse import ArgumentParser
from collections import defaultdict
import csv
import fileinput
from common import fileoutput

def summarize_timing(lines):
    rows = []
    for line in lines:
        if any(line.startswith(prog) for prog in ('atropos', 'skewer', 'seqpurge')):
            profile = line.rstrip().split("_")
            line = next(lines)
            assert line.startswith('real')
            time = float(line.rstrip()[5:])
            if profile[0] == 'atropos':
                if len(profile) == 7:
                    prog, threads, err1, err2, qcut, aligner, writer = profile
                    err = '{}_{}'.format(err1, err2)
                else:
                    prog, threads, err, qcut, aligner, writer = profile
                prog = "{} ({} + {})".format(profile[0], aligner, writer)
            elif len(profile) == 5:
                prog, threads, err1, err2, qcut = profile
                err = '{}_{}'.format(err1, err2)
            else:
                prog, threads, err, qcut = profile
            err = '0.' + err
            qcut = qcut[1:]
            rows.append((prog, threads, err, qcut, time))
            
            # throw away user and sys times - not super informative, especially with multi-threaded programs
            line = next(lines)
            assert line.startswith('user')
            line = next(lines)
            assert line.startswith('sys')
    
    return rows

def format_table(raw_rows):
    table = defaultdict(lambda: defaultdict(lambda: []))
    threads = set()
    for i, row in enumerate(raw_rows):
        t = int(row[1])
        table[row[0]][t].append(row[4])
        threads.add(t)
    rows = []
    threads = list(sorted(threads))
    for key in sorted(table.keys()):
        row = [key]
        for t in threads:
            if t in table[key]:
                times = table[key][t]
                min_time = min(times)
                max_time = max(times)
                row.append("{} - {}".format(min_time, max_time))
            else:
                row.append("NA")
        rows.append(row)
    return rows, threads

def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', default='-')
    parser.add_argument('-o', '--output', default='-')
    parser.add_argument('-f', '--output-format', choices=('simple', 'latex'), default='simple')
    parser.add_argument('-n', '--table-name', default="timing")
    parser.add_argument('-c', '--table-caption', default="Timing results")
    args = parser.parse_args()
    
    with fileinput.input(args.input) as i:
        rows = summarize_timing(i)
    
    if len(rows) == 0:
        print("No results to summarize")
    
    with fileoutput(args.output) as o:
        w = csv.writer(o, delimiter="\t")
        if args.output_format == 'simple':
            w.writerow(('Program', 'Threads', 'Error Rate', 'Quality Cutoff', 'Time'))
            w.writerows(rows)
        else:
            from mako.template import Template
            template = Template(filename='timing_table_template.latex')
            rows, threads = format_table(rows)
            o.write(template.render(rows=rows, threads=threads, name=args.table_name, caption=args.table_caption))

if __name__ == "__main__":
    main()
