#!/usr/bin/env python
# Summarize the timiing information output by the command scripts into a nicely formatted tsv.
#
# For analyses run locally, there will be a input file:
#
#   python summarize_timing_info.py -i timing_log.txt -o timing_summary.txt
#
# For analyses run on the cluster, there will be one stderr log for each process:
#
#   cat commands_t*.e* | python summarize_timing_info.py -i - -o timing_summary.txt

from argparse import ArgumentParser
import csv
import fileinput

class fileoutput(object):
    def __init__(path, mode='wt'):
        if path == '-':
            self.fh = sys.stdout
            close = False
        else:
            self.fh = open(path, mode)
            close = True
    
    def __enter__(self):
        return self.fh
    
    def __exit__(self, type, value, traceback):
        if self.close:
            self.fh.close()

def summarize_timing(i, w):
    for line in i:
        if any(i.startswith(prog) for prog in ('atropos', 'skewer', 'seqpurge')):
            profile = line.rstrip().split("_")
            line = next(i)
            assert line.startswith('real')
            time = float(line.rstrip()[5:])
            w.writerow((profile[0], profile[4] if len(profile) > 4 else "") + profile[1:4] + (time,))
            # throw away user and sys times - not super informative, especially with multi-threaded programs
            line = next(i)
            assert line.startswith('user')
            line = next(i)
            assert line.startswith('sys')

def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', default='-')
    parser.add_argument('-o', '--output', default='-')
    args = parser.parse_args()
    
    with fileinput.input(args.input) as i, fileoutput(args.output) as o:
        w = csv.writer(o, separator="\t")
        summarize_timing(i, w)

if __name__ == "__main__":
    main()
