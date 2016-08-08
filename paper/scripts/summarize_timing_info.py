#!/usr/bin/env python
# Summarize the timiing information output by the command scripts into a nicely formatted tsv.

from argparse import ArgumentParser
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

def summarize_timing(i, o):
    pass

def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', default='-')
    parser.add_argument('-o', '--output', default='-')
    args = parser.parse_args()
    
    with fileinput.input(args.input) as i, fileoutput(args.output) as o:
        summarize_timing(i, o)

if __name__ == "__main__":
    main()
