#!/usr/bin/env python
import argparse
import csv
import os
import pandas as pd
from common import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
sb.set(style="whitegrid")
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "-t", "--tool-name-file",
        help="File that maps profile names to display names for tools")
    parser.add_argument(
        "-T", "--threads",
        type=int, default=4, help="Thread value to keep.")
    parser.add_argument(
        "--exclude-discarded", 
        action="store_true", default=False)
    args = parser.parse_args()
    
    # Replace tool names with display versions
    tool_name_table = None
    if args.tool_name_file:
        tool_name_table = {}
        with open(args.tool_name_file, 'rt') as inp:
            for row in csv.DictReader(inp, delimiter="\t"):
                tool_name_table[row['ProfileName']] = row['DisplayName']
    
    txt_file = "{}.txt".format(args.output)
    
    if not os.path.exists(txt_file):
        with open(args.input, 'rt') as inp, open(txt_file, 'wt') as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(('Read', 'Program', 'MAPQ', 'In Region?'))
            read_name = None
            read1 = None
            read2 = None
            
            def add_entry(row, read, read_dict):
                if read is not None:
                    q = int(row['read{}_quality'.format(read)])
                    if args.exclude_discarded and q1 < 0:
                        return True
                    else:
                        entry = (
                            "{}.{}".format(read_name, read), prog, q,
                            bool(int(row['read{}_in_region'.format(read)])))
                        if prog in read_dict:
                            assert read_dict[prog] == entry
                        else:
                            read_dict[prog] = entry
                
            for row in csv.DictReader(inp, delimiter="\t"):
                # Since neither threads nor compression change the outcome of trimming,
                # we don't need to worry about stratifying by those metric.
                prog = row['prog']
                if not (prog == 'untrimmed' or int(row['threads']) == args.threads):
                    continue
                
                if read_name != row['read_name']:
                    if read_name:
                        if read1:
                            writer.writerows(read1.values())
                        if read2:
                            writer.writerows(read2.values())
                    read_name = row['read_name']
                    read1 = {}
                    read2 = {}
                
                if add_entry(row, 1, read1):
                    read1 = None
                if add_entry(row, 2, read2):
                    read2 = None
    
    with fileopen(txt_file, 'rt') as inp:
        table = pd.read_csv(inp, sep="\t")
    
    table = (table
        .groupby(['Program', 'MAPQ', 'In Region?'])
        .size()
        .reset_index())
    table = table.rename(columns={ 0 : 'Count' })
    
    untrimmed = table[
        (table.Program == 'untrimmed') & 
        (table.MAPQ >= 0) & 
        (table['In Region?'] == True)]
    trimmed = table[
        (table.Program != 'untrimmed') & 
        (table.MAPQ >= 0) & (table['In Region?'] == True)].copy()
    trimmed.Delta = np.NaN
    for MAPQ in trimmed.MAPQ.unique():
            u = untrimmed.loc[untrimmed.MAPQ==MAPQ, 'Count'].values
            for prog in trimmed.Program.unique():
                t = trimmed.loc[
                    (trimmed.Program==prog) & (trimmed.MAPQ==MAPQ), 
                    'Count'].values
                if len(t) > 0:
                    trimmed.loc[
                        (trimmed.Program==prog) & (trimmed.MAPQ==MAPQ) , 
                        'Delta'] = t[0] - u[0]
    
    plot = sb.factorplot(
        x='MAPQ', y='Delta', hue='Program', col='In Region?', 
        data=trimmed, kind="bar", ci=None, sharey=True)
    plot.set_xlabels('Mapping Quality Score (MAPQ) Cutoff')
    plot.set_ylabels('Difference versus Untrimmed Reads')
    svg_file = args.output + ".svg"
    plot.savefig(svg_file)

if __name__ == "__main__":
    main()
