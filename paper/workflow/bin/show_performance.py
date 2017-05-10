#!/usr/bin/env python
"""Given rows from parse_gtime, generate a figure and a latex table of results.
"""
import argparse
from common import fileopen
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default="-")
parser.add_argument("-o", "--output", default="-")
parser.add_argument(
    "-f", "--formats", choice=('txt', 'tex', 'svg', 'pickle'), nargs='*',
    default=['tex', 'svg'])
args = parser.parse_args()

# read raw data
with fileopen(args.input, "rt") as inp:
    table = pd.read_csv(inp, sep='\t', names=(
        'Program', 'Threads', 'Dataset', 'Quality', 'Duration', 'CPU', 'Memory'))

# write to file if input was from stdin
if 'txt' in args.formats:
    table.to_csv(args.output + ".txt", sep="\t", index=False)

if 'pickle' in args.format:
    import pickle
    pickle_file = args.output + '.pickle'
    pickle.dump(table, pickle_file)

# generate latex table
if 'tex' in args.formats:
    from mako.template import Template
    table_template = Template(filename=os.path.join(
        os.path.dirname(__file__), "table_template.tex"))
    tex_file = args.output + ".tex"
    with fileopen(tex_file, "wt") as o:
        o.write(table_template.render(table.
            groupby(['Threads', 'Program']).
            agg([min, max]).
            sortlevel()))

# generate figure
if 'svg' in args.formats:
    import seaborn as sb
    threads = table.Threads.unique()
    if len(threads) == 1:
        plot = sb.PairGrid(
            table, x_vars=("Program"), y_vars=('Duration', 'Memory', 'CPU'))
    else:
        plot = sb.PairGrid(
            table, x_vars=("Program"), y_vars=('Duration', 'Memory', 'CPU'),
            hue="Program")
    svg_file = args.output + ".svg"
    plot.map(sb.barplot).savefig(svg_file)