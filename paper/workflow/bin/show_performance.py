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
parser.add_argument("-n", "--name", default="table")
parser.add_argument("-c", "--caption", default="")
parser.add_argument(
    "-f", "--formats", choices=('txt', 'tex', 'svg', 'pickle'), nargs='*',
    default=['tex', 'svg'])
args = parser.parse_args()

# read raw data
with fileopen(args.input, "rt") as inp:
    table = pd.read_csv(inp, sep='\t', names=(
        'Program', 'Threads', 'Dataset', 'Quality', 'Duration', 'DurationSecs', 
        'CPUPct', 'MemoryMB'))

# write to file if input was from stdin
if 'txt' in args.formats:
    table.to_csv(args.output + ".txt", sep="\t", index=False)

if 'pickle' in args.formats:
    import pickle
    pickle_file = args.output + '.pickle'
    with fileopen(pickle_file, 'wb') as out:
        pickle.dump(table, out)

# generate latex table
if 'tex' in args.formats:
    from mako.template import Template
    table_template = Template(filename=os.path.join(
        os.path.dirname(__file__), "table_template.tex"))
    tex_file = args.output + ".tex"
    with fileopen(tex_file, "wt") as o:
        o.write(table_template.render(
            name=args.name, caption=args.caption,
            table=table.
                groupby(['Threads', 'Program']).
                agg([min, max]).
                sort_index()))

# generate figure
if 'svg' in args.formats:
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sb
    threads = table.Threads.unique()
    if len(threads) == 1:
        plot = sb.PairGrid(
            table, x_vars=("Program"), y_vars=('DurationSecs', 'MemoryMB', 'CPUPct'))
    else:
        plot = sb.PairGrid(
            table, x_vars=("Program"), y_vars=('DurationSecs', 'MemoryMB', 'CPUPct'),
            hue="Program")
    svg_file = args.output + ".svg"
    plot = plot.map(sb.barplot)
    for i, lab in enumerate(["Duration (sec)", "Memory (MB)", "CPU (%)"]):
        plot.axes[i][0].set_ylabel(lab)
    plot.axes[1][0].set_yscale('log', basey=10)
    plot = plot.add_legend()
    plot.savefig(svg_file)
