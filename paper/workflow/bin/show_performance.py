#!/usr/bin/env python
"""Given rows from parse_gtime, generate a figure and a latex table of results.
"""
import argparse
from common import fileopen
import os
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("-n", "--name", default="table")
    parser.add_argument("-c", "--caption", default="")
    parser.add_argument(
        "-t", "--threads", type=int, default=None,
        help="Set all rows to have the same value for the Threads column.")
    parser.add_argument(
        "-f", "--formats", choices=('txt', 'tex', 'svg', 'pickle'), nargs='*',
        default=['tex', 'svg'])
    args = parser.parse_args()

    # read raw data
    with fileopen(args.input, "rt") as inp:
        table = pd.read_csv(
            inp, sep='\t', names=(
                'Program', 'Program2', 'Threads', 'Dataset', 'Quality', 
                'DurationSecs', 'CPUPct', 'MemoryMB'),
            dtype={ 'Program' : 'category', 'Dataset' : 'category'})
    
    if args.threads:
        table.Threads = args.threads
    
    # save table (useful if input was stdin)

    if 'txt' in args.formats:
        table.to_csv(args.output + ".txt", sep="\t", index=False)

    if 'pickle' in args.formats:
        import pickle
        pickle_file = args.output + '.pickle'
        with fileopen(pickle_file, 'wb') as out:
            pickle.dump(table, out)

    # generate latex table
    if 'tex' in args.formats:
        texdat = table.melt(
            id_vars=['Program2', 'Threads', 'Dataset', 'Quality'], 
            value_vars=['DurationSecs','CPUPct','MemoryMB'])
        texdat = (texdat.
            groupby(['Threads', 'Program2', 'variable']).
            agg({ 'value' : [min, max] }).
            sort_index())
        texdat.columns = texdat.columns.droplevel()
        from mako.template import Template
        table_template = Template(filename=os.path.join(
            os.path.dirname(__file__), "performance_table.tex"))
        tex_file = args.output + ".tex"
        with fileopen(tex_file, "wt") as o:
            o.write(table_template.render(
                name=args.name, caption=args.caption, table=texdat))

    # generate figure
    if 'svg' in args.formats:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sb
        sb.set(style="whitegrid")
        import numpy as np
        
        progs = list(zip(
            ('adapterremoval', 'atropos (adapter)', 'atropos (insert)', 'seqpurge', 'skewer'),
            ('AdapterRemoval', 'Atropos (adapter)', 'Atropos (insert)', 'SeqPurge', 'Skewer'),
            sb.color_palette(n_colors=5)
        ))
        pal = dict((p[1], p[2]) for p in progs)
        
        svgdat = table.melt(
            id_vars=['Program', 'Threads', 'Dataset', 'Quality'], 
            value_vars=['DurationSecs','CPUPct','MemoryMB'])
        svgdat['Program'] = svgdat['Program'].map(dict((p[0], p[1]) for p in progs))
        svgdat['Program'] = svgdat['Program'].astype('category')
        svgdat['Dataset'] = svgdat['Dataset'].astype('category')
        svgdat['variable'] = pd.Categorical(
            svgdat['variable'], categories=['DurationSecs', 'MemoryMB', 'CPUPct'])
        
        threads = svgdat.Threads.unique()
        plot = sb.factorplot(
            x='Threads', y="value", col="variable", hue="Program", 
            data=svgdat, kind="bar", sharey=False, estimator=np.mean)
        if len(threads) == 1:
            plot.set_xticklabels('')
            plot.set_xlabels('')
        else:
            plot.set_xlabels('Threads')
        plot.axes[0,0].set_ylabel('Runtime (sec)')
        plot.axes[0,1].set_ylabel('Memory (Mb)')
        plot.axes[0,2].set_ylabel('CPU (%)')
        plot.fig.subplots_adjust(wspace=0.35)
        plot.set_titles('')
        svg_file = args.output + ".svg"
        plot.savefig(svg_file)

if __name__ == "__main__":
    main()
