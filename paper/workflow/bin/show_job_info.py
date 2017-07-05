#!/usr/bin/env python
"""Given rows from a parse job script, generates tables/figures for certain
performance metrics. Each output from the parse job script is expected to have 
one header row and one metrics row. The first 5 columns are expected to be the 
task profile information.
"""
import argparse
from common import fileopen, parse_size
import csv
import os
import pandas as pd

class Metrics():
    def __init__(self, header, rows, output, formats, template_path):
        self.header = header
        self.rows = rows
        self.output = output
        self.formats = formats
        self.template_path = template_path
    
    def show(self, name, column, **kwargs):
        table = self._get_table(column)
        prefix = "{}.{}".format(self.output, name)
        
        for fmt in self.formats:
            outfile = "{}.{}".format(prefix, fmt)
            if fmt == 'txt':
                table.to_csv(outfile, sep="\t", index=False)
            elif fmt == 'pickle':
                import pickle
                with fileopen(outfile, 'wb') as out:
                    pickle.dump(table, out)
            else:
                fn = getattr(self, "{}_{}".format(name, fmt))
                fn(table, column, outfile)
        
    def mem_tex(self, table, column, outfile, name=None, caption=None):
        texdat = (table.
            drop('Program', 1).
            rename(columns={
                'Program2' : 'Program',
                column : 'Memory'
            }).
            groupby(['Dataset', 'Threads', 'Program']).
            agg({ 'Memory' : max }))
        texdat = texdat.assign(MemoryMB=round(texdat['Memory'] / 1000000, 1))
        
        from mako.template import Template
        table_template = Template(filename=os.path.join(
            self.template_path, "job_memory_table.tex"))
        with fileopen(outfile, "wt") as o:
            o.write(table_template.render(
                name=name, caption=caption, table=texdat))
    
    def mem_svg(self, table, column, outfile):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sb
        sb.set(style="whitegrid")
        
        svgdat = (table.
            rename(columns={ column : 'Memory' }).
            groupby(['Dataset', 'Threads', 'Program']).
            agg({ 'Memory' : max }).
            reset_index())
        svgdat = svgdat.assign(MemoryMB=svgdat['Memory'] / 1000000)
        
        threads = svgdat.Threads.unique()
        if len(threads) == 1:
            plot = sb.factorplot(
                x='Program', y='MemoryMB', col="Dataset", 
                data=svgdat, kind="bar", ci=None, sharey=True)
        else:
            plot = sb.factorplot(
                x='Threads', y='MemoryMB', col="Dataset", hue="Program", 
                data=svgdat, kind="bar", ci=None, sharey=True)

        if len(threads) == 1:
            plot = plot.set_titles('')

        plot = plot.set_xlabels('Threads')
        plot = plot.set_ylabels('Memory (MB)')
        plot = plot.set_xticklabels(rotation=90)
        plot.fig.subplots_adjust(wspace=0.35)
        plot.savefig(outfile)
    
    def _get_table(self, column, is_size=True):
        cols = list(range(5))
        cols.append(self.header.index(column))
        header = [self.header[c] for c in cols]
        rows = [
            [row[c] for c in cols]
            for row in self.rows
        ]
        if is_size:
            for row in rows:
                row[5] = parse_size(row[5])
        table = pd.DataFrame.from_records(rows, columns=header)
        table = table.rename(columns={ 
            'prog' : 'Program',
            'prog2' : 'Program2',
            'threads' : 'Threads',
            'dataset' : 'Dataset',
            'qcut' : 'Quality',
        })
        table['Threads'] = pd.to_numeric(table['Threads'])
        table['Dataset'] = pd.Categorical(table['Dataset'])
        table['Program'] = pd.Categorical(table['Program'])
        table['Program2'] = pd.Categorical(table['Program2'])
        return table

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("-m", "--metric", nargs='*', default=None)
    parser.add_argument(
        "-f", "--formats", 
        choices=('txt', 'tex', 'pickle', 'svg'), nargs='*', default=['txt'])
    args = parser.parse_args()
    
    header = None
    rows = []
    
    with fileopen(args.input, 'rt') as inp:
        for i, line in enumerate(csv.reader(inp, delimiter="\t")):
            if i == 0:
                header = line
            elif i % 2 == 0:
                if len(line) == 0:
                    break
                else:
                    assert header == line
            else:
                rows.append(line)
    
    metrics = Metrics(
        header, rows, args.output, args.formats, os.path.dirname(__file__))
    
    for metric in args.metric:
        name, column = metric.split('=')
        metrics.show(name, column)

if __name__ == '__main__':
    main()
