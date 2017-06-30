#!/usr/bin/env python
"""Given rows from a parse job script, generates tables/figures for certain
performance metrics. Each output from the parse job script is expected to have 
one header row and one metrics row. The first 5 columns are expected to be the 
task profile information.
"""
import argparse
from common import fileopen
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
            elif fmt == 'pickle'
                import pickle
                with fileopen(outfile, 'wb') as out:
                    pickle.dump(table, out)
            else:
                fn = getattr(self, "{}_{}".format(name, fmt))
                fn(table, column, outfile)
        
    def mem_tex(self, table, column, outfile, name=None, caption=None):
        texdat = (table.
            groupby(['Dataset', 'Threads', 'Program2']).
            agg({ column : max }))
        texdat = texdat.rename(columns={ 
            'Program2' : 'Program',
            column : 'Memory'
        })
        texdat.columns = texdat.columns.droplevel()
        from mako.template import Template
        table_template = Template(filename=os.path.join(
            self.template_path, "job_memory.tex"))
        with fileopen(outfile, "wt") as o:
            o.write(table_template.render(
                name=name, caption=caption, table=texdat))
    
    def mem_svg(self, table, column, outfile):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sb
        sb.set(style="whitegrid")
        
        table.assign(log10_mem_mb=np.log10(table[column] / 1000000))
        
        threads = svgdat.Threads.unique()
        if len(threads) == 1:
            plot = sb.factorplot(
                x='Program', y='log10_mem_mb', col="Dataset", 
                data=table, kind="bar", ci=None, sharey=True,
                estimator=max)
        else:
            plot = sb.factorplot(
                x='Threads', y='log10_mem_mb', col="Dataset", hue="Program", 
                data=table, kind="bar", ci=None, sharey=True,
                estimator=max)
        plot.set_xlabels('')
        plot.set_ylabels('Log10(Mb)')
        plot.fig.subplots_adjust(wspace=0.35)
        plot.set_titles('')
        plot.set_xticklabels(rotation=90)
        plot = plot.add_legend()
        plot.savefig(outfile)
    
    def _get_table(self, column):
        cols = list(range(5))
        cols.append(self.header.index(column))
        header = [self.header[c] for c in cols]
        rows = [
            [row[c] for c in cols]
            for row in self.rows
        ]
        return pd.from_records(rows, columns=header)

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
                assert header == line
            else:
                rows.append(line)
    
    metrics = Metrics(
        header, rows, args.output, args.formats, os.path.dirname(__file__))
    
    for metric in args.metric:
        name, key = metric.split('=')
        metrics.show(name, column)

if __name__ == '__main__':
    main()
