#!/usr/bin/env python
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumenetParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "-f", "--formats", choices=('txt', 'svg', 'pickle'), nargs="+",
        default=['tex'])
    parser.add_argument(
        "-t", "--tool-name-file",
        help="File that maps profile names to display names for tools")
    parser.add_argument("--exclude-discarded", action="store_true", default=False)
    args = parser.parse_args()
    
    with fileopen(args.input, 'rt') as inp:
        table = pd.read_csv(inp, sep="\t")
    
    if 'txt' in args.formats:
        with fileopen(args.output + '.txt', 'wt') as out:
            table.to_csv(out, sep="\t", index=False)
        
    if 'pickle' in args.formats:
        import pickle
        with fileopen(args.output + '.pickle', 'wb') as out:
            pickle.dump(table, out, protocol=pickle.HIGHEST_PROTOCOL)
    
    if 'svg' in args.formats:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sb
        sb.set(style="whitegrid")
        
        # Since neither threads nor compression change the outcome of trimming,
        # we don't need to worry about stratifying by those metric.
        table = table[(table.threads==4) | (table.prog == 'untrimmed')]
        table = table.drop(['prog2', 'threads', 'dataset', 'qcut'], 1)
        
        r1 = table[['prog','read1_in_region','read1_quality']]
        r1 = r1.rename(columns={ 
            'read1_in_region' : 'in_region',
            'read1_quality' : 'quality'
        })
        r2 = table[['prog','read2_in_region','read2_quality']]
        r2 = r2.rename(columns={ 
            'read2_in_region' : 'in_region',
            'read2_quality' : 'quality'
        })
        rna_reads = (r1
            .append(r2)
            .groupby(['prog', 'quality', 'in_region'])
            .size()
            .reset_index())
        rna_reads.in_region = rna_reads.in_region.map(bool)
        rna_reads = rna_reads.rename(columns={
            'prog' : 'Program',
            'quality' : 'MAPQ',
            'in_region' : 'In Region?',
            0 : 'Count'
        })
        
        # rna_reads['Delta'] = np.NaN
        # for mapq in set(rna_reads.MAPQ.values):
        #     for inregion in (True, False):
        #         untrimmed_row = ((rna_reads['Program']==prog) & 
        #             (rna_reads['MAPQ']==mapq) & 
        #             (rna_reads['In Region?']==inregion))
        #         for prog in set(rna_reads.Program.values):
        #             prog_row = ((rna_reads['Program']==prog) &
        #                 (rna_reads['MAPQ']==mapq) & 
        #                 (rna_reads['In Region?']==inregion))
        #             a = rna_reads[prog_row]['Count']
        #             b = rna_reads[untrimmed_row]['Count']
        #             b = 0 if b.size == 0 else b.values[0]
        #             if a.size == 0:
        #                 diff = -b
        #                 new_rows.append([
        #                     prog, mapq, inregion, 0, diff
        #                 ])
        #             else:
        #                 diff = a.values[0] - b
        #                 rna_reads.loc[prog_row, 'Delta'] = diff
        
        rna_reads = rna_reads.append(
            pd.DataFrame(new_rows, columns=rna_reads.columns))
        
        # Replace tool names with display versions
        if args.tool_name_file:
            with open(args.tool_name_file, 'rt') as inp:
                tool_name_table = pd.read_csv(inp, sep="\t", index_col='ProfileName')
            rna_reads.Program = rna_reads.Program.map(
                lambda x: tool_name_table.loc[x, 'DisplayName']).values
        
        plot = sb.factorplot(
            x='MAPQ', y='Count', hue='Program', col='In Region?', data)
        plot.set_xlabels('Mapping Quality Score (MAPQ) Cutoff')
        plot.set_ylabels('Number of Aligned Reads')
        svg_file = args.output + ".svg"
        plot.savefig(svg_file)

if __name__ == "__main__":
    main()