#!/usr/bin/env python
"""Given rows from parse_gtime, generate a figure and a latex table of results.
"""
import sys
import numpy as np
import pandas as pd
import seaborn as sb
from mako.template import Template

table_template = Template(r"""
<%
threads = table.Threads.unique()
progs = table.Program.unique()
%>
\begin{table}[ht]
\centering
\begin{tabular}{l${'r|r' * len(threads)}}
\sisetup{detect-weight=true,detect-inline-weight=math}
    % if len(threads) == 1:
    Program & Min & Max \\
    \hfill{} & \multicolumn{2}{c}{Execution Time (sec.)} \\\hline
    % else:
    \hfill{} & ${" & ".join("\multicolumn{{2}}{{c}}{{{} Threads}}".format(t) for t in threads)} \\\hline
    Program & \multicolumn{${len(threads)*2}}{c}{Execution Time (Min$|$Max sec.)} \\\hline
    % endif
    % for p in progs:
        % for t in threads:
        <% row = table.loc[t][p] %>
        ${" & ".join(str(r) for r in row)} 
        
        % endfor 
        \\
    % endfor
    \hline
    
    % if len(threads) == 1:
    \hfill{} & \multicolumn{2}{c}{Memory Usage (GB)} \\\hline
    % else
    \hfill{} & \multicolumn{${len(threads)*2}}{c}{Memory Usage (Min$|$Max GB)} \\\hline
    % endif
    % for p in progs:
        ${p}
        <% row = table.loc[p], 'Memory'] %>
        ${" & ".join(r for r in row)} 
        
        % endfor 
        \\
    % endfor
    \hline
    
    % if len(threads) == 1:
    \hfill{} & \multicolumn{2}{c}{CPU Usage (GB)} \\\hline
    % else
    \hfill{} & \multicolumn{${len(threads)*2}}{c}{CPU Usage (Min$|$Max %)} \\\hline
    % endif
    % for p in progs:
        ${p}
        <% row = table.loc[p], 'CPU'] %>
        ${" & ".join(r for r in row)} 
        
        % endfor 
        \\
    % endfor
    
\end{tabular}
\caption{\label{tab:${name}}${caption}}
\end{table}
""")

# master table
table = pd.read_csv(sys.stdin, sep='\t', names=(
    'Program', 'Threads', 'Dataset', 'Quality', 'Duration', 'CPU', 'Memory'))

# generate table
with open("${output.1}", "wt") as o:
    o.write(table_template.render(table.
        groupby(['Threads', 'Program']).
        agg([min, max]).
        sortlevel()))

# generate figure
threads = table.Threads.unique()
if len(threads) == 1:
    plot = sb.PairGrid(
        table, x_vars=("Program"), y_vars=('Duration', 'Memory', 'CPU'))
else:
    plot = sb.PairGrid(
        table, x_vars=("Program"), y_vars=('Duration', 'Memory', 'CPU'),
        hue="Program")
plot.map(sb.barplot).savefig("{}.svg".format())