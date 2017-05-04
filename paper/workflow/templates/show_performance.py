#!/usr/bin/env python
"""Given rows from parse_gtime, generate a figure and a latex table of results.
"""
import sys
import numpy as np
import pandas as pd
import seaborn as sb
from mako.template import Template

TABLE_TEMPLATE = Template(r"""
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
        <% row = table[lambda df: df.Program == p and df.Threads == t, :] %>
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
        <% row = table[table.loc['Program'] == p, 'Duration'] %>
        ${" & ".join(r for r in row)} 
        
        % endfor 
        \\
    % endfor
    \hline
    
\end{tabular}
\caption{\label{tab:${name}}${caption}}
\end{table}
""")



# master table
table = pd.read_csv(sys.stdin, sep='\t', names=(
    'Program', 'Threads', 'Dataset', 'Quality', 'Duration', 'CPU', 'Memory'))

# generate table
with open("${output.1}", "wt") as o:
    latex = template.render(table.
        groupby(['Program', 'Threads']).
        agg([min, max]).
        sort_values(['Threads', 'Program']))
    o.write(latex)

