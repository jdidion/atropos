# Detect adapter sequences directly from reads based on kmer frequency.

import math
import re
import statistics as stats
import sys
from urllib.request import urlopen
import khmer
from khmer import khmer_args
from .util import reverse_complement, sequence_complexity
from .xopen import open_output

def main(options):
    known_contaminants = load_known_contaminants_from_url()
    infiles = [f for f in (
        options.single_input, options.input1, options.input2, options.interleaved_input
    ) if f is not None]
    with open_output(options.adapter_report) as o:
        print("Detecting adapters and other potential contaminant sequences based on\n"
              "{}-mers in {} reads".format(options.kmer_size, options.max_reads), file=o)
        for i, f in enumerate(infiles, 1):
            file_header = "File {}: {}".format(i, f)
            print("", file=o)
            print(file_header, file=o)
            print('-' * len(file_header), file=o)
            contam = detect_contaminants(f, k=options.kmer_size, n_reads=options.max_reads,
                known_contaminants=known_contaminants)
            contam = filter_contaminants(contam, known_only=True)
            summarize_contaminants(contam, o)

def detect_contaminants(fq, k=12, n_reads=10000, overrep_cutoff=100, known_contaminants=None):
    parser = khmer.ReadParser(fq)
    read = next(parser)
    n_win = len(read.sequence) - k + 1 # assuming all sequences are same length
    tablesize = n_reads * n_win
    
    countgraph = khmer.Countgraph(k, tablesize, khmer_args.DEFAULT_N_TABLES)
    countgraph.set_use_bigcount(True)
    
    countgraph.consume_and_tag(read.sequence)
    for i, read in enumerate(parser):
        if i >= n_reads:
            break
        countgraph.consume_and_tag(read.sequence)
    
    n_expected = math.ceil(tablesize / float(4**k))
    min_count = n_expected * overrep_cutoff
    if min_count >= 2**16:
        raise Exception("The minimum count for an over-represented k-kmer {} "
                        "is greater than the max khmer count (2^16)".format(min_count))
    candidates = {}
    
    for tag in countgraph.get_tagset():
        count = countgraph.get(tag)
        if count >= min_count:
            candidates[tag] = count
    
    if known_contaminants:
        contam = []
        seen = set()
        
        def match(kmer):
            n = candidates.get(kmer, 0)
            if n > 0:
                seen.add(kmer)
            return n
        
        for name, seq in known_contaminants.items():
            l = len(seq)
            if l < k:
                print("Cannot check {}; sequence is shorter than {}".format(name, k))
                continue
            
            n_kmers = l - k + 1
            matches = 0
            match_counts = []
            for i in range(n_kmers):
                kmer = seq[i:(i+k)]
                kmer_count = max(
                    match(kmer),
                    match(reverse_complement(kmer))
                )
                if kmer_count > 0:
                    matches += 1
                    match_counts.append(kmer_count)
            
            if matches > 0:
                # not sure what the correct metric is to use here
                overall_count = sum(match_counts) / float(n_kmers)
                contam.append((seq, overall_count / float(tablesize), name, float(matches) / n_kmers))
        
        contam.sort(key=lambda x: x[3], reverse=True)
        
        # Add remaining tags
        extra = [(tag, candidates[tag] / float(tablesize)) for tag in set(candidates.keys()) - seen]
        extra.sort(key=lambda x: x[1], reverse=True)
        
        return contam + extra
    else:
        contam = [(tag, count / float(tablesize)) for tag, count in candidates.items()]
        contam.sort(key=lambda x: x[1], reverse=True)
        return contam

def filter_contaminants(contam, min_complexity=1.1, min_len=None, min_freq=0.0001, known_only=False, limit=20):
    def _filter(row):
        if min_len is not None and len(row[0]) < min_len:
            return False
        if min_freq is not None and row[1] < min_freq:
            return False
        if min_complexity is not None and sequence_complexity(row[0]) < min_complexity:
            return False
        if known_only and len(row) == 2:
            return False
        return True
    contam = list(filter(_filter, contam))
    if limit is not None:
        contam = contam[:limit]
    return contam

def summarize_contaminants(contam, outstream):
    print("Detected {} adapters/contaminants:".format(len(contam)), file=outstream)
    pad = len(str(len(contam)))
    for i, row in enumerate(contam):
        print(("{:>" + str(pad) + "}. {}").format(i+1, row[0]), file=outstream)
        if len(row) > 2:
            print("    Name: {}".format(row[2]), file=outstream)
            print("    K-mers that are frequent: {:.2%}".format(row[3]), file=outstream)
        print("    Median frequency of k-mers: {:.2%}".format(row[1]), file=outstream)

def load_known_contaminants(path):
    with open(path, "rt") as i:
        return parse_known_contaminants(i)

def load_known_contaminants_from_url(url="https://gist.githubusercontent.com/jdidion/ba7a83c0934abe4bd040d1bfc5752d5f/raw/a6372f21281705ac9031697fcaed3d1f64cea9a5/sequencing_adapters.txt"):
     return parse_known_contaminants(urlopen(url).read().decode().split("\n"))

def parse_known_contaminants(line_iter):
    regex = re.compile('([^\t]+)\t+(.+)')
    contam = {}
    for line in line_iter:
        line = line.rstrip()
        if len(line) == 0 or line.startswith('#'):
            continue
        m = regex.match(line)
        if m:
            contam[m.group(1)] = m.group(2)
    return contam
