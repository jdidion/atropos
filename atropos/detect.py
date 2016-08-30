# Detect adapter sequences directly from reads based on kmer frequency.

from collections import defaultdict
import math
import re
import statistics as stats
import sys
from urllib.request import urlopen
from .seqio import FastqReader
from .util import reverse_complement, sequence_complexity, enumerate_range
from .xopen import open_output, xopen

def main(options):
    known_contaminants = load_known_contaminants_from_url()
    if options.known_contaminant:
        merge_contaminants(
            known_contaminants,
            load_known_contaminants_from_option_strings(options.known_contaminant))
    if options.known_contaminants_file:
        merge_contaminants(
            known_contaminants,
            load_known_contaminants_from_file(options.known_contaminants_file))
    
    infiles = [f for f in (
        options.single_input, options.input1, options.input2, options.interleaved_input
    ) if f is not None]
    
    n_reads = options.max_reads or 10000
    
    with open_output(options.adapter_report) as o:
        print("Detecting adapters and other potential contaminant sequences based on\n"
              "{}-mers in {} reads".format(options.kmer_size, n_reads), file=o)
        for i, f in enumerate(infiles, 1):
            fq = FastqReader(f)
            if options.progress:
                from .progress import create_progress_reader
                fq = create_progress_reader(fq, max_items=n_reads, counter_magnitude="k")
            file_header = "File {}: {}".format(i, f)
            print("", file=o)
            print(file_header, file=o)
            print('-' * len(file_header), file=o)
            contam = detect_contaminants(fq, k=options.kmer_size, n_reads=n_reads,
                known_contaminants=known_contaminants, include_unknown=options.include_unknown)
            summarize_contaminants(contam, o)

def detect_contaminants(fq, k=12, n_reads=10000, overrep_cutoff=100, known_contaminants=None, include_unknown=True):
    #if known_contaminants and not include_unknown:
    #    contam = detect_known_contaminants(fq, n_reads, overrep_cutoff, known_contaminants)
    #elif n_reads <= 50000:
    #    contam = detect_contaminants_heuristic(fq, k, n_reads, overrep_cutoff, known_contaminants)
    #else:
    contam = detect_contaminants_khmer(fq, k, n_reads, overrep_cutoff, known_contaminants)
    return filter_contaminants(contam, not include_unknown)

def detect_known_contaminants(fq, n_reads, overrep_cutoff, known_contaminants):
    """Test known contaminants against reads."""
    pass

def detect_contaminants_heuristic(fq, k, n_reads, overrep_cutoff, known_contaminants):
    """Use a heuristic iterative algorithm to arrive at likely contaminants."""
    polyA = re.compile('A{8,}.*|A{2,}$')

    def add_to_tree(seq, tree, k):
        for i in range(len(seq) - k + 1):
            kmer = seq[i:(i+k)]
            tree[kmer][0] += 1
            tree[kmer][1].add(seq)

    tree = defaultdict(lambda: [0, set()])
    k = 12
    min_count = 100
    nreads = 20000
    readlen = 125
    overrep_cutoff = 5000

    def expected(k, overrep_cutoff=1):
        num_win = readlen - k + 1
        stat_exp = math.ceil(nreads * num_win * overrep_cutoff / float(4**k))
        return max(0.001 * nreads, stat_exp)

    min_count = math.ceil(expected(k, overrep_cutoff))

    with FastqReader(sys.argv[1]) as fq:
        for i, read in tqdm.tqdm(enumerate(fq)):
            if i == nreads:
                break
            seq = read.sequence
            if sequence_complexity(seq) > 1.0:
                match = polyA.search(seq)
                if match:
                    seq = seq[:match.start()]
                if len(seq) >= k:
                    add_to_tree(seq, tree, k)

    prev = None
    cur = {}
    results = {}

    while True:
        all_seqs = set()
        for kmer, (count, seqs) in tree.items():
            if count > min_count:
                cur[kmer] = (count, seqs)
                all_seqs.update(seqs)
        print("{} {} {}".format(k, min_count, len(all_seqs)))
        if len(all_seqs) == 0:
            break
        if prev:
            for kmer, (count, seqs) in prev.items():
                if not any(seq in cur for seq in seqs) and sequence_complexity(kmer) > 1.0:
                    results[kmer] = count
        prev = cur
        cur = {}
        tree = defaultdict(lambda: [0, set()])
        k += 1
        min_count = math.ceil(expected(k, overrep_cutoff))
        for seq in all_seqs:
            add_to_tree(seq, tree, k)

    results = list(results.items())

    # First merge, adjacent sequences
    results.sort(key=lambda i: i[0])
    cur = results[0]
    merged = []
    for i in range(1, len(results)):
        if results[i][0].startswith(cur[0]):
            cur = [results[i][0], results[i][1] + cur[1]]
        else:
            merged.append(cur)
            cur = results[i]
    merged.append(cur)
    results = merged

    def align(seq1, seq2, min_overlap_frac=0.9):
        aligner = Aligner(
            seq1, 0.0,
            SEMIGLOBAL,
            False, False)
        aligner.min_overlap = math.ceil(min(len(seq1), len(seq2)) * min_overlap_frac)
        aligner.indel_cost = 100000
        match = aligner.locate(seq2)
        if match:
            return seq1[match[0]:match[1]]
        else:
            return None

    # Second merge, by frequency
    results.sort(key=lambda i: i[1], reverse=True)
    cur = results[0]
    merged = []
    unmerged = []
    while len(results) > 1:
        seq1, count1 = results[0]
        for j in range(1, len(results)):
            seq2, count2 = results[j]
            if len(seq1) >= len(seq2) and seq2 in seq1:
                count1 += count2
            elif seq1 in seq2:
                # if they are close in count, keep the longer sequence
                if count1 < (2 * count2):
                    seq1 = seq2
                count1 += count2
            else:
                unmerged.append(results[j])
        merged.append([seq1, count1])
        results = unmerged
        unmerged = []
    merged = merged + results
    merged.sort(key=lambda i: i[1], reverse=True)

    # keep anything that's within 50% of the top hit
    if len(merged) == 0:
        print("No hits :(")
    else:
        min_count = int(merged[0][1] * 0.5)
        merged = filter(lambda x: x[1] >= min_count, merged)

    # match to known sequences
    # TODO

def detect_contaminants_khmer(fq, k, n_reads, overrep_cutoff, known_contaminants):
    import khmer
    from khmer import khmer_args
    
    read = next(fq)
    n_win = len(read.sequence) - k + 1 # assuming all sequences are same length
    tablesize = n_reads * n_win
    countgraph = khmer.Countgraph(k, tablesize, khmer_args.DEFAULT_N_TABLES)
    countgraph.set_use_bigcount(True)
    
    countgraph.consume_and_tag(read.sequence)
    for i, read in enumerate(fq):
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
        
        for seq, names in known_contaminants.items():
            l = len(seq)
            if l < k:
                print("Cannot check {}; sequence is shorter than {}".format(list(names)[0], k))
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
                contam.append((seq, overall_count / float(tablesize), names, float(matches) / n_kmers))
        
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
            print("    Name(s): {}".format(",\n             ".join(row[2])), file=outstream)
            print("    K-mers that are frequent: {:.2%}".format(row[3]), file=outstream)
        print("    Median frequency of k-mers: {:.2%}".format(row[1]), file=outstream)

def load_known_contaminants_from_file(path):
    with open(path, "rt") as i:
        return parse_known_contaminants(i)

def load_known_contaminants_from_option_strings(opt_strings):
    """Parse contaminants from list of name=seq options supplied on command line."""
    return parse_known_contaminants(opt_strings, delim='=')

def load_known_contaminants_from_url(url="https://gist.githubusercontent.com/jdidion/ba7a83c0934abe4bd040d1bfc5752d5f/raw/a6372f21281705ac9031697fcaed3d1f64cea9a5/sequencing_adapters.txt"):
     return parse_known_contaminants(urlopen(url).read().decode().split("\n"))

def parse_known_contaminants(line_iter, delim='\t'):
    regex = re.compile("([^{0}]+){0}+(.+)".format(delim))
    contam = defaultdict(lambda: set())
    for line in line_iter:
        line = line.rstrip()
        if len(line) == 0 or line.startswith('#'):
            continue
        m = regex.match(line)
        if m:
            contam[m.group(2)].add(m.group(1))
    return contam

def merge_contaminants(c1, c2):
    for k, v in c2.items():
        if k in c1:
            c1[k] = c1[k] + c2[k]
        else:
            c1[k] = c2[k]
