# Detect adapter sequences directly from reads based on kmer frequency.

from collections import defaultdict
import logging
import math
import re
import statistics as stats
import sys
from urllib.request import urlopen
from .seqio import FastqReader
from .util import reverse_complement, sequence_complexity, enumerate_range
from .xopen import open_output, xopen

# TODO: Test whether using rc=True in parse_known_contaminants is as fast
# and as accurate as testing both the forward and reverse complement
# read sequence against only the forward-orientation known contaminants.
# Also, offer an option of whether to test the reverse complement, with
# the default being false.

# TODO: whether and where to cache list of known contaminants?

def main(options):
    known_contaminants = None
    if options.include_contaminants != 'unknown':
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
        print("\nDetecting adapters and other potential contaminant sequences based on "
              "{}-mers in {} reads".format(options.kmer_size, n_reads), file=o)
        for i, f in enumerate(infiles, 1):
            fq = FastqReader(f)
            fq_iter = iter(fq)
            try:
                if options.progress:
                    from .progress import create_progress_reader
                    fq_iter = create_progress_reader(fq_iter, max_items=n_reads, counter_magnitude="K",
                        values_have_size=False)
                file_header = "File {}: {}".format(i, f)
                print("", file=o)
                print(file_header, file=o)
                print('-' * len(file_header), file=o)
                contam = detect_contaminants(fq_iter, k=options.kmer_size, n_reads=n_reads,
                    known_contaminants=known_contaminants, include=options.include_contaminants)
                summarize_contaminants(contam, o)
            finally:
                fq.close()

def detect_contaminants(fq, k=12, n_reads=10000, overrep_cutoff=100, known_contaminants=None, include="all"):
    if known_contaminants and include == 'known':
        logging.getLogger().debug("Detecting contaminants using the known-only algorithm")
        contam = detect_known_contaminants(fq, n_reads, overrep_cutoff, known_contaminants)
        min_freq = 0.1
    elif n_reads <= 50000:
        logging.getLogger().debug("Detecting contaminants using the heuristic algorithm")
        contam = detect_contaminants_heuristic(fq, k, n_reads, overrep_cutoff, known_contaminants)
        min_freq = 0.1 * n_reads
    else:
        logging.getLogger().debug("Detecting contaminants using the kmer-based algorithm")
        contam = detect_contaminants_khmer(fq, k, n_reads, overrep_cutoff, known_contaminants)
        min_freq = 0.0001
    return filter_and_sort_contaminants(contam, include=include, min_freq=min_freq, min_len=k)

def detect_known_contaminants(fq, n_reads, overrep_cutoff, known_contaminants, min_match_frac=0.5):
    """
    Test known contaminants against reads.
    This has linear complexity and is more specific than the khmer matcher, but
    less specific than the heuristic matcher. It's also less sensitive since
    it does not try to detect unknown contaminants.
    """
    k = min(len(s) for s in known_contaminants.keys())
    contaminants = create_contaminant_matchers(known_contaminants, k)
    
    def ingest_read(read, counts):
        seq = read.sequence
        seqrc = reverse_complement(seq)
        for contam in contaminants:
            if (contam.match(seq)[0] > min_match_frac):
                counts[contam] += 1
    
    counts = defaultdict(lambda: 0)
    read = next(fq)
    readlen = len(read.sequence)
    ingest_read(read, counts)
    for i, read in enumerate_range(fq, 1, n_reads):
        ingest_read(read, counts)
    
    num_win = readlen - k + 1
    min_count = math.ceil(n_reads * num_win * overrep_cutoff / float(4**k))
    contaminants = filter(lambda x: x[1] >= min_count, counts.items())
    return list(Match(c[0], match_frac=float(c[1]) / n_reads) for c in contaminants)

def detect_contaminants_heuristic(fq, k_start, n_reads, overrep_cutoff, known_contaminants,
                                  min_freq=0.001, min_contaminant_match_frac=0.9):
    """
    Use a heuristic iterative algorithm to arrive at likely contaminants.
    This is the most accurate algorithm overall, but it has quadratic complexity
    and becomes too slow/memory-intenstive when n_reads > 50k.
    """
    k = k_start
    kmers = defaultdict(lambda: [0, set()])
    read = next(fq)
    readlen = len(read)
    polyA = re.compile('A{8,}.*|A{2,}$')
    
    def expected(k, overrep_cutoff=1):
        num_win = readlen - k + 1
        stat_exp = math.ceil(n_reads * num_win * overrep_cutoff / float(4**k))
        return max(min_freq * n_reads, stat_exp)
    
    def add_kmers(seq, kmers, k):
        for i in range(len(seq) - k + 1):
            kmer = seq[i:(i+k)]
            kmers[kmer][0] += 1
            kmers[kmer][1].add(seq)
            
    def add_kmers_initial(read, kmers):
        seq = read.sequence
        if sequence_complexity(seq) <= 1.0:
            return
        match = polyA.search(seq)
        if match:
            seq = seq[:match.start()]
        if len(seq) >= k:
            add_kmers(seq, kmers, k)
    
    # Build the initial kmer profile
    add_kmers_initial(read, kmers)
    for i, read in enumerate_range(fq, 1, n_reads):
        add_kmers_initial(read, kmers)

    # Identify candidate kmers for increasing values of k
    
    prev = None
    cur = {}
    results = {}
    min_count = math.ceil(expected(k, overrep_cutoff))

    while True:
        all_seqs = set()
        for kmer, (count, seqs) in kmers.items():
            if count > min_count:
                cur[kmer] = (count, seqs)
                all_seqs.update(seqs)
        
        if len(all_seqs) == 0:
            break
        
        if prev:
            for kmer, (count, seqs) in prev.items():
                if not any(seq in cur for seq in seqs) and sequence_complexity(kmer) > 1.0:
                    results[kmer] = count
        
        k += 1
        kmers = defaultdict(lambda: [0, set()])
        for seq in all_seqs:
            add_kmers(seq, kmers, k)
        
        min_count = math.ceil(expected(k, overrep_cutoff))
        prev = cur
        cur = {}
    
    results = list(results.items())
    
    # Now merge overlapping sequences to eliminate redundancy in the
    # set of candidate kmers.
    
    # First merge by length and frequency. This sorting
    # gives preferences to longer sequences.
    results.sort(key=lambda i: len(i[0]) * math.log(i[1]), reverse=True)
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
    results = merged + results
    
    # Merge adjacent sequences
    # This doesn't seem to help, so commenting out
    # results.sort(key=lambda i: i[0])
    # cur = results[0]
    # merged = []
    # for i in range(1, len(results)):
    #     if results[i][0].startswith(cur[0]):
    #         cur = [results[i][0], results[i][1] + cur[1]]
    #     else:
    #         merged.append(cur)
    #         cur = results[i]
    # merged.append(cur)
    # results = merged
    
    if len(results) == 0:
        return []
        
    # Re-sort by frequency
    results.sort(key=lambda i: i[1], reverse=True)
    # Keep anything that's within 50% of the top hit
    min_count = int(results[0][1] * 0.5)
    results = filter(lambda x: x[1] >= min_count, results)
    # Convert to matches
    results = list(Match(x[0], x[1]) for x in results)
    
    if not known_contaminants:
        return results
    
    # Match to known sequences
    contaminants = create_contaminant_matchers(known_contaminants, k_start)
    known = {}
    unknown = []
    
    for match in results:
        best_matches = {}
        best_match_frac = (min_contaminant_match_frac, 0)
        for contam in contaminants:
            match_frac_fw = contam.match(match.seq)
            seqrc = reverse_complement(match.seq)
            match_frac_rv = contam.match(seqrc)
            if match_frac_fw[0] >= match_frac_rv[0]:
                match_frac = match_frac_fw
                compare_seq = match.seq
            else:
                match_frac = match_frac_rv
                compare_seq = seqrc
            if match_frac[0] < best_match_frac[0]:
                continue
            if (contam.seq in compare_seq or
                    align(compare_seq, contam.seq, min_contaminant_match_frac)):
                if (match_frac[0] > best_match_frac[0] or (
                        match_frac[0] == best_match_frac[0] and
                        match_frac[1] > best_match_frac[1])):
                    best_matches = {}
                    best_match_frac = match_frac
                best_matches[contam] = (match, match_frac)
        
        if best_matches:
            for contam, match in best_matches.items():
                if contam not in known or match[1] > known[contam][1]:
                    known[contam] = match
        else:
            unknown.append(match)
    
    # resolve many-many relationships
    
    matches = defaultdict(lambda: [])
    for contam, (match, match_frac) in known.items():
        matches[match].append((contam, match_frac))
    
    known = []
    for match, contams in matches.items():
        if len(contams) == 1:
            contam, match_frac = contams[0]
            match.set_contaminant(contam, *match_frac)
        else:
            contams.sort(key=lambda x: x[1], reverse=True)
            contam, match_frac = equiv[0]
            equiv = [c for c in contams[1:] if c[1] == match_frac]
            if len(equiv) == 0:
                match.set_contaminant(contam, *match_frac)
            else:
                names = list(contam.names)
                seqs = [contam.seq]
                for e in equiv:
                    names.extend(e[0].names)
                    seqs.append(e[0].seq)
                match.set_known(names, seqs, *match_frac)
        known.append(match)
    
    return known + unknown
        
def detect_contaminants_khmer(fq, k, n_reads, overrep_cutoff, known_contaminants):
    """
    Identify contaminants based on kmer frequency using a fast kmer counting
    approach (as implemented in the khmer library). This approach is fast but
    not as accurate as the other two.
    """
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
                contam.append(Match(
                    seq, overall_count / float(tablesize), names, float(matches) / n_kmers))
        
        # Add remaining tags
        extra = [
            Match(tag, candidates[tag] / float(tablesize))
            for tag in set(candidates.keys()) - seen
        ]
        
        return contam + extra
    else:
        return [Match(tag, count / float(tablesize)) for tag, count in candidates.items()]

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

def filter_and_sort_contaminants(matches, include='all', min_freq=0.01, min_len=0,
                                 min_complexity=1.1, min_match_frac=0.1, limit=20):
    def _filter(match):
        if min_len and len(match) < min_len:
            return False
        if min_freq and match.count < min_freq:
            return False
        if min_complexity and match.seq_complexity < min_complexity:
            return False
        if include == 'known' and not match.is_known:
            return False
        elif include == 'unknown' and match.is_known:
            return False
        if min_match_frac and match.is_known and match.match_frac < min_match_frac:
            return False
        return True
    matches = list(filter(_filter, matches))
    
    matches.sort(key=lambda x: len(x) * math.log(x.count), reverse=True)
    
    if limit is not None:
        matches = matches[:limit]
        
    return matches

def summarize_contaminants(matches, outstream):
    print("Detected {} adapters/contaminants:".format(len(matches)), file=outstream)
    pad = len(str(len(matches)))
    for i, match in enumerate(matches):
        print(("{:>" + str(pad) + "}. {}").format(i+1, match.seq), file=outstream)
        if match.is_known:
            print("    Name(s): {}".format(",\n             ".join(match.names)), file=outstream)
            print("    Known sequence(s): {}".format(",\n             ".join(match.known_seqs)), file=outstream)
            print("    K-mers that match: {:.2%}".format(match.match_frac), file=outstream)
        if match.count_is_frequency:
            print("    Frequency of k-mers: {:.2%}".format(match.count), file=outstream)
        else:
            print("    Number of k-mer matches: {}".format(match.count), file=outstream)

class Match(object):
    def __init__(self, seq_or_contam, count=0, names=None, match_frac=None, match_frac2=None):
        if isinstance(seq_or_contam, ContaminantMatcher):
            self.seq = seq_or_contam.seq
            self.count = seq_or_contam.matches
            self.names = tuple(seq_or_contam.names)
            self.known_seqs = [seq_or_contam.seq]
        else:
            self.seq = seq_or_contam
            self.count = count
            self.names = names
            self.known_seqs = None
        self.match_frac = match_frac
        self.match_frac2 = match_frac2
    
    def __len__(self):
        return len(self.seq)
    
    def __repr__(self):
        if self.is_known:
            return "{} => {} ({}))".format(self.seq, self.names, self.known_seqs)
        else:
            return self.seq
    
    @property
    def seq_complexity(self):
        return sequence_complexity(self.seq)
    
    @property
    def count_is_frequency(self):
        return isinstance(self.count, float)
    
    def set_contaminant(self, contam, match_frac, match_frac2=None):
        self.set_known(contam.names, [contam.seq], match_frac, match_frac2)
    
    def set_known(self, names, seqs, match_frac, match_frac2=None):
        self.names = names
        self.known_seqs = seqs
        self.match_frac = match_frac
        self.match_frac2 = match_frac2
    
    @property
    def is_known(self):
        return self.known_seqs is not None

class ContaminantMatcher(object):
    def __init__(self, seq, names, k):
        self.seq = seq
        self.names = names
        self.kmers = set(seq[i:(i+k)] for i in range(len(seq) - k + 1))
        self.n_kmers = len(self.kmers)
        self.k = k
        self.matches = 0
        
    def match(self, seq):
        """
        Returns (num_matches, num_contam_kmers, num_seq_kmers).
        """
        kmers = set(seq[i:(i+self.k)] for i in range(len(seq) - self.k + 1))
        n_matches = float(len(self.kmers & kmers))
        self.matches += n_matches
        return (n_matches / self.n_kmers, n_matches / len(kmers))

def create_contaminant_matchers(contaminants, k):
    return [
        ContaminantMatcher(seq, names, k)
        for seq, names in contaminants.items()
    ]
        
def load_known_contaminants_from_file(path):
    with open(path, "rt") as i:
        return parse_known_contaminants(i)

def load_known_contaminants_from_option_strings(opt_strings):
    """Parse contaminants from list of name=seq options supplied on command line."""
    return parse_known_contaminants(opt_strings, delim='=')

def load_known_contaminants_from_url(url="https://gist.githubusercontent.com/jdidion/ba7a83c0934abe4bd040d1bfc5752d5f/raw/a6372f21281705ac9031697fcaed3d1f64cea9a5/sequencing_adapters.txt"):
    logging.getLogger().info("\nDownloading list of known contaminants from {}".format(url))
    return parse_known_contaminants(urlopen(url).read().decode().split("\n"))

def parse_known_contaminants(line_iter, delim='\t', rc=False):
    regex = re.compile("([^{0}]+){0}+(.+)".format(delim))
    contam = defaultdict(lambda: set())
    for line in line_iter:
        line = line.rstrip()
        if len(line) == 0 or line.startswith('#'):
            continue
        m = regex.match(line)
        if m:
            seq = m.group(2)
            name = m.group(1)
            contam[seq].add(name)
            if rc:
                contam[reverse_complement(seq)] = "{}_rc".format(name)
    return contam

def merge_contaminants(c1, c2):
    for k, v in c2.items():
        if k in c1:
            c1[k] = c1[k] + c2[k]
        else:
            c1[k] = c2[k]
