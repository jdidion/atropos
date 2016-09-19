# Detect adapter sequences directly from reads based on kmer frequency.

from collections import defaultdict
import logging
import math
import os
import pickle
import re
import statistics as stats
import sys
from urllib.request import urlopen
from .align import Aligner, SEMIGLOBAL
from .seqio import open_reader
from .util import reverse_complement, sequence_complexity, enumerate_range
from .xopen import open_output, xopen

# TODO: Test whether using rc=True in parse_known_contaminants is as fast
# and as accurate as testing both the forward and reverse complement
# read sequence against only the forward-orientation known contaminants.
# Also, offer an option of whether to test the reverse complement, with
# the default being false.

def load_known_contaminants(options):
    cache_file = ".contaminants"
    if not options.no_cache_contaminants and os.path.exists(cache_file):
        with open(cache_file, "rb") as cache:
            return pickle.load(cache)
    else:
        known_contaminants = load_known_contaminants_from_url()
        if options.known_contaminant:
            merge_contaminants(
                known_contaminants,
                load_known_contaminants_from_option_strings(options.known_contaminant))
        if options.known_contaminants_file:
            merge_contaminants(
                known_contaminants,
                load_known_contaminants_from_file(options.known_contaminants_file))
        if not options.no_cache_contaminants:
            # need to make it pickleable
            temp = {}
            for seq, names in known_contaminants.items():
                temp[seq] = list(names)
            with open(cache_file, "wb") as cache:
                pickle.dump(temp, cache)
        return known_contaminants

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
                contam[reverse_complement(seq)].add("{}_rc".format(name))
    return contam

def merge_contaminants(c1, c2):
    for k, v in c2.items():
        if k in c1:
            c1[k] = c1[k] + c2[k]
        else:
            c1[k] = c2[k]

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
        self.abundance = 0
    
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
    
    def match_to(self, seq):
        """
        Determine whether this match's sequence is within 'seq'
        by simple exact string comparison.
        """
        if self.seq in seq:
            self.abundance += 1
            return True
        return False

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
        fw_kmers = set(seq[i:(i+self.k)] for i in range(len(seq) - self.k + 1))
        fw_matches = float(len(self.kmers & fw_kmers))
        
        seqrc = reverse_complement(match.seq)
        rv_kmers = set(seqrc[i:(i+self.k)] for i in range(len(seqrc) - self.k + 1))
        rv_matches = float(len(self.kmers & rv_kmers))
        
        if fw_matches >= rv_matches:
            n_matches = fw_matches
            kmers = fw_kmers
            compare_seq = seq
        else:
            n_matches = rv_matches
            kmers = rv_kmers
            compare_seq = seqrc
        
        self.matches += n_matches
        return (n_matches / self.n_kmers, n_matches / len(kmers), compare_seq)

def create_contaminant_matchers(contaminants, k):
    return [
        ContaminantMatcher(seq, names, k)
        for seq, names in contaminants.items()
    ]

class Detector(object):
    def __init__(self, k=12, n_reads=10000, read_length=None, overrep_cutoff=100,
                 known_contaminants=None, estimate_abundance=True):
        self.k = k
        self.n_reads = n_reads
        self.read_length = read_length
        self.overrep_cutoff = overrep_cutoff
        self.known_contaminants = known_contaminants
        self.estimate_abundance = estimate_abundance
    
    def consume_all(self, reader):
        """
        Consume up to self.n_reads reads from the reader.
        """
        read = next(reader)
        self.consume_first(read)
        for read in enumerate_range(reader, 1, n_reads):
            self.consume(read)
    
    def consume_first(self, read):
        self.read_length = len(read.sequence)
        self.consume(read)
    
    def consume(self, read):
        """
        Consumes a single read.
        """
        raise NotImplemented()
    
    def filter_and_sort(self, include="all", min_len=None, min_complexity=1.1, min_match_frac=0.1, limit=20):
        """
        Returns list of contaminants.
        """
        if min_len is None:
            min_len = self.k
        
        matches = self.get_matches()
        
        def _filter(match):
            if min_len and len(match) < min_len:
                return False
            if min_freq and match.count < self.min_freq:
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
    
    def estimate_abundance(matches, seqs):
        """
        Add abundance information to matches,
        """
        num_contaminated = 0
        for read_seq in seqs:
            contaminated = False
            for match in matches:
                if match.match_to(read_seq):
                    contaminated = True
            if contaminated:
                num_contaminated += 1
        return num_contaminated

class KnownContaminantDetector(Detector):
    """
    Test known contaminants against reads.
    This has linear complexity and is more specific than the khmer matcher, but
    less specific than the heuristic matcher. It's also less sensitive since
    it does not try to detect unknown contaminants.
    """
    def __init__(self, known_contaminants, min_match_frac=0.5, **kwargs):
        super(KnownContaminantDetector, self).__init__(known_contaminants=known_contaminants, **kwargs)
        self.min_k = min(len(s) for s in known_contaminants.keys())
        self.contaminant_matchers = create_contaminant_matchers(known_contaminants, k)
        self.min_match_frac = min_match_frac
        self.counts = defaultdict(lambda: 0)
    
    @property
    def min_freq(self):
        return 0.1
    
    def consume(self, read):
        seq = read.sequence
        seqrc = reverse_complement(seq)
        for contam in self.contaminant_matchers:
            if contam.match(seq)[0] > self.min_match_frac:
                self.counts[contam] += 1
    
    def get_matches(self):
        num_win = self.read_length - self.min_k + 1
        min_count = math.ceil(self.n_reads * num_win * self.overrep_cutoff / float(4**self.min_k))
        contaminants = filter(lambda x: x[1] >= min_count, self.counts.items())
        return list(Match(c[0], match_frac=float(c[1]) / self.n_reads) for c in contaminants)


class HeuristicDetector(Detector):
    @property
    def min_freq(self):
        return 0.1 * self.n_reads

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
    
    # Now merge overlapping sequences by length and frequency to eliminate
    # redundancy in the set of candidate kmers.
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
            match_frac, compare_seq = contam.match(match.seq)
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
            contam, match_frac = contams[0]
            equiv = [c for c in contams[1:] if c[1] == match_frac]
            if len(equiv) == 0:
                match.set_contaminant(contam, *match_frac)
            else:
                names = set(contam.names)
                seqs = set((contam.seq,))
                for e in equiv:
                    names.update(e[0].names)
                    seqs.add(e[0].seq)
                match.set_known(list(names), list(seqs), *match_frac)
        known.append(match)
    
    return known + unknown

class KhmerDetector
    @property
    def min_freq(self):
        return 0.0001

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

def summarize_contaminants(outstream, matches, n_reads, num_contaminated=None, header=None):
    print("", file=outstream)
    
    if header:
        print(header, file=outstream)
        print('-' * len(header), file=outstream)
    
    print("Detected {} adapters/contaminants:".format(len(matches)), file=outstream)
    
    if num_contaminated:
        print("Full-length contaminant(s) found in {} out of {} reads ({:.1%})".format(
            num_contaminated, n_reads, num_contaminated / n_reads), file=outstream)
    
    pad = len(str(len(matches)))
    for i, match in enumerate(matches):
        print(("{:>" + str(pad) + "}. {}").format(i+1, match.seq), file=outstream)
        if match.is_known:
            print("    Name(s): {}".format(",\n             ".join(match.names)), file=outstream)
            print("    Known sequence(s): {}".format(",\n             ".join(match.known_seqs)), file=outstream)
            print("    Known sequence K-mers that match detected contaminant: {:.2%}".format(match.match_frac), file=outstream)
        if match.abundance:
            print("    Abundance (full-length) in {} reads: {} ({:.1%})".format(n_reads, match.abundance, match.abundance / n_reads))
        if match.match_frac2:
            print("    Detected contaminant kmers that match known sequence: {:.2%}".format(match.match_frac2), file=outstream)
        if match.count_is_frequency:
            print("    Frequency of k-mers: {:.2%}".format(match.count), file=outstream)
        else:
            print("    Number of k-mer matches: {}".format(match.count), file=outstream)
