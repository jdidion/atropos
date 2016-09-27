# Detect adapter sequences directly from reads based on kmer frequency.

from collections import defaultdict
import logging
import math
import os
import re
import statistics as stats
import sys
from .align import Aligner, SEMIGLOBAL
from .seqio import open_reader
from .util import reverse_complement, sequence_complexity, enumerate_range
from .xopen import open_output, xopen

# TODO: Test whether using rc=True in parse_known_contaminants is as fast
# and as accurate as testing both the forward and reverse complement
# read sequence against only the forward-orientation known contaminants.
# Also, offer an option of whether to test the reverse complement, with
# the default being false.

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
        self.abundance = None
    
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
    
    def estimate_abundance(self, read_sequences):
        """
        Determine whether this match's sequence is within 'seq'
        by simple exact string comparison.
        """
        self.abundance = sum(1 for read_seq in read_sequences if self.seq in read_seq)

class ContaminantMatcher(object):
    def __init__(self, seq, names, k):
        self.seq = seq
        self.names = names
        self.kmers = set(seq[i:(i+k)] for i in range(len(seq) - k + 1))
        self.n_kmers = len(self.kmers)
        self.k = k
        self.matches = 0
        
    def match(self, seq, seqrc):
        """
        Returns (num_matches, num_contam_kmers, num_seq_kmers).
        """
        fw_kmers = set(seq[i:(i+self.k)] for i in range(len(seq) - self.k + 1))
        fw_matches = float(len(self.kmers & fw_kmers))
        
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
        for seq, names in contaminants.iter_sequences()
    ]

POLY_A = re.compile('A{8,}.*|A{2,}$')

class Detector(object):
    def __init__(self, k=12, n_reads=10000, overrep_cutoff=100, known_contaminants=None):
        self.k = k
        self.n_reads = n_reads
        self.overrep_cutoff = overrep_cutoff
        self.known_contaminants = known_contaminants
        self._read_length = None
        self._read_sequences = set()
        self._matches = None
    
    def consume_all(self, reader):
        """
        Consume up to self.n_reads reads from the reader.
        """
        read = next(reader)
        self.consume_first(read)
        for read in enumerate_range(reader, 1, n_reads):
            self.consume(read)
    
    def consume_all_batches(self, batch_iterator):
        """
        Consume all reads from the specified batch_iterator.
        It is expected that the iterator was constructed
        with max_reads == n_reads.
        """
        for batch_num, (batch_size, batch) in enumerate(batch_iterator):
            if batch_num == 0:
                self.consume_first(batch[0])
                batch = batch[1:]
            for read in batch:
                self.consume(read)
        
    def consume_first(self, read):
        assert self._read_length is None
        self._read_length = len(read.sequence)
        self.consume(read)
    
    def consume(self, read):
        """
        Consumes a single read.
        """
        seq = self._filter_seq(read.sequence)
        if seq:
            self._read_sequences.add(seq)
    
    def _filter_seq(self, seq):
        if sequence_complexity(seq) <= 1.0:
            return None
        match = POLY_A.search(seq)
        if match:
            seq = seq[:match.start()]
        if len(seq) < self.k:
            return None
        return seq
    
    def matches(self, **kwargs):
        if self._matches is None or len(kwargs) > 0:
            self._filter_and_sort(**kwargs)
        return self._matches
    
    def _filter_and_sort(self, include="all", min_len=None, min_complexity=1.1, min_match_frac=0.1, limit=20):
        """
        Identify, filter, and sort contaminants.
        """
        if min_len is None:
            min_len = self.k
        
        matches = self._get_contaminants()
        
        for match in matches:
            match.estimate_abundance(self._read_sequences)
        
        def _filter(match):
            if match.count < self.min_report_freq:
                return False
            if min_len and len(match) < min_len:
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
        
        self._matches = matches
    
    def summarize(self, outstream, name=None, **kwargs):
        header = "File: {}".format(name) if name else None
        summarize_contaminants(outstream, self.matches(**kwargs), self.n_reads, header)

class PairedDetector(object):
    def __init__(self, detector_class, **kwargs):
        self.d1 = detector_class(**kwargs)
        self.d2 = detector_class(**kwargs)
    
    def consume_all(self, reader):
        read1, read2 = next(reader)
        self.d1.consume_first(read1)
        self.d2.consume_first(read2)
        for read1, read2 in reader:
            self.d1.consume(read1)
            self.d2.consume(read2)
    
    def consume_all_batches(self, batch_iterator):
        for batch_num, (batch_size, batch) in enumerate(batch_iterator):
            if batch_num == 0:
                read1, read2 = batch[0]
                self.d1.consume_first(read1)
                self.d2.consume_first(read2)
                batch = batch[1:]
            for read1, read2 in batch:
                self.d1.consume(read1)
                self.d2.consume(read2)
    
    def matches(self, **kwargs):
        return (
            self.d1.matches(**kwargs),
            self.d2.matches(**kwargs))
    
    def summarize(self, outstream, names=[None, None], **kwargs):
        name1, name2 = names
        self.d1.summarize(outstream, name1, **kwargs)
        self.d2.summarize(outstream, name2, **kwargs)

class KnownContaminantDetector(Detector):
    """
    Test known contaminants against reads.
    This has linear complexity and is more specific than the khmer matcher, but
    less specific than the heuristic matcher. It's also less sensitive since
    it does not try to detect unknown contaminants.
    """
    def __init__(self, known_contaminants, min_match_frac=0.5, **kwargs):
        super(KnownContaminantDetector, self).__init__(known_contaminants=known_contaminants, **kwargs)
        self.min_match_frac = min_match_frac
        self._min_k = min(len(s) for s in known_contaminants.sequences)
    
    @property
    def min_report_freq(self):
        return 0.1
    
    def _filter_seq(self, seq):
        seq = super(KnownContaminantDetector, self)._filter_seq(seq)
        if seq and len(seq) >= self._min_k:
            return seq
        return None
    
    def _get_contaminants(self, read_seqs):
        contaminant_matchers = create_contaminant_matchers(known_contaminants, k)
        counts = defaultdict(lambda: 0)

        for seq in read_seqs:
            seqrc = reverse_complement(seq)
            for contam in contaminant_matchers:
                match = contam.match(seq, seqrc)
                if match[0] > self.min_match_frac:
                    counts[contam] += 1
        
        min_count = math.ceil(
            self.n_reads * (self._read_length - self._min_k + 1) * self.overrep_cutoff /
            float(4**self._min_k))
        
        return [
            Match(c[0], match_frac=float(c[1]) / self.n_reads)
            for c in filter(
                lambda x: x[1] >= min_count,
                counts.items()
            )
        ]

class HeuristicDetector(Detector):
    """
    Use a heuristic iterative algorithm to arrive at likely contaminants.
    This is the most accurate algorithm overall, but it has quadratic complexity
    and becomes too slow/memory-intenstive when n_reads > 50k.
    """
    def __init__(self, min_freq=0.001, min_contaminant_match_frac=0.9, **kwargs):
        super(HeuristicDetector, self).__init__(**kwargs)
        self.min_freq = min_freq
        self.min_contaminant_match_frac = min_contaminant_match_frac
    
    @property
    def min_report_freq(self):
        return 0.1 * self.n_reads
    
    def _get_contaminants(self):
        def _min_count(k):
            return math.ceil(self.n_reads * max(
                self.min_freq,
                (self._read_length - k + 1) * self.overrep_cutoff / float(4**k)))
        
        k = self.k
        kmers = defaultdict(lambda: [0, set()])
        
        for seq in self._read_sequences:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:(i+k)]
                kmers[kmer][0] += 1
                kmers[kmer][1].add(seq)
        
        prev = None
        cur = {}
        results = {}
        min_count = _min_count(k)
        
        # Identify candidate kmers for increasing values of k
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
                for i in range(len(seq) - k + 1):
                    kmer = seq[i:(i+k)]
                    kmers[kmer][0] += 1
                    kmers[kmer][1].add(seq)
            
            min_count = _min_count(k)
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
        
        if len(results) == 0:
            return []
            
        # Re-sort by frequency
        results.sort(key=lambda i: i[1], reverse=True)
        # Keep anything that's within 50% of the top hit
        # TODO: make this user-configurable?
        min_count = int(results[0][1] * 0.5)
        results = filter(lambda x: x[1] >= min_count, results)
        # Convert to matches
        matches = [Match(x[0], x[1]) for x in results]
        
        if self.known_contaminants:
            # Match to known sequences
            contaminants = create_contaminant_matchers(self.known_contaminants, self.k)
            known = {}
            unknown = []
            
            for match in matches:
                seq = match.seq
                seqrc = reverse_complement(seq)
                best_matches = {}
                best_match_frac = (self.min_contaminant_match_frac, 0)
                for contam in contaminants:
                    match_frac1, match_frac2, compare_seq = contam.match(seq, seqrc)
                    if match_frac1 < best_match_frac[0]:
                        continue
                    if (contam.seq in compare_seq or
                            align(compare_seq, contam.seq, self.min_contaminant_match_frac)):
                        if (match_frac1 > best_match_frac[0] or (
                                match_frac1 == best_match_frac[0] and
                                match_frac2 > best_match_frac[1])):
                            best_matches = {}
                            best_match_frac = (match_frac1, match_frac2)
                        best_matches[contam] = (match, (match_frac1, match_frac2))
                
                if best_matches:
                    for contam, match in best_matches.items():
                        if contam not in known or match[1] > known[contam][1]:
                            known[contam] = match
                else:
                    unknown.append(match)
            
            # resolve many-many relationships
            
            new_matches = defaultdict(lambda: [])
            for contam, (match, match_frac) in known.items():
                new_matches[match].append((contam, match_frac))
            
            known = []
            for match, contams in new_matches.items():
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
            
            matches = known + unknown
        
        return matches

class KhmerDetector(Detector):
    """
    Identify contaminants based on kmer frequency using a fast kmer counting
    approach (as implemented in the khmer library). This approach is fast but
    not as accurate as the other two.
    """
    def __init__(self, **kwargs):
        # make sure khmer is installed
        import khmer
        super(KhmerDetector, self).__init__(**kwargs)
        
    @property
    def min_report_freq(self):
        return 0.0001
    
    def _get_contaminants(self):
        from khmer import khmer_args
        n_win = self._read_length - self.k + 1 # assuming all sequences are same length
        tablesize = n_reads * n_win
        countgraph = khmer.Countgraph(self.k, tablesize, khmer_args.DEFAULT_N_TABLES)
        countgraph.set_use_bigcount(True)
        
        for seq in self._read_sequences:
            countgraph.consume_and_tag(seq)
        
        n_expected = math.ceil(tablesize / float(4**k))
        min_count = n_expected * self.overrep_cutoff
        if min_count >= 2**16:
            raise Exception("The minimum count for an over-represented k-kmer {} "
                            "is greater than the max khmer count (2^16)".format(min_count))
    
        candidates = {}
        
        for tag in countgraph.get_tagset():
            count = countgraph.get(tag)
            if count >= min_count:
                candidates[tag] = count
        
        if self.known_contaminants:
            matches = []
            seen = set()
            
            def match(kmer):
                n = candidates.get(kmer, 0)
                if n > 0:
                    seen.add(kmer)
                return n
            
            for seq, names in self.known_contaminants.iter_sequences():
                l = len(seq)
                if l < k:
                    print("Cannot check {}; sequence is shorter than {}".format(list(names)[0], k))
                    continue
                
                n_kmers = l - self.k + 1
                num_matches = 0
                match_counts = []
                for i in range(n_kmers):
                    kmer = seq[i:(i+k)]
                    kmer_count = max(
                        match(kmer),
                        match(reverse_complement(kmer))
                    )
                    if kmer_count > 0:
                        num_matches += 1
                        match_counts.append(kmer_count)
                
                if num_matches > 0:
                    # not sure what the correct metric is to use here
                    overall_count = sum(match_counts) / float(n_kmers)
                    matches.append(Match(
                        seq, overall_count / float(tablesize), names, float(num_matches) / n_kmers))
            
            # Add remaining tags
            for tag in set(candidates.keys()) - seen:
                matches.append(Match(tag, candidates[tag] / float(tablesize)))

        else:
            matches = [
                Match(tag, count / float(tablesize))
                for tag, count in candidates.items()
            ]
        
        return matches

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

def summarize_contaminants(outstream, matches, n_reads, header=None):
    print("", file=outstream)
    
    if header:
        print(header, file=outstream)
        print('-' * len(header), file=outstream)
    
    print("Detected {} adapters/contaminants:".format(len(matches)), file=outstream)
    
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
