# coding: utf-8
"""
Alignment module.
"""
from collections import namedtuple
from ._align import Aligner, compare_prefixes, locate
from .util import reverse_complement

# flags for global alignment

# The interpretation of the first flag is:
# An initial portion of seq1 may be skipped at no cost.
# This is equivalent to saying that in the alignment,
# gaps in the beginning of seq2 are free.
#
# The other flags have an equivalent meaning.
START_WITHIN_SEQ1 = 1
START_WITHIN_SEQ2 = 2
STOP_WITHIN_SEQ1 = 4
STOP_WITHIN_SEQ2 = 8

# Use this to get regular semiglobal alignment
# (all gaps in the beginning or end are free)
SEMIGLOBAL = START_WITHIN_SEQ1 | START_WITHIN_SEQ2 | STOP_WITHIN_SEQ1 | STOP_WITHIN_SEQ2

def compare_suffixes(s1, s2, wildcard_ref=False, wildcard_query=False):
    """
    Find out whether one string is the suffix of the other one, allowing
    mismatches. Used to find an anchored 3' adapter when no indels are allowed.
    """
    s1 = s1[::-1]
    s2 = s2[::-1]
    _, length, _, _, matches, errors = compare_prefixes(s1, s2, wildcard_ref, wildcard_query)
    return (len(s1) - length, len(s1), len(s2) - length, len(s2), matches, errors)

# Common match-result object returned by aligners

class Match(object):
    """
    TODO creating instances of this class is relatively slow and responsible for quite some runtime.
    """
    __slots__ = ['astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'front', 'adapter', 'read', 'length']
    def __init__(self, astart, astop, rstart, rstop, matches, errors, front=None, adapter=None, read=None):
        self.astart = astart
        self.astop = astop
        self.rstart = rstart
        self.rstop = rstop
        self.matches = matches
        self.errors = errors
        self.front = self._guess_is_front() if front is None else front
        self.adapter = adapter
        self.read = read
        # Number of aligned characters in the adapter. If there are indels,
        # this may be different from the number of characters in the read.
        self.length = self.astop - self.astart
        assert self.length > 0
        assert self.length - self.errors > 0
        if self.adapter:
            assert self.errors / self.length <= self.adapter.max_error_rate

    def __str__(self):
        return 'Match(astart={0}, astop={1}, rstart={2}, rstop={3}, matches={4}, errors={5})'.format(
            self.astart, self.astop, self.rstart, self.rstop, self.matches, self.errors)
    
    def copy(self):
        return Match(self.astart, self.astop, self.rstart, self.rstop, self.matches, self.errors,
                     self.front, self.adapter, self.read)
    
    def _guess_is_front(self):
        """
        Return whether this is guessed to be a front adapter.

        The match is assumed to be a front adapter when the first base of
        the read is involved in the alignment to the adapter.
        """
        return self.rstart == 0

    def wildcards(self, wildcard_char='N'):
        """
        Return a string that contains, for each wildcard character,
        the character that it matches. For example, if the adapter
        ATNGNA matches ATCGTA, then the string 'CT' is returned.

        If there are indels, this is not reliable as the full alignment
        is not available.
        """
        wildcards = [ self.read.sequence[self.rstart + i:self.rstart + i + 1] for i in range(self.length)
            if self.adapter.sequence[self.astart + i] == wildcard_char and self.rstart + i < len(self.read.sequence) ]
        return ''.join(wildcards)

    def rest(self):
        """
        Return the part of the read before this match if this is a
        'front' (5') adapter,
        return the part after the match if this is not a 'front' adapter (3').
        This can be an empty string.
        """
        if self.front:
            return self.read.sequence[:self.rstart]
        else:
            return self.read.sequence[self.rstop:]
    
    # TODO: write test
    def get_info_record(self):
        seq = self.read.sequence
        qualities = self.read.qualities
        if qualities is None:
            qualities = ''
        rsize = rsize_total = self.rstop - self.rstart
        if self.front and self.rstart > 0:
            rsize_total = self.rstop
        elif not self.front and self.rstop < len(seq):
            rsize_total = len(seq) - self.rstart
        return MatchInfo(
            self.read.name,
            self.errors,
            self.rstart,
            self.rstop,
            seq[0:self.rstart],
            seq[self.rstart:self.rstop],
            seq[self.rstop:],
            self.adapter.name,
            qualities[0:self.rstart],
            qualities[self.rstart:self.rstop],
            qualities[self.rstop:],
            self.front,
            self.astop - self.astart,
            rsize, rsize_total
        )

MatchInfo = namedtuple("MatchInfo", (
    "read_name", "errors", "rstart", "rstop", "seq_before", "seq_adapter", "seq_after",
    "adapter_name", "qual_before", "qual_adapter", "qual_after", "is_front", "asize",
    "rsize_adapter", "rsize_total"
))

# Alternative semi-global alignment (
# http://www.bioinf.uni-freiburg.de/Lehre/Courses/2013_SS/V_Bioinformatik_1/lecture4.pdf)
# strategies designed to improve insert matching of paired-end reads.
#
# Note: these are currently just prototype implementations. They will need to be optimized
# using numpy and/or re-written in cython.
#
# 1. SeqPrep algorithm: insert match algorithm that performs thresholded exhaustive
#    comparison to minimize probability of incorrect alignment. Relies on the fact that
#    overlapping reads share alleles and indels (i.e. no gaps are required) (in C++).
#    https://github.com/imgag/ngs-bits/tree/master/src/SeqPurge.
#   * Speed up sequence comparison:
#     * Between adapters and overhangs:
#       * http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080745/pdf/btu177.pdf
#     * Between reads:
#       * http://bioinformatics.oxfordjournals.org/content/30/14/2000.abstract
# 2. Skewer algorithm: bit-masked k-difference matching (in C++).
#    https://github.com/relipmoc/skewer
# 3. Quality-aware overlap alignment (in Haskell).
#    https://hackage.haskell.org/package/bio-0.5.3/docs/Bio-Alignment-QAlign.html
# 4. FOGSAA, modified for semi-global alignment.
#    http://www.nature.com/articles/srep01746
#    http://www.isical.ac.in/~bioinfo_miu/FOGSAA.7z

from bitarray import bitarray
import math

class FactorialCache(object):
    def __init__(self, init_size=150):
        self.factorials = [1] * init_size
        self.max_n = 1
        self.cur_array_size = init_size
        
    def factorial(self, n):
        if n > self.max_n:
            self._fill_upto(n)
        return self.factorials[n]
        
    def _fill_upto(self, n):
        if n >= self.cur_array_size:
            extension_size = n - self.cur_array_size + 1
            self.factorials += [1] * extension_size
        i = self.max_n
        next_i = i + 1
        while i < n:
            self.factorials[next_i] = next_i * self.factorials[i]
            i = next_i
            next_i += 1
        self.max_n = i

BASE_ENCODING = dict(
    A=bitarray('1000'),
    C=bitarray('0100'),
    G=bitarray('0010'),
    T=bitarray('0001'),
    N=bitarray('1111')
)

class OneHotEncoded:
    def __init__(self, seq=None, reverse_complement=False):
        """Create a new one-hot encoded DNA sequence from a DNA sequence."""
        if seq is not None:
            seq = seq.upper()
            self.size = len(seq)
            self.bits = bitarray()
            self.bits.encode(BASE_ENCODING, seq)
            ambig = ('1' if b == 'N' else '0' for b in seq)
            if reverse_complement:
                self.bits.reverse()
                self.ambig = bitarray(''.join(reversed(tuple(ambig))))
            else:
                self.ambig = bitarray(''.join(ambig))
    
    def __len__(self):
        return self.size
    
    def __getitem__(self, key):
        """Returns a new OneHotEncoded with the specified subsequence."""
        bits, ambig = self._subseq(key.start, key.stop)
        new_ohe = OneHotEncoded()
        new_ohe.size = len(bits) // 4
        new_ohe.bits = bits
        new_ohe.ambig = ambig
        return new_ohe
    
    def copy(self, reverse_complement=False):
        new_ohe = OneHotEncoded()
        new_ohe.size = self.size
        new_ohe.bits = self.bits.copy()
        new_ohe.ambig = self.ambig.copy()
        if reverse_complement:
            new_ohe.reverse_complement()
        return new_ohe
        
    def compare(self, other, offset=0, ambig_match=None):
        """
        Returns the number of matching bases between this sequence and `other`.
        Currently, this can handle Ns but not other IUPAC ambiguous characters.
        
        offset: If positive: the number of bases to move self forward from the
                left relative to other before comparing. If negative: the number
                of bases to move self backward from the right relative to other
                before comparing.
        ambig_match: whether to count ambiguous characters as matches (True),
                     mismatches (False) or ignore them (None).
        
        Offset and overalap are mutually exclusive. The diagram below shows
        how they work. Comparsions are made on the sequences between the pipes.
        
        offset > 0
        
        ..offset..|-------seq1--------|----------
        ----------|-------seq2--------|..offset..
        
        offset < 0
                      |--offset--|-----seq1------
        ------seq2----|--offset--|
        """
        
        l1 = len(self)
        l2 = len(other)
        if offset == 0 and (l1 == l2):
            eq = self.bits & other.bits
            ambig1 = self.ambig
            ambig2 = other.ambig
            size = l1
        else:
            start1 = start2 = 0
            if offset > 0:
                assert offset < l2
                start2 = offset
                size = min(l1, l2 - offset)
            elif offset < 0:
                offset = abs(offset)
                assert offset < max(l1, l2)
                if offset > l2:
                    start1 = offset - l2
                    size = min(l2, l1 - start1)
                else:
                    start2 = l2 - offset
                    size = min(l1, offset)
            else:
                size = min(l1, l2)
            bits1, ambig1 = self._subseq(start1, start1 + size)
            bits2, ambig2 = other._subseq(start2, start2 + size)
            eq = bits1 & bits2
        
        # when ambiguous bases match, they contribute 4 counts rather than 1
        matches = eq.count() - ((ambig1 & ambig2).count() * 3)
        if ambig_match is True:
            return (matches, size)
        else:
            num_ambig = (ambig1 | ambig2).count()
            matches -= num_ambig
            if ambig_match is False:
                return (matches, size)
            else:
                return (matches, size - num_ambig)
    
    def _subseq(self, start, stop):
        return (
            self.bits[(start*4):(stop*4)],
            self.ambig[start:stop]
        )
    
    def reverse_complement(self):
        """Reverse-complements the sequence in place."""
        self.bits.reverse()
        self.ambig.reverse()
                
    def to_acgt(self):
        """Returns the sequence as a list of DNA base characters."""
        return self.bits.decode(BASE_ENCODING)
    
    def translate(self):
        """
        Returns a list of tuples where, at position i, the first element is the
        four-bit encoding of the base at i, and the second element is the character
        corresponding to that base.
        """
        binstr = self.bits.to01()
        basestr = self.to_acgt()
        return ((binstr[(i*4):((i+1)*4)], basestr[i])
            for i in range(len(self)))
    
    def __repr__(self):
        """Returns the original DNA base string."""
        return "".join(self.to_acgt())
    
    def __str__(self):
        """
        Returns a string representation of the one-hot encoding in matrix form,
        along with the equivalent base character for each row.
        """
        return "ACGT\n----\n" + "\n".join("{} : {}".format(*t) for t in self.translate())

class SeqPurgeAligner(object):
    """
    Implementation of the SeqPurge insert matching algorithm.
    This only works with paired-end reads with 3' adapters.
    """
    def __init__(self, max_error_prob=1E-6, min_insert_overlap=1, min_insert_match_frac=0.8,
                 min_adapter_overlap=1, min_adapter_match_frac=0.8, adapter_check_cutoff=9):
        self.max_error_prob = max_error_prob
        self.min_insert_overlap = min_insert_overlap
        self.min_insert_match_frac = float(min_insert_match_frac)
        self.max_insert_mismatch_frac = 1.0 - min_insert_match_frac
        self.min_adapter_overlap = min_adapter_overlap
        self.min_adapter_match_frac = float(min_adapter_match_frac)
        self.max_adapter_mismatch_frac = 1.0 - self.min_adapter_match_frac
        self.adapter_check_cutoff = adapter_check_cutoff
        self.factorial_cache = FactorialCache()
    
    def match_insert1(self, seq1, seq2, adapter1, adapter2):
        """Use one-hot encoding for insert and adpater matching"""
        l1 = len(seq1)
        l2 = len(seq2)
        seq_len = min(l1, l2)
        if l1 > l2:
            seq1 = seq1[:l2]
        elif l2 > l1:
            seq2 = seq1[:l1]
        
        seq1 = OneHotEncoded(seq1)
        seq2 = OneHotEncoded(seq2, reverse_complement=True)
            
        best_offset = best_a1_match = best_a2_match = None
        best_p = self.max_error_prob
        insert_match = None
        
        for offset in range(self.min_insert_overlap, seq_len):
            # Count number of matches between the two sequences.
            # size = total number of bases compared (excluding Ns)
            matches, size = seq1.compare(seq2, offset)
        
            # Check that at least min_insert_overlap base were compared and that
            # there are fewer than the maximium number of allowed mismatches
            if size < self.min_insert_overlap or matches < (size * self.min_insert_match_frac):
                continue
        
            # Compute probability of seeing this number of matches at random and
            # check that it's less than the maximum allowed
            p = self.match_probability(matches, size)
            if p >= best_p:
                continue
        
            # Otherwise, we've found an acceptable insert match
            if insert_match is None:
                insert_match = size
        
            if (seq_len - size) < self.min_adapter_overlap:
                continue
            
            # Match the overhang against the adapter sequence
            a1_match = seq1.compare(adapter1, offset=-offset)
            a2_match = seq2.compare(adapter2, offset=-offset)
            min_adapter_matches = math.ceil(offset * self.min_adapter_match_frac)
            if a1_match[0] < min_adapter_matches and a2_match[0] < min_adapter_matches:
                continue
            a1_prob = self.match_probability(*a1_match)
            a2_prob = self.match_probability(*a2_match)
            if offset > self.adapter_check_cutoff and (a1_prob * a2_prob) > self.max_error_prob:
                continue
        
            best_p = p
            best_offset = offset
            # This is a shortcut - if the match isn't symmetrical (e.g. if bases are trimmed off
            # the 5' end of one read but not the other), this will screw up the adapter stats. We
            # can solve this by always aligning adapters, at the expense of performance
            best_match = a1_match if a1_prob < a2_prob else a2_match
        
        if best_offset is None or best_offset < self.min_adapter_overlap:
            return (None, None, insert_match)
        
        insert_match_size = seq_len - best_offset
        adapter_len1 = min(len(adapter1), l1 - insert_match_size)
        adapter_len2 = min(len(adapter2), l2 - insert_match_size)
        best_matches = best_match[0]
        best_mismatches = best_match[1] - best_matches
        match1 = Match(0, adapter_len1, insert_match_size, l1, best_matches, best_mismatches)
        match2 = Match(0, adapter_len2, insert_match_size, l2, best_matches, best_mismatches)
        return (match1, match2, insert_match)
    
    def match_insert2(self, seq1, seq2, adapter1, adapter2):
        """Use cutadapt aligner for insert and adapter matching"""
        alen1 = len(adapter1)
        alen2 = len(adapter2)
        l1 = len(seq1)
        l2 = len(seq2)
        seq_len = min(l1, l2)
        if l1 > l2:
            seq1 = seq1[:l2]
        elif l2 > l1:
            seq2 = seq1[:l1]
        
        seq2_rc = reverse_complement(seq2)
        result = [None, None, 0]
        
        aligner = Aligner(seq2_rc, 1 - self.min_insert_match_frac,
                          START_WITHIN_SEQ1 | STOP_WITHIN_SEQ2, False, False)
        aligner.min_overlap = self.min_insert_overlap
        aligner.indel_cost = 100000
        insert_match = aligner.locate(seq1)
        
        if insert_match:
            offset = min(insert_match[0], seq_len - insert_match[3])
            insert_match_size = seq_len - offset
            result[2] = insert_match_size
            
            if (offset < self.min_adapter_overlap or
                    self.match_probability(insert_match[4], insert_match_size) > self.max_error_prob):
                return result
            
            a1_match = compare_prefixes(seq1[insert_match_size:], adapter1)
            a2_match = compare_prefixes(seq2[insert_match_size:], adapter2)
            adapter_len = min(offset, alen1, alen2)
            min_adapter_matches = math.ceil(adapter_len * self.min_adapter_match_frac)
            if a1_match[4] < min_adapter_matches and a2_match[4] < min_adapter_matches:
                return result
            a1_prob = self.match_probability(a1_match[4], adapter_len)
            a2_prob = self.match_probability(a2_match[4], adapter_len)
            if adapter_len > self.adapter_check_cutoff and (a1_prob * a2_prob) > self.max_error_prob:
                return result
            
            adapter_len1 = min(alen1, l1 - insert_match_size)
            adapter_len2 = min(alen2, l2 - insert_match_size)
            best_adapter_matches, best_adapter_mismatches = (a1_match if a1_prob < a2_prob else a2_match)[4:6]
            result[0] = Match(0, adapter_len1, insert_match_size, l1, best_adapter_matches, best_adapter_mismatches)
            result[1] = Match(0, adapter_len2, insert_match_size, l2, best_adapter_matches, best_adapter_mismatches)
        
        return result
    
    def match_insert3(self, seq1, seq2, adapter1, adapter2):
        """Use the original seqpurge implementation (iterative comparison)"""
        alen1 = len(adapter1)
        alen2 = len(adapter2)
        slen1 = len(seq1)
        slen2 = len(seq2)
        seq_len = min(slen1, slen2)
        if slen1 > seq_len:
            seq1 = seq1[:seq_len]
        elif slen2 > seq_len:
            seq2 = seq1[:seq_len]
        seq2_rc = reverse_complement(seq2)
        
        best_offset = best_matches = best_mismatches = None
        best_p = self.max_error_prob
        best_insert_match = None
        
        def compare_insert(offset):
            size = seq_len - offset
            max_mismatches = self.max_insert_mismatch_frac * size
            mismatches = 0
            invalid = 0
            for b1, b2 in zip(seq1[:size], seq2_rc[offset:]):
                if b1 == 'N' or b2 == 'N':
                    invalid += 1
                elif b1 != b2:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        return None
            size -= invalid
            return (size - mismatches, size)
        
        for offset in range(self.min_insert_overlap, seq_len):
            insert_match = compare_insert(offset)
            if not insert_match:
                continue
            
            matches, size = insert_match
            if size < self.min_insert_overlap:
                continue
        
            # Compute probability of seeing this number of matches at random and
            # check that it's less than the maximum allowed
            p = self.match_probability(matches, size)
            if p >= best_p:
                continue
        
            # Otherwise, we've found an acceptable insert match
            if best_insert_match is None:
                best_insert_match = size
        
            if (seq_len - size) < self.min_adapter_overlap:
                continue
            
            a1_match = compare_prefixes(seq1[size:], adapter1)
            a2_match = compare_prefixes(seq2[size:], adapter2)
            adapter_len = min(seq_len - size, alen1, alen2)
            min_adapter_matches = math.ceil(adapter_len * self.min_adapter_match_frac)
            if a1_match[4] < min_adapter_matches and a2_match[4] < min_adapter_matches:
                continue
            a1_prob = self.match_probability(a1_match[4], adapter_len)
            a2_prob = self.match_probability(a2_match[4], adapter_len)
            if adapter_len > self.adapter_check_cutoff and (a1_prob * a2_prob) > self.max_error_prob:
                continue
            
            best_p = p
            best_offset = offset
            # This is a shortcut - if the match isn't symmetrical (e.g. if bases are trimmed off
            # the 5' end of one read but not the other), this will screw up the adapter stats. We
            # can solve this by always aligning adapters, at the expense of performance
            best_adapter_matches, best_adapter_mismatches = (a1_match if a1_prob < a2_prob else a2_match)[4:6]
        
        if best_offset is None or best_offset < self.min_adapter_overlap:
            return (None, None, best_insert_match)
        
        insert_match_size = seq_len - best_offset
        adapter_len1 = min(alen1, slen1 - insert_match_size)
        adapter_len2 = min(alen2, slen2 - insert_match_size)
        return (
            Match(0, adapter_len1, insert_match_size, slen1, best_adapter_matches, best_adapter_mismatches),
            Match(0, adapter_len2, insert_match_size, slen2, best_adapter_matches, best_adapter_mismatches),
            insert_match_size
        )

    match_insert = match_insert2
    
    def match_adapter(self, seq, adapter):
        """
        Try to match the adapter to the read sequence, starting at the 5' end.
        Note: this does not attempt to deal with indels.
        Returns the position at which the adapter starts within the read, or
        None if no match was found.
        """
        seq_len = len(seq)
        seq = OneHotEncoded(seq)
        for offset in range(0, seq_len - self.min_adapter_overlap):
            matches, size = adapter.compare(seq, offset)
            if matches < (size * self.min_adapter_match_frac):
                continue
            if self.match_probability(matches, size) > self.max_error_prob:
                continue
            adapter_len = min(len(adapter), len(seq) - offset)
            return Match(0, adapter_len, offset, offset + adapter_len, matches, size - matches)
        return None
    
    def match_probability(self, matches, size):
        nfac = self.factorial_cache.factorial(size)
        p = 0.0
        i = matches
        
        while i <= size:
            p += (
                (0.75 ** (size - i)) *
                (0.25 ** i) *
                nfac /
                self.factorial_cache.factorial(i) /
                self.factorial_cache.factorial(size - i)
            )
            i += 1

        return p
