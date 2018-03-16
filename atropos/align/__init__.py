# coding: utf-8
"""
Alignment module.
"""
from collections import namedtuple
from atropos.align._align import Aligner, MultiAligner, compare_prefixes, locate
from atropos.util import RandomMatchProbability, reverse_complement

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
SEMIGLOBAL = (
    START_WITHIN_SEQ1 | START_WITHIN_SEQ2 | 
    STOP_WITHIN_SEQ1 | STOP_WITHIN_SEQ2)

def compare_suffixes(
        suffix_ref, suffix_query, wildcard_ref=False, wildcard_query=False):
    """Find out whether one string is the suffix of the other one, allowing
    mismatches. Used to find an anchored 3' adapter when no indels are allowed.
    
    Args:
        suffix_ref, suffix_query: The suffices to compare.
        wildcard_ref, wildcard_query: Whether wildcards are valid in either of
            the suffices.
    """
    suffix_ref = suffix_ref[::-1]
    suffix_query = suffix_query[::-1]
    _, length, _, _, matches, errors = compare_prefixes(
        suffix_ref, suffix_query, wildcard_ref, wildcard_query)
    return (
        len(suffix_ref) - length, len(suffix_ref), len(suffix_query) - length,
        len(suffix_query), matches, errors)

# Common match-result object returned by aligners

# TODO creating instances of this class is relatively slow and responsible for
# quite some runtime.

class Match(object):
    """An alignment match.
    
    Args:
        astart: Starting position of the match within the adapter.
        astop: Ending position of the match within the adapter.
        rstart: Starting position of the match within the read.
        rstop: Ending position of the match within the read.
        matches: Number of matching bases.
        errors: Number of mismatching bases.
        front: Whether the match is to the front of the read.
        adapter: The :class:`Adapter`.
        read: The :class:`Sequence`.
        length: The length of the match.
    """
    __slots__ = [
        'astart', 'astop', 'rstart', 'rstop', 'matches', 'errors', 'front',
        'adapter', 'read', 'length']
    
    def __init__(
            self, astart, astop, rstart, rstop, matches, errors,
            front=None, adapter=None, read=None):
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
        # TODO: this assertion may not always hold now that we use both error
        # rates and probabilities
        #if self.adapter:
        #    assert self.errors / self.length <= self.adapter.max_error_rate

    def __str__(self):
        return (
            'Match(astart={0}, astop={1}, rstart={2}, rstop={3}, matches={4}, '
            'errors={5})').format(
                self.astart, self.astop, self.rstart, self.rstop, self.matches,
                self.errors)
    
    def copy(self):
        """Create a copy of this Match.
        """
        return Match(
            self.astart, self.astop, self.rstart, self.rstop, self.matches,
            self.errors, self.front, self.adapter, self.read)
    
    def _guess_is_front(self):
        """Return whether this is guessed to be a front adapter.
        
        The match is assumed to be a front adapter when the first base of
        the read is involved in the alignment to the adapter.
        """
        return self.rstart == 0
    
    def wildcards(self, wildcard_char='N'):
        """Return a string that contains, for each wildcard character,
        the character that it matches. For example, if the adapter
        ATNGNA matches ATCGTA, then the string 'CT' is returned.

        If there are indels, this is not reliable as the full alignment
        is not available.
        """
        wildcards = [
            self.read.sequence[self.rstart + i]
            for i in range(self.length)
            if (self.adapter.sequence[self.astart + i] == wildcard_char and
                self.rstart + i < len(self.read.sequence))
        ]
        return ''.join(wildcards)

    def rest(self):
        """Returns the part of the read before this match if this is a
        'front' (5') adapter, or the part after the match if this is not
        a 'front' adapter (3'). This can be an empty string.
        """
        if self.front:
            return self.read.sequence[:self.rstart]
        else:
            return self.read.sequence[self.rstop:]
    
    # TODO: write test
    def get_info_record(self):
        """Returns a :class:`MatchInfo`, which contains information about the
        match to write into the info file.
        """
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
            rsize, rsize_total)

MatchInfo = namedtuple("MatchInfo", (
    "read_name", "errors", "rstart", "rstop", "seq_before", "seq_adapter",
    "seq_after", "adapter_name", "qual_before", "qual_adapter", "qual_after",
    "is_front", "asize", "rsize_adapter", "rsize_total"))

# Alternative semi-global alignment (
# http://www.bioinf.uni-freiburg.de/Lehre/Courses/2013_SS/V_Bioinformatik_1/lecture4.pdf)
# strategies designed to improve insert matching of paired-end reads.
#
# Note: these are currently just prototype implementations. They will need to be optimized
# using numpy and/or re-written in cython.
#
# 1. SeqPurge algorithm: insert match algorithm that performs thresholded exhaustive
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
# 5. EDLIB: edit distance-based alignment
#    https://github.com/Martinsos/edlib
# 6. Phred-adjusted ML for error probability: 
# https://biosails.github.io/pheniqs/glossary.html#phred_adjusted_maximum_likelihood_decoding
# 7. Adaptive banded alignment: https://github.com/ocxtal/libgaba
# Also think about different sequence encodings that might enable faster alignment
# https://github.com/hammerlab/kerseq/blob/master/kerseq/sequence_encoding.py
# 8. https://github.com/yamada-kd/nepal
# 9. The SeqAn C++ library implements several alignment algorithms: 
# http://www.sciencedirect.com/science/article/pii/S0168165617315420
# 10. Could we treat paired end read + adapter alignment as an MSA problem?
# 11. Look at alignment-free tools for pairwise sequence comparison:
# * http://www.combio.pl/alfree/tools/
# * http://www.combio.pl/alfree
# * http://bioinformatics.org.au/tools/decaf+py/

class InsertAligner(object):
    """Implementation of an insert matching algorithm.
    
    If the inserts align, the overhangs are searched for the adapter sequences.
    Otherwise, each read is search for its adapter separately.
    
    This only works with paired-end reads with 3' adapters.
    
    Args:
        adapter1, adapter2: read1, read2 adapters.
        match_probability: Callable that calculates random match probability
            given arguments (num_matches, match_len).
        insert_max_rmp: Max random match probability for the insert match.
        adapter_max_rmp: Max random match probability for the adapter match.
        min_insert_overlap: Minimum number of bases the inserts must overlap
            to be considered matching.
        max_insert_mismatch_frac: Maximum fraction of mismatching bases between
            the inserts to be considered matching.
        min_adapter_overlap: Minimum number of bases the adapter must overlap
            the read to be considered matching.
        max_adapter_mismatch_frac: Maximum fraction of mismatching bases between
            the adapter and the read to be considered matching.
        adapter_check_cutoff: Threshold number of matching bases required before
            adapter matches are checked against random match probability.
        base_probs: Dict of (match_prob, mismatch_prob), which are the
            probabilities passed to
            :method:`atropos.util.RandomMatchProbability.__call__()`.
    """
    def __init__(
            self, adapter1, adapter2,
            match_probability=RandomMatchProbability(),
            insert_max_rmp=1E-6, adapter_max_rmp=0.001,
            min_insert_overlap=1, max_insert_mismatch_frac=0.2,
            min_adapter_overlap=1, max_adapter_mismatch_frac=0.2,
            adapter_check_cutoff=9, base_probs=None,
            adapter_wildcards=True, read_wildcards=False):
        self.adapter1 = adapter1
        self.adapter1_len = len(adapter1)
        self.adapter2 = adapter2
        self.adapter2_len = len(adapter2)
        self.match_probability = match_probability
        self.insert_max_rmp = insert_max_rmp
        self.adapter_max_rmp = adapter_max_rmp
        self.min_insert_overlap = min_insert_overlap
        self.max_insert_mismatch_frac = float(max_insert_mismatch_frac)
        self.min_adapter_overlap = min_adapter_overlap
        self.max_adapter_mismatch_frac = float(max_adapter_mismatch_frac)
        self.adapter_check_cutoff = adapter_check_cutoff
        self.base_probs = base_probs or dict(
            match_prob=0.25, mismatch_prob=0.75)
        self.adapter_wildcards = adapter_wildcards
        self.read_wildcards = read_wildcards
        self.aligner = MultiAligner(
            max_insert_mismatch_frac,
            START_WITHIN_SEQ1 | STOP_WITHIN_SEQ2,
            min_insert_overlap)
    
    def match_insert(self, seq1, seq2):
        """Use cutadapt aligner for insert and adapter matching.
        
        Args:
            seq1, seq2: Sequences to match.
        
        Returns:
            A :class:`Match` object, or None if there is no match.
        """
        len1 = len(seq1)
        len2 = len(seq2)
        seq_len = min(len1, len2)
        if len1 > len2:
            seq1 = seq1[:len2]
        elif len2 > len1:
            seq2 = seq2[:len1]

        seq2_rc = reverse_complement(seq2)

        def _match(insert_match, offset, insert_match_size, prob): # pylint disable=unused-argument
            if offset < self.min_adapter_overlap:
                # The reads are mostly overlapping, to the point where
                # there's not enough overhang to do a confident adapter
                # match. We return just the insert match to signal that
                # error correction can be done even though no adapter
                # trimming is required.
                return (insert_match, None, None)
            
            # TODO: this is very sensitive to the exact correct choice of
            # adapter. For example, if you specifiy GATCGGAA... and the correct
            # adapter is AGATCGGAA..., the prefixes will not match exactly and
            # the alignment will fail. We need to use a comparison that is a bit
            # more forgiving.
            
            a1_match = compare_prefixes(
                seq1[insert_match_size:], self.adapter1,
                wildcard_ref=self.adapter_wildcards,
                wildcard_query=self.read_wildcards)
            a2_match = compare_prefixes(
                seq2[insert_match_size:], self.adapter2,
                wildcard_ref=self.adapter_wildcards,
                wildcard_query=self.read_wildcards)
            adapter_len = min(offset, self.adapter1_len, self.adapter2_len)
            max_adapter_mismatches = round(
                adapter_len * self.max_adapter_mismatch_frac)
            if (
                    a1_match[5] > max_adapter_mismatches and
                    a2_match[5] > max_adapter_mismatches):
                return None
            
            a1_prob = self.match_probability(a1_match[4], adapter_len)
            a2_prob = self.match_probability(a2_match[4], adapter_len)
            if (
                    (adapter_len > self.adapter_check_cutoff) and
                    ((a1_prob * a2_prob) > self.adapter_max_rmp)):
                return None

            adapter_len1 = min(self.adapter1_len, len1 - insert_match_size)
            adapter_len2 = min(self.adapter2_len, len2 - insert_match_size)
            best_adapter_matches, best_adapter_mismatches = (
                a1_match if a1_prob < a2_prob else a2_match)[4:6]
            
            return (
                insert_match,
                Match(
                    0, adapter_len1, insert_match_size, len1,
                    best_adapter_matches, best_adapter_mismatches),
                Match(
                    0, adapter_len2, insert_match_size, len2,
                    best_adapter_matches, best_adapter_mismatches))
        
        # # This is the old way of doing things, where we use the built-in
        # # Aligner to do a single match.
        # aligner = Aligner(
        #     seq2_rc,
        #     self.max_insert_mismatch_frac,
        #     START_WITHIN_SEQ1 | STOP_WITHIN_SEQ2,
        #     False, False)
        # aligner.min_overlap = self.min_insert_overlap
        # aligner.indel_cost = 100000
        #
        # insert_match = aligner.locate(seq1)
        #
        # if not insert_match:
        #     return None
        #
        # offset = min(insert_match[0], seq_len - insert_match[3])
        # insert_match_size = seq_len - offset
        # prob = self.match_probability(insert_match[4], insert_match_size)
        #
        # if prob > self.insert_max_rmp:
        #     return None
        #
        # return _match(insert_match, offset, insert_match_size, prob)
        
        # Use an aligner that returns all matches that satisfy the
        # overlap and error rate thresholds. We sort by matches and
        # then mismatches, and then check each in turn until we find
        # one with an adapter match (if any).

        insert_matches = self.aligner.locate(seq2_rc, seq1)

        if insert_matches:
            # Filter by random-match probability
            filtered_matches = []
            for insert_match in insert_matches:
                offset = min(insert_match[0], seq_len - insert_match[3])
                insert_match_size = seq_len - offset
                prob = self.match_probability(insert_match[4], insert_match_size, **self.base_probs)
                if prob <= self.insert_max_rmp:
                    filtered_matches.append((insert_match, offset, insert_match_size, prob))
            
            if filtered_matches:
                if len(filtered_matches) == 1:
                    return _match(*filtered_matches[0])
                else:
                    # Test matches in order of random-match probability.
                    # TODO: compare against sorting by length (which is how
                    # SeqPurge essentially does it).
                    #filtered_matches.sort(key=lambda x: x[2], reverse=True)
                    filtered_matches.sort(key=lambda x: x[3])
                    for match_args in filtered_matches:
                        match = _match(*match_args)
                        if match:
                            return match
            
            return None
