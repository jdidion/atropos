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
        if self.length <= 0:
            raise ValueError('Match length must be >= 0')
        if self.length - self.errors <= 0:
            raise ValueError('A Match requires at least one matching position.')
        # TODO: this assertion may not always hold now that we use both error
        # rates and probabilities
        #if self.adapter:
        #    assert self.errors / self.length <= self.adapter.max_error_rate

    def __repr__(self):
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

    # def __repr__(self) -> str:
    #     return f"InsertAligner<adapter1={self.adapter1}, adapter2={self.adapter2}, " \
    #            f"match_probability={self.match_probability}, " \
    #            f"insert_max_rmp={self.insert_max_rmp}, " \
    #            f"adapter_max_rmp={self.adapter_max_rmp}, " \
    #            f"min_insert_overlap={self.min_insert_overlap}, " \
    #            f"max_insert_mismatch_frac={self.max_insert_mismatch_frac}, " \
    #            f"min_adapter_overlap={self.min_adapter_overlap}, " \
    #            f"max_adapter_mismatch_frac={self.max_adapter_mismatch_frac}, " \
    #            f"adapter_check_cutoff={self.adapter_check_cutoff}, " \
    #            f"base_probs={self.base_probs}, " \
    #            f"adapter_wildcards={self.adapter_wildcards}, " \
    #            f"read_wildcards={self.read_wildcards}, " \
    #            f"aligner={self.aligner}>"

    def match_insert(self, seq1, seq2):
        """Use cutadapt aligner for insert and adapter matching.

        Args:
            seq1, seq2: Sequences to match.

        Returns:
            A :class:`Match` object, or None if there is no match.
        """
        seq_len1 = len(seq1)
        seq_len2 = len(seq2)
        seq_len = min(seq_len1, seq_len2)
        if seq_len1 > seq_len2:
            seq1 = seq1[:seq_len2]
        elif seq_len2 > seq_len1:
            seq2 = seq2[:seq_len1]

        seq2_rc = reverse_complement(seq2)

        def _match(_insert_match, _offset, _insert_match_size, _):
            if _offset < self.min_adapter_overlap:
                # The reads are mostly overlapping, to the point where
                # there's not enough overhang to do a confident adapter
                # match. We return just the insert match to signal that
                # error correction can be done even though no adapter
                # trimming is required.
                return (_insert_match, None, None)

            # TODO: this is very sensitive to the exact correct choice of
            # adapter. For example, if you specifiy GATCGGAA... and the correct
            # adapter is AGATCGGAA..., the prefixes will not match exactly and
            # the alignment will fail. We need to use a comparison that is a bit
            # more forgiving.

            def _adapter_match(insert_seq, adapter_seq, adapter_len):
                amatch = compare_prefixes(
                    insert_seq[_insert_match_size:], adapter_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards)
                alen = min(_offset, adapter_len)
                return amatch, alen, round(alen * self.max_adapter_mismatch_frac)

            a1_match, a1_length, a1_max_mismatches = _adapter_match(
                seq1, self.adapter1, self.adapter1_len)
            a2_match, a2_length, a2_max_mismatches = _adapter_match(
                seq2, self.adapter2, self.adapter2_len)

            if (
                    a1_match[5] > a1_max_mismatches and
                    a2_match[5] > a2_max_mismatches):
                return None

            if min(a1_length, a2_length) > self.adapter_check_cutoff:
                a1_prob = self.match_probability(a1_match[4], a1_length)
                a2_prob = self.match_probability(a2_match[4], a2_length)
                if (a1_prob * a2_prob) > self.adapter_max_rmp:
                    return None

            mismatches = min(a1_match[5], a2_match[5])

            def _create_match(alen, slen):
                alen = min(alen, slen - _insert_match_size)
                _mismatches = min(alen, mismatches)
                _matches = alen - _mismatches
                return Match(0, alen, _insert_match_size, slen, _matches, _mismatches)

            return (
                _insert_match,
                _create_match(a1_length, seq_len1),
                _create_match(a2_length, seq_len2)
            )

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
                    # filtered_matches.sort(key=lambda x: x[2], reverse=True)
                    filtered_matches.sort(key=lambda x: x[3])
                    for match_args in filtered_matches:
                        match = _match(*match_args)
                        if match:
                            return match

            return None
