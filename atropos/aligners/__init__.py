from enum import IntFlag
from typing import Dict, Optional, Tuple, TypeVar, Union

from atropos.aligners._align import Aligner, MultiAligner, compare_prefixes, locate
from atropos.utils.ngs import reverse_complement
from atropos.utils.statistics import RandomMatchProbability


DEFAULT_INSERT_MAX_RMP = 1e-6
"""Default value for InsertAligner max insert random match probability"""
DEFAULT_ADAPTER_MAX_RMP = 0.001
"""Default value for InsertAligner max adapter random match probability"""
DEFAULT_MIN_INSERT_OVERLAP = 1
"""Default value for InsertAligner min insert overlap"""
DEFAULT_MAX_INSERT_MISMATCH_FRAC = 0.2
"""Default value for InsertAligner max insert mismatch fraction"""
DEFAULT_MIN_ADAPTER_OVERLAP = 1
"""Default value for InsertAligner min adapter overlap"""
DEFAULT_MAX_ADAPTER_MISMATCH_FRAC = 0.2
"""Default value for InsertAligner max adapter mismatch fraction"""
DEFAULT_ADAPTER_CHECK_CUTOFF = 9
"""Default value for InsertAligner adapter check cutoff"""


class GapRule(IntFlag):
    START_WITHIN_SEQ1 = 1
    """Gaps at begining of seq2 have no penalty."""
    START_WITHIN_SEQ2 = 2
    """Gaps at begining of seq1 have no penalty."""
    STOP_WITHIN_SEQ1 = 4
    """Gaps at end of seq2 have no penalty."""
    STOP_WITHIN_SEQ2 = 8
    """Gaps at end of seq1 have no penalty"""
    SEMIGLOBAL = (
        START_WITHIN_SEQ1
        | START_WITHIN_SEQ2
        | STOP_WITHIN_SEQ1
        | STOP_WITHIN_SEQ2
    )
    """Typical semiglobal alignment (all gaps in the beginning or end are free)"""


MatchTuple = Tuple[int, int, int, int, int, int]


M = TypeVar("M", bound="Match")


class Match:
    """
    Common match-result object returned by aligners.
    """

    __slots__ = [
        "astart",
        "astop",
        "rstart",
        "rstop",
        "matches",
        "errors",
        "front",
        "length",
    ]

    def __init__(
        self,
        astart: int,
        astop: int,
        rstart: int,
        rstop: int,
        matches: int,
        errors: int,
        front: Optional[bool] = None,
    ):
        """
        Args:
            astart: Starting position of the match within the adapter.
            astop: Ending position of the match within the adapter.
            rstart: Starting position of the match within the read.
            rstop: Ending position of the match within the read.
            matches: Number of matching bases.
            errors: Number of mismatching bases.
            front: Whether the match is to the front of the read.
        """
        self.astart = astart
        self.astop = astop
        self.rstart = rstart
        self.rstop = rstop
        self.matches = matches
        self.errors = errors
        self.front = self._guess_is_front() if front is None else front
        # Number of aligned characters in the adapter. If there are indels,
        # this may be different from the number of characters in the read.
        self.length = self.astop - self.astart
        if self.length <= 0:
            raise ValueError("Match length must be >= 0")
        if self.length - self.errors <= 0:
            raise ValueError("A Match requires at least one matching position.")
        # TODO: this assertion may not always hold now that we use both error
        #  rates and probabilities
        # if self.adapter:
        #    assert self.errors / self.length <= self.adapter.max_error_rate

    def __str__(self) -> str:
        return (
            f"Match(astart={self.astart}, astop={self.astop}, "
            f"rstart={self.rstart}, rstop={self.rstop}, matches={self.matches}, "
            f"errors={self.errors})"
        )

    def copy(self: M) -> M:
        """
        Create a copy of this Match.
        """
        return Match(
            self.astart,
            self.astop,
            self.rstart,
            self.rstop,
            self.matches,
            self.errors,
            self.front,
        )

    def _guess_is_front(self) -> bool:
        """
        Guesses whether this is match is for a front adapter.

        The match is assumed to be a front adapter when the first base of
        the read is involved in the alignment to the adapter.
        """
        return self.rstart == 0


InsertMatchResult = Tuple[MatchTuple, Optional[Match], Optional[Match]]


class InsertAligner:
    """
    Implementation of an insert matching algorithm.

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
        self,
        adapter1: str,
        adapter2: str,
        match_probability: Optional[RandomMatchProbability] = None,
        insert_max_rmp: float = DEFAULT_INSERT_MAX_RMP,
        adapter_max_rmp: float = DEFAULT_ADAPTER_MAX_RMP,
        min_insert_overlap: int = DEFAULT_MIN_INSERT_OVERLAP,
        max_insert_mismatch_frac: Union[str, float] = DEFAULT_MAX_INSERT_MISMATCH_FRAC,
        min_adapter_overlap: int = DEFAULT_MIN_ADAPTER_OVERLAP,
        max_adapter_mismatch_frac: Union[str, float] =
        DEFAULT_MAX_ADAPTER_MISMATCH_FRAC,
        adapter_check_cutoff: int = DEFAULT_ADAPTER_CHECK_CUTOFF,
        base_probs: Optional[Dict[str, float]] = None,
        adapter_wildcards: bool = True,
        read_wildcards: bool = False,
    ):
        self.adapter1 = adapter1
        self.adapter1_len = len(adapter1)
        self.adapter2 = adapter2
        self.adapter2_len = len(adapter2)
        self.match_probability = match_probability or RandomMatchProbability()
        self.insert_max_rmp = insert_max_rmp
        self.adapter_max_rmp = adapter_max_rmp
        self.min_insert_overlap = min_insert_overlap
        self.max_insert_mismatch_frac = float(max_insert_mismatch_frac)
        self.min_adapter_overlap = min_adapter_overlap
        self.max_adapter_mismatch_frac = float(max_adapter_mismatch_frac)
        self.adapter_check_cutoff = adapter_check_cutoff
        self.base_probs = base_probs or dict(match_prob=0.25, mismatch_prob=0.75)
        self.adapter_wildcards = adapter_wildcards
        self.read_wildcards = read_wildcards
        self.aligner = MultiAligner(
            max_insert_mismatch_frac,
            GapRule.START_WITHIN_SEQ1 | GapRule.STOP_WITHIN_SEQ2,
            min_insert_overlap,
        )

    def match_insert(self, seq1: str, seq2: str) -> Optional[InsertMatchResult]:
        """
        Uses cutadapt aligner for insert and adapter matching.

        Args:
            seq1, seq2: Sequences to match.

        Returns:
            A :class:`Match` object, or None if there is no match.
        """
        seq_len1 = len(seq1)
        seq_len2 = len(seq2)

        if seq_len1 > seq_len2:
            seq1 = seq1[:seq_len2]
            seq_len = seq_len2
        else:
            seq_len = seq_len1
            if seq_len2 > seq_len1:
                seq2 = seq2[:seq_len1]

        seq2_rc = reverse_complement(seq2)

        def _match(
            _insert_match: MatchTuple, _offset: int, _insert_match_size: int
        ) -> Optional[InsertMatchResult]:
            if offset < self.min_adapter_overlap:
                # The reads are mostly overlapping, to the point where
                # there's not enough overhang to do a confident adapter
                # match. We return just the insert match to signal that
                # error correction can be done even though no adapter
                # trimming is required.
                return _insert_match, None, None

            # TODO: this is very sensitive to the exact correct choice of adapter.
            #  For example, if you specifiy GATCGGAA... and the correct adapter is
            #  AGATCGGAA..., the prefixes will not match exactly and the alignment
            #  will fail. We need to use a comparison that is a bit more forgiving.
            def _adapter_match(
                insert_seq: str, adapter_seq: str, adapter_len: int
            ) -> Tuple[MatchTuple, int, float]:
                amatch = compare_prefixes(
                    insert_seq[_insert_match_size:],
                    adapter_seq,
                    wildcard_ref=self.adapter_wildcards,
                    wildcard_query=self.read_wildcards
                )
                alen = min(_offset, adapter_len)
                return amatch, alen, round(alen * self.max_adapter_mismatch_frac)

            a1_match, a1_length, a1_max_mismatches = _adapter_match(
                seq1, self.adapter1, self.adapter1_len
            )
            a2_match, a2_length, a2_max_mismatches = _adapter_match(
                seq2, self.adapter2, self.adapter2_len
            )

            if (
                a1_match[5] > a1_max_mismatches and
                a2_match[5] > a2_max_mismatches
            ):
                return None

            if min(a1_length, a2_length) > self.adapter_check_cutoff:
                a1_prob = self.match_probability(a1_match[4], a1_length)
                a2_prob = self.match_probability(a2_match[4], a2_length)
                if (a1_prob * a2_prob) > self.adapter_max_rmp:
                    return None

            mismatches = min(a1_match[5], a2_match[5])

            def _create_match(alen: int, slen: int) -> Match:
                alen = min(alen, slen - _insert_match_size)
                _mismatches = min(alen, mismatches)
                _matches = alen - _mismatches
                return Match(0, alen, _insert_match_size, slen, _matches, _mismatches)

            return (
                _insert_match,
                _create_match(a1_length, seq_len1),
                _create_match(a2_length, seq_len2)
            )

        insert_matches = self.aligner.locate(seq2_rc, seq1)
        if insert_matches:
            # Filter by random-match probability
            filtered_matches = []

            for insert_match in insert_matches:
                offset = min(insert_match[0], seq_len - insert_match[3])
                insert_match_size = seq_len - offset
                prob = self.match_probability(
                    insert_match[4], insert_match_size, **self.base_probs
                )
                if prob <= self.insert_max_rmp:
                    filtered_matches.append(
                        (insert_match, offset, insert_match_size, prob)
                    )

            if filtered_matches:
                if len(filtered_matches) == 1:
                    return _match(*filtered_matches[0][0:3])
                else:
                    # Test matches in order of random-match probability.
                    # TODO: compare against sorting by length (which is how
                    #  SeqPurge essentially does it).
                    #  filtered_matches.sort(key=lambda x: x[2], reverse=True)
                    filtered_matches.sort(key=lambda x: x[3])
                    for match_args in filtered_matches:
                        match = _match(*match_args[0:3])
                        if match:
                            return match

            return None


def compare_suffixes(
    suffix_ref: str,
    suffix_query: str,
    wildcard_ref: bool = False,
    wildcard_query: bool = False,
) -> MatchTuple:
    """
    Find out whether one string is the suffix of the other one, allowing mismatches.
    Used to find an anchored 3' adapter when no indels are allowed.

    Args:
        suffix_ref, suffix_query: The suffices to compare.
        wildcard_ref, wildcard_query: Whether wildcards are valid in either of
            the suffices.
    """
    suffix_ref = suffix_ref[::-1]
    suffix_query = suffix_query[::-1]
    _, length, _, _, matches, errors = compare_prefixes(
        suffix_ref, suffix_query, wildcard_ref, wildcard_query
    )
    return (
        len(suffix_ref) - length,
        len(suffix_ref),
        len(suffix_query) - length,
        len(suffix_query),
        matches,
        errors,
    )
