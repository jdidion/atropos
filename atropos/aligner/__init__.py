from enum import IntFlag
from typing import Tuple

from ._aligner import Aligner, MultiAligner, compare_prefixes, locate


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
