# coding: utf-8
import math
from .utils import approx_equal
from atropos.adapters import BACK
from atropos.align import (
    locate, compare_prefixes, compare_suffixes, Aligner, InsertAligner
)
from atropos.util import RandomMatchProbability


class TestAligner():

    def test(self):
        reference = 'CTCCAGCTTAGACATATC'
        aligner = Aligner(reference, 0.1, flags=BACK)
        aligner.locate('CC')

    def test_100_percent_error_rate(self):
        reference = 'GCTTAGACATATC'
        aligner = Aligner(reference, 1.0, flags=BACK)
        aligner.locate('CAA')


def test_polya():
    s = 'AAAAAAAAAAAAAAAAA'
    t = 'ACAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    result = locate(s, t, 0.0, BACK)
    # start_s, stop_s, start_t, stop_t, matches, cost = result
    assert result == (0, len(s), 4, 4 + len(s), len(s), 0)


# Sequences with IUPAC wildcards
# R=A|G, Y=C|T, S=G|C, W=A|T, K=G|T, M=A|C, B=C|G|T, D=A|G|T, H=A|C|T, V=A|C|G,
# N=A|C|G|T, X={}
WILDCARD_SEQUENCES = [
    'CCCATTGATC',  # original sequence without wildcards
    'CCCRTTRATC',  # R=A|G
    'YCCATYGATC',  # Y=C|T
    'CSSATTSATC',  # S=G|C
    'CCCWWWGATC',  # W=A|T
    'CCCATKKATC',  # K=G|T
    'CCMATTGMTC',  # M=A|C
    'BCCATTBABC',  # B=C|G|T
    'BCCATTBABC',  # B
    'CCCDTTDADC',  # D=A|G|T
    'CHCATHGATC',  # H=A|C|T
    'CVCVTTVATC',  # V=A|C|G
    'CCNATNGATC',  # N=A|C|G|T
    'CCCNTTNATC',
    # N
    #   'CCCXTTXATC',  # X
]


def test_compare_prefixes():
    assert compare_prefixes('AAXAA', 'AAAAATTTTTTTTT') == (0, 5, 0, 5, 4, 1)
    assert compare_prefixes('AANAA', 'AACAATTTTTTTTT', wildcard_ref=True) == (
        0, 5, 0, 5, 5, 0
    )
    assert compare_prefixes('AANAA', 'AACAATTTTTTTTT', wildcard_ref=True) == (
        0, 5, 0, 5, 5, 0
    )
    assert compare_prefixes('XAAAAA', 'AAAAATTTTTTTTT') == (0, 6, 0, 6, 4, 2)
    a = WILDCARD_SEQUENCES[0]
    for s in WILDCARD_SEQUENCES:
        r = s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
        result = compare_prefixes(a, r, wildcard_query=True)
        assert result == (0, 10, 0, 10, 10, 0), result
        result = compare_prefixes(r, a, wildcard_ref=True)
        assert result == (0, 10, 0, 10, 10, 0)
    for s in WILDCARD_SEQUENCES:
        for t in WILDCARD_SEQUENCES:
            r = s + 'GCCAGGG'
            result = compare_prefixes(s, r)
            assert result == (0, 10, 0, 10, 10, 0)
            result = compare_prefixes(r, s, wildcard_ref=True, wildcard_query=True)
            assert result == (0, 10, 0, 10, 10, 0)
    r = WILDCARD_SEQUENCES[0] + 'GCCAGG'
    for wildc_ref in (False, True):
        for wildc_query in (False, True):
            result = compare_prefixes(
                'CCCXTTXATC', r, wildcard_ref=wildc_ref, wildcard_query=wildc_query
            )
            assert result == (0, 10, 0, 10, 8, 2)


def test_compare_suffixes():
    assert compare_suffixes('AAXAA', 'TTTTTTTAAAAA') == (0, 5, 7, 12, 4, 1)
    assert compare_suffixes('AANAA', 'TTTTTTTAACAA', wildcard_ref=True) == (
        0, 5, 7, 12, 5, 0
    )
    assert compare_suffixes('AANAA', 'TTTTTTTAACAA', wildcard_ref=True) == (
        0, 5, 7, 12, 5, 0
    )
    assert compare_suffixes('AAAAAX', 'TTTTTTTAAAAA') == (0, 6, 6, 12, 4, 2)


def test_wildcards_in_adapter():
    r = 'CATCTGTCC' + WILDCARD_SEQUENCES[0] + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
    for a in WILDCARD_SEQUENCES:
        result = locate(a, r, 0.0, BACK, wildcard_ref=True)
        assert result == (0, 10, 9, 19, 10, 0), result
    a = 'CCCXTTXATC'
    result = locate(a, r, 0.0, BACK, wildcard_ref=True)
    assert result is None


def test_wildcards_in_read():
    a = WILDCARD_SEQUENCES[0]
    for s in WILDCARD_SEQUENCES:
        r = 'CATCTGTCC' + s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
        result = locate(a, r, 0.0, BACK, wildcard_query=True)
        if 'X' in s:
            assert result is None
        else:
            assert result == (0, 10, 9, 19, 10, 0), result


def test_wildcards_in_both():
    for a in WILDCARD_SEQUENCES:
        for s in WILDCARD_SEQUENCES:
            if 'X' in s or 'X' in a:
                continue

            r = 'CATCTGTCC' + s + 'GCCAGGGTTGATTCGGCTGATCTGGCCG'
            result = locate(a, r, 0.0, BACK, wildcard_ref=True, wildcard_query=True)
            assert result == (0, 10, 9, 19, 10, 0), result


def test_no_match():
    a = locate('CTGATCTGGCCG', 'AAAAGGG', 0.1, BACK)
    assert a is None, a


def test_factorial():
    f = RandomMatchProbability()
    # simple test
    assert f.factorial(0) == 1
    assert f.factorial(1) == 1
    assert f.factorial(3) == 6
    assert int(f.factorial(27)) == int(math.factorial(27))
    # test big number
    assert int(f.factorial(150)) == int(math.factorial(150))


def test_match_probability():
    a = InsertAligner('TTAGACATAT', 'CAGTGGAGTA')
    k = 3
    n = 5
    i3 = (120 / (6 * 2)) * (0.25 ** 3) * (0.75 ** 2)
    i4 = (120 / 24) * (0.25 ** 4) * 0.75
    i5 = 0.25 ** 5
    assert approx_equal(a.match_probability(k, n), i3 + i4 + i5, 0.0001)


def test_insert_align():
    a1_seq = 'TTAGACATATGG'
    a2_seq = 'CAGTGGAGTATA'
    aligner = InsertAligner(a1_seq, a2_seq)
    r1 = 'AGTCGAGCCCATTGCAGACT' + a1_seq[0:10]
    r2 = 'AGTCTGCAATGGGCTCGACT' + a2_seq[0:10]
    insert_match, match1, match2 = aligner.match_insert(r1, r2)
    assert match1.rstart == 20
    assert match1.length == 10
    assert match2.rstart == 20
    assert match2.length == 10


def test_short_adapter_overlap():
    a1_seq = 'TTAGACATAT'
    a2_seq = 'CAGTGGAGTA'
    seq1 = 'GACAGGCCGTTTGAATGTTGACGGGATGTT'
    seq2 = 'CATCCCGTCAACATTCAAACGGCCTGTCCA'
    aligner = InsertAligner(a1_seq, a2_seq)
    insert_match, match1, match2 = aligner.match_insert(seq1, seq2)
    assert match1.rstart == 28
    assert match1.length == 2
    assert match2.rstart == 28
    assert match2.length == 2


def test_insert_align_different_lengths():
    # Test for #51. This is the alignment:
    # TTAGGCTATGGCTTCTCGGGTTGAGGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAA
    #                                                                                                                       AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
    # ATAGGCTATGGCTTCTCGAGTTGAAGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAAAGATCGGAAGAGCGTCTCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGACGTATCATTAAAAAAAAAAACACATCACATCAACAAGATAACACGACTTCTCCATCCACAGTACCGATGACCTCAACATTAGT
    #
    # The inserts should align even though there is a gap at the 5' end of read1.
    a1_seq = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG'
    a2_seq = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    seq1 = 'TTTGCAGCTTTTGTAGACAAGTGCTGTGCAGCTGATGTCAAAGAGACCTGCTTTGCTCTGGAGGGTCCAAAACTTGTAGCCTCAACCCGAGAAGCCATAGCCTAA'
    seq2 = 'ATAGGCTATGGCTTCTCGAGTTGAAGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAAAGATCGGAAGAGCGTCTCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGACGTATCATTAAAAAAAAAAACACATCACATCAACAAGATAACACGACTTCTCCATCCACAGTACCGATGACCTCAACATTAGT'


def test_multi_aligner_no_mismatches():
    from atropos.align._align import MultiAligner

    a = MultiAligner(max_error_rate=0, min_overlap=3)
    matches = a.locate('AGAGATCAGATGACAGATC', 'GATCA')
    assert len(matches) == 2
    matches.sort(key=lambda x: x[4], reverse=True)
    assert matches[0][0] == 3
    assert matches[0][1] == 8
    assert matches[0][2] == 0
    assert matches[0][3] == 5
    assert matches[0][4] == 5
    assert matches[0][5] == 0
    assert matches[1][0] == 15
    assert matches[1][1] == 19
    assert matches[1][2] == 0
    assert matches[1][3] == 4
    assert matches[1][4] == 4
    assert matches[1][5] == 0


def test_multi_aligner_with_mismatches():
    from atropos.align._align import MultiAligner

    a = MultiAligner(max_error_rate=0.1, min_overlap=10)
    matches = a.locate('GATATCAGATGACAGATCAGAGATCAGAT', 'GAGATCAGATGA')
    assert len(matches) == 2
    matches.sort(key=lambda x: x[5])
    assert matches[0][0] == 19
    assert matches[0][1] == 29
    assert matches[0][2] == 0
    assert matches[0][3] == 10
    assert matches[0][4] == 10
    assert matches[0][5] == 0
    assert matches[1][0] == 0
    assert matches[1][1] == 12
    assert matches[1][2] == 0
    assert matches[1][3] == 12
    assert matches[1][4] == 11
    assert matches[1][5] == 1
