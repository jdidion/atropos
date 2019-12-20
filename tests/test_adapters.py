from pytest import raises

from atropos.adapters import (
    Adapter,
    AdapterMatch,
    AdapterType,
    ColorspaceAdapter,
    InsertAligner,
    LinkedAdapter,
    parse_braces,
)
from atropos.io.sequence import Sequence

from .utils import approx_equal


def test_issue_52():
    adapter = Adapter(
        sequence="GAACTCCAGTCACNNNNN",
        where=AdapterType.BACK,
        max_error_rate=0.12,
        min_overlap=5,
        read_wildcards=False,
        adapter_wildcards=True,
    )
    read = Sequence(name="abc", sequence="CCCCAGAACTACAGTCCCGGC")
    am = AdapterMatch(
        astart=0,
        astop=17,
        rstart=5,
        rstop=21,
        matches=15,
        errors=2,
        front=None,
        adapter=adapter,
        read=read,
    )
    assert am.wildcards() == "GGC"
    """
    The result above should actually be 'CGGC' since the correct
    alignment is this one:

    adapter         GAACTCCAGTCACNNNNN
    mismatches           X     X
    read       CCCCAGAACTACAGTC-CCGGC

    Since we do not keep the alignment, guessing 'GGC' is the best we
    can currently do.
    """


def test_issue_80():
    # This issue turned out to not be an actual issue with the alignment
    # algorithm. The following alignment is found because it has more matches
    # than the 'obvious' one:
    #
    # TCGTATGCCGTCTTC
    # =========X==XX=
    # TCGTATGCCCTC--C
    #
    # This is correct, albeit a little surprising, since an alignment without
    # indels would have only two errors.
    adapter = Adapter(
        sequence="TCGTATGCCGTCTTC",
        where=AdapterType.BACK,
        max_error_rate=0.2,
        min_overlap=3,
        read_wildcards=False,
        adapter_wildcards=False,
    )
    read = Sequence(name="seq2", sequence="TCGTATGCCCTCC")
    result = adapter.match_to(read)
    assert read.original_length == 13, result
    assert result.errors == 3, result
    assert result.astart == 0, result
    assert result.astop == 15, result


def test_str():
    a = Adapter("ACGT", where=AdapterType.BACK, max_error_rate=0.1)
    str(a)
    str(a.match_to(Sequence(name="seq", sequence="TTACGT")))
    ca = ColorspaceAdapter("0123", where=AdapterType.BACK, max_error_rate=0.1)
    str(ca)


def test_color():
    with raises(ValueError):
        ColorspaceAdapter("0123", where=AdapterType.FRONT, max_error_rate=0.1)


def test_parse_braces():
    assert parse_braces("") == ""
    assert parse_braces("A") == "A"
    assert parse_braces("A{0}") == ""
    assert parse_braces("A{1}") == "A"
    assert parse_braces("A{2}") == "AA"
    assert parse_braces("A{2}C") == "AAC"
    assert parse_braces("ACGTN{3}TGACCC") == "ACGTNNNTGACCC"
    assert parse_braces("ACGTN{10}TGACCC") == "ACGTNNNNNNNNNNTGACCC"
    assert parse_braces("ACGTN{3}TGA{4}CCC") == "ACGTNNNTGAAAACCC"
    assert parse_braces("ACGTN{0}TGA{4}CCC") == "ACGTTGAAAACCC"


def test_parse_braces_fail():
    for expression in [
        "{",
        "}",
        "{}",
        "{5",
        "{1}",
        "A{-7}",
        "A{",
        "A{1",
        "N{7",
        "AN{7",
        "A{4{}",
        "A{4}{3}",
        "A{b}",
        "A{6X}",
        "A{X6}",
    ]:
        with raises(ValueError):
            parse_braces(expression)


def test_linked_adapter():
    linked_adapter = LinkedAdapter("AAAA", "TTTT")
    sequence = Sequence(name="seq", sequence="AAAACCCCCTTTT")
    match = linked_adapter.match_to(sequence)
    trimmed = linked_adapter.trimmed(match)
    assert trimmed.name == "seq"
    assert trimmed.sequence == "CCCCC"


def test_random_match_probabilities():
    a = Adapter("A", AdapterType.BACK)
    rmp = a.random_match_probabilities()
    assert rmp == [1.0, 0.25]
    a = Adapter("AC", AdapterType.BACK, gc_content=0.4)
    rmp = a.random_match_probabilities()
    assert rmp == [1.0, 0.3, 0.06]


def test_match_probability():
    a = InsertAligner("TTAGACATAT", "CAGTGGAGTA")
    k = 3
    n = 5
    i3 = (120 / (6 * 2)) * (0.25 ** 3) * (0.75 ** 2)
    i4 = (120 / 24) * (0.25 ** 4) * 0.75
    i5 = 0.25 ** 5
    assert approx_equal(a.match_probability(k, n), i3 + i4 + i5, 0.0001)


def test_insert_align():
    a1_seq = "TTAGACATATGG"
    a2_seq = "CAGTGGAGTATA"
    aligner = InsertAligner(a1_seq, a2_seq)
    r1 = "AGTCGAGCCCATTGCAGACT" + a1_seq[0:10]
    r2 = "AGTCTGCAATGGGCTCGACT" + a2_seq[0:10]
    insert_match, match1, match2 = aligner.match_insert(r1, r2)
    assert match1.rstart == 20
    assert match1.length == 10
    assert match2.rstart == 20
    assert match2.length == 10


def test_short_adapter_overlap():
    a1_seq = "TTAGACATAT"
    a2_seq = "CAGTGGAGTA"
    seq1 = "GACAGGCCGTTTGAATGTTGACGGGATGTT"
    seq2 = "CATCCCGTCAACATTCAAACGGCCTGTCCA"
    aligner = InsertAligner(a1_seq, a2_seq)
    insert_match, match1, match2 = aligner.match_insert(seq1, seq2)
    assert match1.rstart == 28
    assert match1.length == 2
    assert match2.rstart == 28
    assert match2.length == 2


# TODO: implement this test
# def test_insert_align_different_lengths():
#     # Test for #51. This is the alignment:
#     # TTAGGCTATGGCTTCTCGGGTTGAGGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAA
#     #                                                                                                                       AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#     # ATAGGCTATGGCTTCTCGAGTTGAAGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAAAGATCGGAAGAGCGTCTCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGACGTATCATTAAAAAAAAAAACACATCACATCAACAAGATAACACGACTTCTCCATCCACAGTACCGATGACCTCAACATTAGT
#     #
#     # The inserts should align even though there is a gap at the 5' end of read1.
#     a1_seq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
#     a2_seq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
#     seq1 = "TTTGCAGCTTTTGTAGACAAGTGCTGTGCAGCTGATGTCAAAGAGACCTGCTTTGCTCTGGAGGGTCCAAAACTTGTAGCCTCAACCCGAGAAGCCATAGCCTAA"
#     seq2 = "ATAGGCTATGGCTTCTCGAGTTGAAGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAAAGATCGGAAGAGCGTCTCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGACGTATCATTAAAAAAAAAAACACATCACATCAACAAGATAACACGACTTCTCCATCCACAGTACCGATGACCTCAACATTAGT"
