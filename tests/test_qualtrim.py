# coding: utf-8
from atropos.commands.trim.qualtrim import nextseq_trim_index
from atropos.io.seqio import Sequence


def test_nextseq_trim():
    s = Sequence('n', '', '')
    assert nextseq_trim_index(s, cutoff=22) == 0
    s = Sequence(
        'n',
        'TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG',
        'AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA',
    )
    assert nextseq_trim_index(s, cutoff=22) == 33
