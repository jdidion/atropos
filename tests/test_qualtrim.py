# coding: utf-8
from atropos.commands.trim.qualtrim import twocolor_trim_index
from atropos.io.sequence import Sequence


def test_twocolor_trim():
    s = Sequence("n", "", "")
    assert twocolor_trim_index(s, cutoff=22) == 0
    s = Sequence(
        "n",
        "TCTCGTATGCCGTCTTATGCTTGAAAAAAAAAAGGGGGGGGGGGGGGGGGNNNNNNNNNNNGGNGG",
        "AA//EAEE//A6///E//A//EA/EEEEEEAEA//EEEEEEEEEEEEEEE###########EE#EA",
    )
    assert twocolor_trim_index(s, cutoff=22) == 33
