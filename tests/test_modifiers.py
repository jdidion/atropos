# coding: utf-8
from atropos.adapters import *
from atropos.align import MatchInfo
from atropos.commands.trim.modifiers import *
from atropos.io.seqio import Sequence
from atropos.util import reverse_complement as rc

DUMMY_ADAPTER = Adapter("ACGT", FRONT)


def front_match(read):
    match = Match(0, 2, 0, 2, 1, 0, True, DUMMY_ADAPTER, read)
    match_info = MatchInfo(
        "read",
        0,
        0,
        2,
        "",
        "AC",
        "ACGTAC",
        "adapter",
        "",
        "##",
        "######",
        True,
        2,
        2,
        2,
    )
    return (match, [match_info])


def back_match(read):
    match = Match(6, 8, 6, 8, 1, 0, False, DUMMY_ADAPTER, read)
    match_info = MatchInfo(
        "read",
        0,
        6,
        8,
        "ACGTAC",
        "GT",
        "",
        "adapter",
        "######",
        "##",
        "",
        False,
        2,
        2,
        2,
    )
    return (match, [match_info])


def test_unconditional_cutter():
    uc = UnconditionalCutter(lengths=[5])
    s = Sequence("read1", 'abcdefg')
    assert UnconditionalCutter(lengths=[2])(s).sequence == 'cdefg'
    assert UnconditionalCutter(lengths=[-2])(s).sequence == 'abcde'
    assert UnconditionalCutter(lengths=[100])(s).sequence == ''
    assert UnconditionalCutter(lengths=[-100])(s).sequence == ''


def test_nend_trimmer():
    trimmer = NEndTrimmer()
    seqs = ['NNNNAAACCTTGGNNN', 'NNNNAAACNNNCTTGGNNN', 'NNNNNN']
    trims = ['AAACCTTGG', 'AAACNNNCTTGG', '']
    for seq, trimmed in zip(seqs, trims):
        _seq = Sequence('read1', seq, qualities='#' * len(seq))
        _trimmed = Sequence('read1', trimmed, qualities='#' * len(trimmed))
        assert trimmer(_seq) == _trimmed


def test_quality_trimmer():
    read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
    qt = QualityTrimmer(10, 10, 33)
    assert qt(read) == Sequence('read1', 'GTTTAC', '456789')
    qt = QualityTrimmer(0, 10, 33)
    assert qt(read) == Sequence('read1', 'ACGTTTAC', '##456789')
    qt = QualityTrimmer(10, 0, 33)
    assert qt(read) == Sequence('read1', 'GTTTACGTA', '456789###')


def test_Modifiers_single():
    m = SingleEndModifiers()
    m.add_modifier(UnconditionalCutter, lengths=[5])
    mod1 = m.get_modifiers(read=1)
    mod2 = m.get_modifiers(read=2)
    assert len(mod1) == 1
    assert isinstance(mod1[0], UnconditionalCutter)
    assert len(mod2) == 0
    # test single-end
    read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
    mod_read = m.modify(read)
    assert len(mod_read) == 1
    assert mod_read[0].sequence == 'TACGTA'


def test_Modifiers_paired_legacy():
    m = PairedEndModifiers(paired="first")
    m.add_modifier(UnconditionalCutter, lengths=[5])
    read1 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
    read2 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
    mod_read1, mod_read2 = m.modify(read1, read2)
    assert mod_read1.sequence == 'TACGTA'
    assert mod_read2.sequence == 'ACGTTTACGTA'


def test_Modifiers_paired_both():
    m = PairedEndModifiers(paired="both")
    m.add_modifier(UnconditionalCutter, read=1 | 2, lengths=[5])
    mod1 = m.get_modifiers(read=1)
    mod2 = m.get_modifiers(read=2)
    assert len(mod1) == 1
    assert len(mod2) == 1
    assert isinstance(mod1[0], UnconditionalCutter)
    assert isinstance(mod2[0], UnconditionalCutter)
    read1 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
    read2 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
    mod_read1, mod_read2 = m.modify(read1, read2)
    assert mod_read1.sequence == 'TACGTA'
    assert mod_read2.sequence == 'TACGTA'


def test_min_cutter_T_T():
    unconditional_before = UnconditionalCutter((2, -2))
    unconditional_after = UnconditionalCutter((1, -1))
    min_trimmer = MinCutter((5, -5), True, True)
    read1 = Sequence('read1', "CAATCGATCGAACGTACCGAT")
    assert read1.clipped == [0, 0, 0, 0], str(read1.clipped)
    read1 = unconditional_before(read1)
    assert read1.sequence == "ATCGATCGAACGTACCG"
    assert read1.clipped == [2, 2, 0, 0], str(read1.clipped)
    # test without adapter trimming
    assert min_trimmer(read1).sequence == "ATCGATCGAACGTACCG"
    # test with adapter trimming
    read2 = read1[:]
    read2.sequence = "ATCGAACGTACCG"
    read2.match, read2.match_info = front_match(read2)
    read3 = min_trimmer(read2)
    assert read3.sequence == "TCGAACGTACCG", read3.sequence
    assert read3.clipped == [2, 2, 1, 0]
    # test with subsequent clipping
    read4 = unconditional_after(read2)
    assert read4.sequence == "TCGAACGTACC", read4.sequence
    assert read4.clipped == [2, 2, 1, 1], read4.clipped
    read5 = min_trimmer(read4)
    assert read5.sequence == "TCGAACGTACC", read5.sequence
    assert read5.clipped == [2, 2, 1, 1], read5.clipped


def test_min_cutter_F_T():
    unconditional_before = UnconditionalCutter((2, -2))
    unconditional_after = UnconditionalCutter((1, -1))
    min_trimmer = MinCutter((5, -5), False, True)
    read1 = Sequence('read1', "CAATCGATCGAACGTACCGAT")
    read1 = unconditional_before(read1)
    assert read1.sequence == "ATCGATCGAACGTACCG"
    assert read1.clipped == [2, 2, 0, 0]
    # test without adapter trimming
    assert min_trimmer(read1).sequence == "ATCGATCGAACGTACCG"
    # test with adapter trimming
    read2 = read1[:]
    read2.match, read2.match_info = front_match(read2)
    read2.sequence = "CGATCGAACGTACCG"
    read3 = min_trimmer(read2)
    assert read3.sequence == "GAACGTACCG", read3.sequence
    assert read3.clipped == [2, 2, 5, 0]
    # test with subsequent clipping
    read4 = unconditional_after(read2)
    assert read4.sequence == "GATCGAACGTACC"
    assert read4.clipped == [2, 2, 1, 1]
    read5 = min_trimmer(read4)
    assert read5.sequence == "GAACGTACC", read5.sequence
    assert read5.clipped == [2, 2, 5, 1]


def test_min_cutter_T_F():
    unconditional_before = UnconditionalCutter((2, -2))
    min_trimmer = MinCutter((4, -4), True, False)
    read1 = Sequence('read1', "CAATCGATCGAACGTACCGAT")
    read1 = unconditional_before(read1)
    assert read1.sequence == "ATCGATCGAACGTACCG"
    assert read1.clipped == [2, 2, 0, 0]
    # test without adapter trimming
    assert min_trimmer(read1).sequence == "CGATCGAACGTAC"


def test_non_directional_bisulfite_trimmer():
    trimmer = NonDirectionalBisulfiteTrimmer(rrbs=True)
    read1 = Sequence('read1', "CAATCGATCGA")
    read2 = Sequence('read2', "CTATCGATC")
    read2.match, read2.match_info = back_match(read2)
    read3 = Sequence('read3', "CTATCGATCCA")
    # assert trimmer(read1).sequence == "ATCGATC"
    assert trimmer(read2).sequence == "CTATCGA"
    assert trimmer(read3).sequence == "CTATCGATCCA"


def test_TruSeq_trimmer():
    trimmer = TruSeqBisulfiteTrimmer()
    read1 = Sequence('read1', "CTATCGATCCACGAGACTAAC")
    assert trimmer(read1).sequence == "ATCCACGAGACTAAC"


def test_Swift_trimmer():
    trimmer = SwiftBisulfiteTrimmer()
    seq = "".join(["ACGT"] * 30)
    read1 = Sequence('read1', seq)
    read2 = Sequence('read2', seq)
    trimmed = trimmer(read1, read2)
    assert trimmed[0].sequence == seq[:-10]
    assert trimmed[1].sequence == seq[10:]


def test_overlapping():
    trimmer = MergeOverlapping(min_overlap=10, error_rate=0.1)
    a1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC'
    a2 = reverse_complement('AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA')
    frag = 'CCAAGCAGACATTCACTCAGATTGCA'
    r1 = (frag + a1)[0:40]
    q1 = '#' * 40
    r2 = reverse_complement(a2 + frag)[0:40]
    q2 = '!' * 40
    parser = AdapterParser()
    adapter1 = parser.parse_from_spec(a1)
    adapter2 = parser.parse_from_spec(a2)
    cutter = AdapterCutter([adapter1, adapter2])
    read1 = Sequence('foo', r1, q1)
    read1 = cutter(read1)
    assert len(read1) == 26
    read2 = Sequence('foo', r2, q2)
    read2 = cutter(read2)
    assert len(read2) == 26
    # complete overlap
    read1_merged, read2_merged = trimmer(read1, read2)
    assert read1_merged.merged
    assert read2_merged is None
    assert read1 == read1_merged
    # partial overlap
    read1.merged = False
    read2 = read2.subseq(0, 24)[2]
    read1_merged, read2_merged = trimmer(read1, read2)
    assert read1_merged.merged
    assert read2_merged is None
    assert read1 == read1_merged
    # partial overlap r1, r2
    read1.merged = False
    read1 = read1.subseq(0, 24)[2]
    read1_merged, read2_merged = trimmer(read1, read2)
    assert read1_merged.merged
    assert read2_merged is None
    assert len(read1_merged) == 26
    assert read1_merged.sequence == 'CCAAGCAGACATTCACTCAGATTGCA'
    assert read1_merged.qualities == ('#' * 24) + ('!' * 2)
    # errors
    # round(0.1 * 24) = 2, so 2 errors should pass but 3 should not
    read1.merged = False
    r1_seq = list(read1.sequence)
    r1_seq[10] = reverse_complement(r1_seq[10])
    r1_seq[20] = reverse_complement(r1_seq[20])
    read1.sequence = "".join(r1_seq)
    read1_merged, read2_merged = trimmer(read1, read2)
    assert read1_merged.merged
    assert read2_merged is None
    assert len(read1_merged) == 26
    assert read1_merged.sequence == 'CCAAGCAGACTTTCACTCAGTTTGCA'
    assert read1_merged.qualities == ('#' * 24) + ('!' * 2)
    # too few overlapping bases
    read1.merged = False
    r1_seq[15] = reverse_complement(r1_seq[15])
    read1.sequence = "".join(r1_seq)
    read1_merged, read2_merged = trimmer(read1, read2)
    assert read1_merged.merged is False
    assert read2 is not None


def test_overlapping_with_error_correction():
    trimmer = MergeOverlapping(
        min_overlap=10, error_rate=0.1, mismatch_action='liberal'
    )
    r1 = 'AGATCGGAAGACCGTCATGTAGGGAAAGAGTGTAGATCTC'
    q1 = 'FFFFFFFFFFF#FFFFFFFFFFFFFFFFFFFFF#######'
    r2 = reverse_complement('AGATCGGTAGAGCGTCGTGTAGGGAAATAGTGTAGATCTC')
    q2 = ''.join(reversed('FFFFFFFFFFFFFFFF#FFFFFFFFFF#FFFFFFFFFFFF'))
    read1 = Sequence('foo', r1, q1)
    read2 = Sequence('foo', r2, q2)
    read1_merged, read2_merged = trimmer(read1, read2)
    assert read1_merged.merged
    assert read2_merged is None
    assert read1_merged.sequence == 'AGATCGGTAGAGCGTCATGTAGGGAAAGAGTGTAGATCTC'
    assert read1_merged.qualities == 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF#######'


def test_mismatched_adapter_overlaps():
    """
    This is a test case from real data. The adapter overlaps 1 less bp
    on the fw read than on the reverse read. We want to make sure that
    the extra 'A' base gets trimmed.
    adapter                                                                                                                                                          GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
    actual                                                               TTGTTTTTATGGAGAGAGTTTTAAGGTTTATTTTAGTTTTAAAGGATATTGTAGGTTAGAGGGAAAGTGTATGATGAAGGTATATATTGGTAGATCGGAAGAGCACACGTCTGAACTTCAGTCAC
    actual rc                          TATGTTCTTTCCCTTCACGTCTCTCTTCGGATCTTTATTGTGATGAGTTGAAAATAAAGGTTAAGTATAGATAAAAAAGTTATTATAGTTTAGAGGGTAAGTGTATGATGGAGTAAAATATTGGT
    adapter rc AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
    """
    r1 = 'TTGTTTTTATGGAGAGAGTTTTAAGGTTTATTTTAGTTTTAAAGGATATTGTAGGTTAGAGGGAAAGTGTATGATGAAGGTATATATTGGTAGATCGGAAGAGCACACGTCTGAACTTCAGTCAC'
    r2 = 'ACCAATATTTTACTCCATCATACACTTACCCTCTAAACTATAATAACTTTTTTATCTATACTTAACCTTTATTTTCAACTCATCACAATAAAGATCCGAAGAGAGACGTGAAGGGAAAGAACATA'
    a1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"  # TruSeq index 7
    a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"  # TruSeq universal
    parser = AdapterParser()
    adapter1 = parser.parse_from_spec(a1)
    adapter2 = parser.parse_from_spec(a2)
    # the data has a fairly high error rate
    cutter = InsertAdapterCutter(
        adapter1, adapter2, max_insert_mismatch_frac=0.3, max_adapter_mismatch_frac=0.3
    )
    read1 = Sequence('foo', r1, '#' * 125)
    read2 = Sequence('foo', r2, '#' * 125)
    new_read1, new_read2 = cutter(read1, read2)
    assert (len(new_read1)) == 91
    assert (len(new_read2)) == 91
    assert (
        new_read1.sequence ==
        'TTGTTTTTATGGAGAGAGTTTTAAGGTTTATTTTAGTTTTAAAGGATATTGTAGGTTAGAGGGAAAGTGTATGATGAAGGTATATATTGGT'
    )


def test_error_correction():
    a1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC'
    a2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA'
    frag = 'CCAAGCAGACATTCACTCAGATTGCA'
    correct_frag = 'CCAAGTAGACATTCGCTCAGATTGCA'
    r1 = list(frag)
    # C>T at pos 6
    r1[5] = 'T'
    q1 = ['#'] * 40
    # quality of read1 > quality of read2 at pos 6
    q1[5] = 'A'
    r1 = (''.join(r1) + a1)[0:40]
    q1 = ''.join(q1)
    r2 = list(frag)
    # A>G at pos 15
    r2[14] = 'G'
    q2 = ['#'] * 40
    # quality of read2 > quality of read1 at pos 11
    q2[len(frag) - 15] = 'A'
    r2 = reverse_complement(reverse_complement(a2) + ''.join(r2))[0:40]
    q2 = ''.join(q2)
    read1 = Sequence('foo', r1, q1)
    read2 = Sequence('foo', r2, q2)
    parser = AdapterParser()
    adapter1 = parser.parse_from_spec(a1)
    adapter2 = parser.parse_from_spec(a2)
    cutter = InsertAdapterCutter(adapter1, adapter2, mismatch_action='liberal')
    new_read1, new_read2 = cutter(read1, read2)
    assert len(new_read1) == 26
    assert new_read1.insert_overlap
    assert new_read1.sequence == correct_frag
    assert len(new_read2) == 26
    assert new_read2.insert_overlap
    assert new_read2.sequence == reverse_complement(correct_frag)


def test_error_correction_no_insert_match_one_adapter_match():
    a1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC'
    a2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA'
    a2_mod = 'ACATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA'
    frag = 'CCAAGCAGACATTCACTCAGATTGCA'
    correct_frag = 'CCAAGTAGACATTCGCTCAGATTGCA'
    r1 = list(frag)
    # C>T at pos 6
    r1[5] = 'T'
    q1 = ['#'] * 40
    # quality of read1 > quality of read2 at pos 6
    q1[5] = 'A'
    r1 = (''.join(r1) + a1)[0:40]
    q1 = ''.join(q1)
    r2 = list(frag)
    # A>G at pos 15
    r2[14] = 'G'
    q2 = ['#'] * 40
    # quality of read2 > quality of read1 at pos 11
    q2[len(frag) - 15] = 'A'
    r2 = reverse_complement(reverse_complement(a2_mod) + ''.join(r2))[0:40]
    q2 = ''.join(q2)
    read1 = Sequence('foo', r1, q1)
    read2 = Sequence('foo', r2, q2)
    parser1 = AdapterParser()
    adapter1 = parser1.parse_from_spec(a1)
    # Allow zero mismatches to prevent adapter alignment
    parser2 = AdapterParser(max_error_rate=0)
    adapter2 = parser2.parse_from_spec(a2)
    # Allow zero mismatches to prevent insert alignment
    cutter = InsertAdapterCutter(
        adapter1, adapter2, mismatch_action='liberal', max_insert_mismatch_frac=0
    )
    new_read1, new_read2 = cutter(read1, read2)
    assert len(new_read1) == 26
    assert not new_read1.insert_overlap
    assert new_read1.sequence == correct_frag
    assert len(new_read2) == 26
    assert not new_read2.insert_overlap
    assert new_read2.sequence == reverse_complement(correct_frag)


def test_error_correction_no_insert_match_two_adapter_matches():
    a1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC'
    a2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA'
    frag = 'CCAAGCAGACATTCACTCAGATTGCA'
    correct_frag = 'CCAAGTAGACATTCGCTCAGATTGCA'
    r1 = list(frag)
    # C>T at pos 6
    r1[5] = 'T'
    q1 = ['#'] * 40
    # quality of read1 > quality of read2 at pos 6
    q1[5] = 'A'
    r1 = (''.join(r1) + a1)[0:40]
    q1 = ''.join(q1)
    r2 = list(frag)
    # A>G at pos 15
    r2[14] = 'G'
    q2 = ['#'] * 40
    # quality of read2 > quality of read1 at pos 11
    q2[len(frag) - 15] = 'A'
    r2 = reverse_complement(reverse_complement(a2) + ''.join(r2))[0:40]
    q2 = ''.join(q2)
    read1 = Sequence('foo', r1, q1)
    read2 = Sequence('foo', r2, q2)
    parser = AdapterParser()
    adapter1 = parser.parse_from_spec(a1)
    adapter2 = parser.parse_from_spec(a2)
    # Allow zero mismatches to prevent insert alignment
    cutter = InsertAdapterCutter(
        adapter1, adapter2, mismatch_action='liberal', max_insert_mismatch_frac=0
    )
    new_read1, new_read2 = cutter(read1, read2)
    assert len(new_read1) == 26
    assert not new_read1.insert_overlap
    assert new_read1.sequence == correct_frag
    assert len(new_read2) == 26
    assert not new_read2.insert_overlap
    assert new_read2.sequence == reverse_complement(correct_frag)


def test_error_correction_unequal_read_lengths():
    # Test case for issue #51
    read1 = Sequence(
        'read1',
        'TTTGCAGCTTTTGTAGACAAGTGCTGTGCAGCTGATGTCAAAGAGACCTGCTTTGCTCTGGAGGGTCCAAAACTTGTAGCCTCAACCCGAGAAGCCATAGCCTAA',
        'CCCCCFCGGGGGBFFAFC<?BEADCCF<FFFFGFFDFDFFGGGGCFGGC?DFFFEC;,===??DG==DDDFFFFG8DDD7+5;;DF*=)))10885D**58>6=0',
    )
    read2 = Sequence(
        'read1',
        'ATAGGCTATGGCTTCTCGAGTTGAAGCTACAAGTTTTGGACCCTCCAGAGCAAAGCAGGTCTCTTTGACATCAGCTGCACAGCACTTGTCTACAAAAGCTGCAAAAGATCGGAAGAGCGTCTCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGACGTATCATTAAAAAAAAAAACACATCACATCAACAAGATAACACGACTTCTCCATCCACAGTACCGATGACCTCAACATTAGT',
        'CCCCCG@FCFGGCFGGGGFEFGFGGFCFGGGFGFGGGGGGGGGGGGGGGGGGGGGGGGGGG9FGGGGGGGFGDFFGGGGGGGGGGGGGGGGG8;>@?@FEGGGGGGGGGGGGGGGGGGGGG=DDFAEFFFGF>B>EA):DFFBDFFB6CDEDDD9=99DD>55)580:A5)*)*;DD>**51:0118):)4))1***0:*)*)((***0*.(((((*)/.)1/(6((()1.)(((6).-----8<:C<73',
    )
    aligner = InsertAligner(
        'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG',
        'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
    )
    insert_match, adapter_match1, adapter_match2 = aligner.match_insert(
        read1.sequence, read2.sequence
    )
    ec = ErrorCorrectorMixin('N')
    ec.correct_errors(read1, read2, insert_match, truncate_seqs=True)
    assert read1.corrected == 3
    assert read2.corrected == 3
    for i in (80, 86, 104):
        assert read1.sequence[i] == 'N', 'Read 1 not corrected to N at {}'.format(i)
        assert read2.sequence[104 - i] == 'N', 'Read 2 not corrected to N at {}'.format(
            104 - i
        )


def ints2quals(ints):
    return ''.join(chr(i + 33) for i in ints)


def test_overwrite_read():
    overwrite = OverwriteRead(20, 40, 10)
    lowseq = 'ACGT' * 5
    highseq = 'TCAG' * 5
    # mean lowq > 20, mean highq > 40
    lowq = (11, 31, 16, 24, 16, 20, 17, 19, 21, 28) * 2
    highq = (22, 62, 32, 48, 32, 40, 34, 38, 42, 56) * 2
    read1 = Sequence('foo', lowseq, ints2quals(lowq))
    read2 = Sequence('foo', highseq, ints2quals(highq))
    new_read1, new_read2 = overwrite(read1, read2)
    assert new_read1.sequence == lowseq
    assert new_read1.qualities == ints2quals(lowq)
    assert new_read2.sequence == highseq
    assert new_read2.qualities == ints2quals(highq)
    assert new_read1.corrected == new_read2.corrected == 0
    # mean lowq < 20, mean highq > 40
    lowq = tuple(i - 1 for i in lowq)
    read1 = Sequence('foo', lowseq, ints2quals(lowq))
    new_read1, new_read2 = overwrite(read1, read2)
    assert new_read1.sequence == rc(highseq)
    assert new_read1.qualities == ints2quals(reversed(highq))
    assert new_read2.sequence == highseq
    assert new_read2.qualities == ints2quals(highq)
    assert new_read1.corrected == new_read2.corrected == 1
    # mean lowq < 20, mean highq < 40
    highq = tuple(i - 1 for i in highq)
    read2 = Sequence('foo', highseq, ints2quals(highq))
    new_read1, new_read2 = overwrite(read1, read2)
    assert new_read1.sequence == lowseq
    assert new_read1.qualities == ints2quals(lowq)
    assert new_read2.sequence == highseq
    assert new_read2.qualities == ints2quals(highq)
    assert new_read1.corrected == new_read2.corrected == 0
