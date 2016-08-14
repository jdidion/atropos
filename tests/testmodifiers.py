# coding: utf-8
from __future__ import print_function, division, absolute_import

from atropos.seqio import Sequence
from atropos.modifiers import *
from atropos.adapters import *
from atropos.align import MatchInfo

DUMMY_ADAPTER = Adapter("ACGT", FRONT)
def front_match(read):
	match = Match(0, 2, 0, 2, 1, 0, True, DUMMY_ADAPTER, read)
	match_info = MatchInfo("read", 0, 0, 2, "", "AC", "ACGTAC", "adapter", "", "##", "######", True, 2, 2, 2)
	return (match, [match_info])
def back_match(read):
	match = Match(6, 8, 6, 8, 1, 0, False, DUMMY_ADAPTER, read)
	match_info = MatchInfo("read", 0, 6, 8, "ACGTAC", "GT", "", "adapter", "######", "##", "", False, 2, 2, 2)
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
		_seq = Sequence('read1', seq, qualities='#'*len(seq))
		_trimmed = Sequence('read1', trimmed, qualities='#'*len(trimmed))
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
	m = Modifiers(paired=False)
	m.add_modifier(UnconditionalCutter, lengths=[5])
	mod1 = m.get_modifiers(read=1)
	mod2 = m.get_modifiers(read=2)
	assert len(mod1) == 1
	assert isinstance(mod1[0], UnconditionalCutter)
	assert len(mod2) == 0
	# test single-end
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	mod_read, mod_bp = m.modify(read)
	assert mod_read[0].sequence == 'TACGTA'

def test_Modifiers_paired_legacy():
	m = Modifiers(paired="first")
	m.add_modifier(UnconditionalCutter, lengths=[5])
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	read2 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	(mod_read, mod_read2), mod_bp = m.modify((read, read2))
	assert mod_read.sequence == 'TACGTA'
	assert mod_read2.sequence == 'ACGTTTACGTA'

def test_Modifiers_paired_both():
	m = Modifiers(paired="both")
	m.add_modifier(UnconditionalCutter, read=1|2, lengths=[5])
	mod1 = m.get_modifiers(read=1)
	mod2 = m.get_modifiers(read=2)
	assert len(mod1) == 1
	assert len(mod2) == 1
	assert isinstance(mod1[0], UnconditionalCutter)
	assert isinstance(mod2[0], UnconditionalCutter)
	read = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	read2 = Sequence('read1', 'ACGTTTACGTA', '##456789###')
	(mod_read, mod_read2), mod_bp = m.modify((read, read2))
	assert mod_read.sequence == 'TACGTA'
	assert mod_read2.sequence == 'TACGTA'

def test_min_cutter_T_T():
	unconditional_before = UnconditionalCutter((2,-2))
	unconditional_after = UnconditionalCutter((1,-1))
	min_trimmer = MinCutter((5,-5), True, True)
	
	read1 = Sequence('read1', "CAATCGATCGAACGTACCGAT")
	assert read1.clipped == [0,0,0,0], str(read1.clipped)
	read1 = unconditional_before(read1)
	assert read1.sequence == "ATCGATCGAACGTACCG"
	assert read1.clipped == [2,2,0,0], str(read1.clipped)
	
	# test without adapter trimming
	assert min_trimmer(read1).sequence == "ATCGATCGAACGTACCG"
	
	# test with adapter trimming
	read2 = read1[:]
	read2.sequence = "ATCGAACGTACCG"
	read2.match, read2.match_info = front_match(read2)
	read3 = min_trimmer(read2)
	assert read3.sequence == "TCGAACGTACCG", read3.sequence
	assert read3.clipped == [2,2,1,0]
	
	# test with subsequent clipping
	read4 = unconditional_after(read2)
	assert read4.sequence == "TCGAACGTACC", read4.sequence
	assert read4.clipped == [2,2,1,1], read4.clipped
	read5 = min_trimmer(read4)
	assert read5.sequence == "TCGAACGTACC", read5.sequence
	assert read5.clipped == [2,2,1,1], read5.clipped

def test_min_cutter_F_T():
	unconditional_before = UnconditionalCutter((2,-2))
	unconditional_after = UnconditionalCutter((1,-1))
	min_trimmer = MinCutter((5,-5), False, True)
	
	read1 = Sequence('read1', "CAATCGATCGAACGTACCGAT")
	read1 = unconditional_before(read1)
	assert read1.sequence == "ATCGATCGAACGTACCG"
	assert read1.clipped == [2,2,0,0]
	
	# test without adapter trimming
	assert min_trimmer(read1).sequence == "ATCGATCGAACGTACCG"
	
	# test with adapter trimming
	read2 = read1[:]
	read2.match, read2.match_info = front_match(read2)
	read2.sequence = "CGATCGAACGTACCG"
	read3 = min_trimmer(read2)
	assert read3.sequence == "GAACGTACCG", read3.sequence
	assert read3.clipped == [2,2,5,0]
	
	# test with subsequent clipping
	read4 = unconditional_after(read2)
	assert read4.sequence == "GATCGAACGTACC"
	assert read4.clipped == [2,2,1,1]
	read5 = min_trimmer(read4)
	assert read5.sequence == "GAACGTACC", read5.sequence
	assert read5.clipped == [2,2,5,1]

def test_min_cutter_T_F():
	unconditional_before = UnconditionalCutter((2,-2))
	min_trimmer = MinCutter((4,-4), True, False)
	
	read1 = Sequence('read1', "CAATCGATCGAACGTACCGAT")
	read1 = unconditional_before(read1)
	assert read1.sequence == "ATCGATCGAACGTACCG"
	assert read1.clipped == [2,2,0,0]
	
	# test without adapter trimming
	assert min_trimmer(read1).sequence == "CGATCGAACGTAC"
	
def test_non_directional_bisulfite_trimmer():
	trimmer = NonDirectionalBisulfiteTrimmer()
	read1 = Sequence('read1', "CAATCGATCGA")
	read2 = Sequence('read2', "CTATCGATC")
	read2.match, read2.match_info = back_match(read2)
	read3 = Sequence('read3', "CTATCGATCCA")
	#assert trimmer(read1).sequence == "ATCGATC"
	assert trimmer(read2).sequence == "CTATCGA"
	assert trimmer(read3).sequence == "CTATCGATCCA"

def test_EpiGnome_trimmer():
	trimmer = EpiGnomeBisulfiteTrimmer()
	read1 = Sequence('read1', "CTATCGATCCACGAGACTAAC")
	assert trimmer(read1).sequence == "ATCCACGAGACTAAC"

def test_Swift_trimmer():
	trimmer = SwiftBisulfiteTrimmer()
	short_read = Sequence('read1', "ACGTACGTACGTACGT")
	long_seq = "".join(["ACGT"] * 30)
	long_read = Sequence('read2', long_seq)
	assert trimmer(short_read).sequence == "GTACGT"
	assert trimmer(long_read).sequence == long_seq[10:-10]

def test_overlapping():
	trimmer = MergeOverlapping(min_overlap=10, error_rate = 0.1)
	a1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTC'
	a2 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTTA'
	frag = 'CCAAGCAGACATTCACTCAGATTGCA'
	r1 = (frag + a1)[0:40]
	q1 = '#' * 40
	r2 = reverse_complement(a2 + frag)[0:40]
	q2 = '!' * 40
	parser = AdapterParser()
	adapter1 = parser.parse(a1)
	adapter2 = parser.parse(reverse_complement(a2))
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
	assert read1_merged.qualities == ('#' * 24) + ( '!' * 2)
	
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
	assert read1_merged.qualities == ('#' * 24) + ( '!' * 2)
	
	# too few overlapping bases
	read1.merged = False
	r1_seq[15] = reverse_complement(r1_seq[15])
	read1.sequence = "".join(r1_seq)
	read1_merged, read2_merged = trimmer(read1, read2)
	assert read1_merged.merged is False
	assert read2 is not None
