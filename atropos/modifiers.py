# coding: utf-8
"""
This module implements all the read modifications that atropos supports.
A modifier must be callable. It is implemented as a function if no parameters
need to be stored, and as a class with a __call__ method if there are parameters
(or statistics).
"""
import re
from .qualtrim import quality_trim_index, nextseq_trim_index
from .align import Aligner, SEMIGLOBAL, InsertAligner
from .util import complement, reverse_complement
from collections import defaultdict
import copy
import logging

# Base classes

class ReadPairModifier(object):
    """Superclass of modifiers that edit a pair of reads simultaneously."""
    def __call__(self, read1, read2):
        raise NotImplemented()

class Trimmer(object):
    def __init__(self):
        self.trimmed_bases = 0
    
    def __call__(self, read):
        raise NotImplemented()
    
    def subseq(self, read, begin=0, end=None):
        if begin or (end is not None):
            front_bases, back_bases, new_read = read.subseq(begin, end)
            self.trimmed_bases += front_bases + back_bases
            return new_read
        else:
            return read
    
    def clip(self, read, front=0, back=0):
        if front or back:
            front_bases, back_bases, new_read = read.clip(front, back)
            self.trimmed_bases += front_bases + back_bases
            return new_read
        else:
            return read

# Modifiers

class AdapterCutter(object):
    """
    Repeatedly find one of multiple adapters in reads.
    The number of times the search is repeated is specified by the
    times parameter.
    """

    def __init__(self, adapters=[], times=1, action='trim'):
        """
        adapters -- list of Adapter objects

        action -- What to do with a found adapter: None, 'trim', or 'mask'
        """
        super(AdapterCutter, self).__init__()
        self.adapters = adapters
        self.times = times
        self.action = action
        self.with_adapters = 0

    def _best_match(self, read):
        """
        Find the best matching adapter in the given read.

        Return either a Match instance or None if there are no matches.
        """
        best = None
        for adapter in self.adapters:
            match = adapter.match_to(read)
            if match is None:
                continue

            # the no. of matches determines which adapter fits best
            if best is None or match.matches > best.matches:
                best = match
        return best

    def __call__(self, read):
        """
        Determine the adapter that best matches the given read.
        Since the best adapter is searched repeatedly, a list
        of Match instances is returned, which
        need to be applied consecutively to the read.
        The list is empty if there are no adapter matches.

        The read is converted to uppercase before it is compared to the adapter
        sequences.

        Cut found adapters from a single read. Return modified read.
        """
        matches = []

        # try at most self.times times to remove an adapter
        trimmed_read = read
        for t in range(self.times):
            match = self._best_match(trimmed_read)
            if match is None:
                # nothing found
                break
            matches.append(match)
            trimmed_read = match.adapter.trimmed(match)
        
        if not matches:
            trimmed_read.match = None
            trimmed_read.match_info = None
            return trimmed_read
        
        if __debug__:
            assert len(trimmed_read) < len(read), "Trimmed read isn't shorter than original"
        
        if self.action == 'trim':
            # read is already trimmed, nothing to do
            pass
        elif self.action == 'mask':
            # add N from last modification
            masked_sequence = trimmed_read.sequence
            for match in sorted(matches, reverse=True, key=lambda m: m.astart):
                ns = 'N' * (len(match.read.sequence) -
                            len(match.adapter.trimmed(match).sequence))
                # add N depending on match position
                if match.front:
                    masked_sequence = ns + masked_sequence
                else:
                    masked_sequence += ns
            # set masked sequence as sequence with original quality
            trimmed_read.sequence = masked_sequence
            trimmed_read.qualities = matches[0].read.qualities

            assert len(trimmed_read.sequence) == len(read)
        elif self.action is None:
            trimmed_read = read
        
        trimmed_read.match = matches[-1]
        trimmed_read.match_info = [match.get_info_record() for match in matches]
        
        self.with_adapters += 1
        return trimmed_read

class InsertAdapterCutter(ReadPairModifier):
    """
    AdapterCutter that uses InsertAligner to first try to identify
    insert overlap before falling back to semi-global adapter alignment.
    """
    def __init__(self, adapter1, adapter2, action='trim', mismatch_action=None, symmetric=True, **aligner_args):
        self.adapter1 = adapter1
        self.adapter2 = adapter2
        self.aligner = InsertAligner(adapter1.sequence, adapter2.sequence, **aligner_args)
        self.action = action
        self.mismatch_action = mismatch_action
        self.symmetric = symmetric
        self.with_adapters = [0, 0]
        self.corrected_pairs = 0
        self.corrected_bp = [0, 0]
    
    def __call__(self, read1, read2):
        insert_match, insert_match_len, adapter_match1, adapter_match2 = self.aligner.match_insert(
            read1.sequence, read2.sequence)
        
        if (self.mismatch_action and insert_match and insert_match[5] > 0 and
                (adapter_match1 or adapter_match2)):
            self.correct_errors(read1, read2, insert_match)
        
        if adapter_match1 is None and adapter_match2 is None:
            adapter_match1 = self.adapter1.match_to(read1)
            adapter_match2 = self.adapter2.match_to(read2)
        
        if self.symmetric:
            if adapter_match1 is None and adapter_match2:
                adapter_match1 = adapter_match2.copy()
            if adapter_match2 is None and adapter_match1:
                adapter_match2 = adapter_match1.copy()
            
        if adapter_match1 is not None:
            adapter_match1.adapter = self.adapter1
            adapter_match1.read = read1
            adapter_match1.front = False
        
        if adapter_match2 is not None:
            adapter_match2.adapter = self.adapter2
            adapter_match2.read = read2
            adapter_match2.front = False
        
        return (self.trim(read1, adapter_match1, 0), self.trim(read2, adapter_match2, 1))
    
    def trim(self, read, match, read_idx):
        # TODO: refactor so that both adapter cutters use the same function
        if not match:
            read.match = None
            read.match_info = None
            return read
        
        if match.rstart < len(read):
            trimmed_read = match.adapter.trimmed(match)
            
            if self.action == 'trim':
                # read is already trimmed, nothing to do
                pass
            elif self.action == 'mask':
                # add N from last modification
                masked_sequence = trimmed_read.sequence
                masked_sequence += ('N' * len(read) - len(trimmed_read))
                # set masked sequence as sequence with original quality
                trimmed_read.sequence = masked_sequence
                trimmed_read.qualities = read.qualities
            elif self.action is None:
                trimmed_read = read
        else:
            trimmed_read = read
        
        trimmed_read.match = match
        trimmed_read.match_info = [match.get_info_record()]
        
        self.with_adapters[read_idx] += 1
        return trimmed_read
    
    def correct_errors(self, read1, read2, insert_match):
        # read2 reverse-complement is the reference, read1 is the query
        r1_seq = list(read1.sequence)
        r2_seq = list(read2.sequence)
        l2 = len(r2_seq)
        
        has_quals = read1.qualities and read2.qualities
        if has_quals:
            r1_qual = list(read1.qualities)
            r2_qual = list(read2.qualities)
        elif self.mismatch_action == 'best':
            raise Exception("Cannot perform quality-based error correction on reads lacking quality information")
        
        r1_changed = 0
        r2_changed = 0
        
        for i, j in zip(
                range(insert_match[2], insert_match[3]),
                range(l2 - insert_match[0] - 1, l2 - insert_match[1] - 1, -1)):
            b1 = r1_seq[i]
            b2 = complement[r2_seq[j]]
            if b1 == b2:
                continue
            if b1 == 'N':
                r1_seq[i] = b2
                if has_quals:
                    r1_qual[i] = r2_qual[j]
                r1_changed += 1
            elif b2 == 'N':
                r2_seq[j] = complement[b1]
                if has_quals:
                    r2_qual[j] = r1_qual[i]
                r2_changed += 1
            elif self.mismatch_action == 'N':
                r1_seq[i] = 'N'
                r2_seq[i] = 'N'
                r1_changed += 1
                r2_changed += 1
            elif r1_qual[i] == r2_qual[j]:
                # default to read1 base
                r2_seq[j] = complement[b1]
                r2_changed += 1
            elif r1_qual[i] > r2_qual[j]:
                r2_seq[j] = complement[b1]
                r2_qual[j] = r1_qual[i]
                r2_changed += 1
            else:
                r1_seq[i] = b2
                r1_qual[i] = r2_qual[j]
                r1_changed += 1
        
        if r1_changed or r2_changed:
            self.corrected_pairs += 1
            if r1_changed:
                self.corrected_bp[0] += r1_changed
                read1.sequence = ''.join(r1_seq)
                read1.corrected = r1_changed
                if has_quals:
                    read1.qualities = ''.join(r1_qual)
            if r2_changed:
                self.corrected_bp[1] += r2_changed
                read2.sequence = ''.join(r2_seq)
                read2.corrected = r2_changed
                if has_quals:
                    read2.qualities = ''.join(r2_qual)

class UnconditionalCutter(Trimmer):
    """
    A modifier that unconditionally removes the first n or the last n bases from a read.
    Equivalent to MinCutter with count_trimmed=False and only_trimmed=False, but
    UnconditionalCutter is applied before adapter trimming.
    
    If the length is positive, the bases are removed from the beginning of the read.
    If the length is negative, the bases are removed from the end of the read.
    """
    
    display_str = "Cut unconditionally"
    
    def __init__(self, lengths=[]):
        super(UnconditionalCutter, self).__init__()
        self.front_length = sum(l for l in lengths if l > 0)
        self.back_length = sum(l for l in lengths if l < 0)

    def __call__(self, read):
        return self.clip(read, self.front_length, self.back_length)

class MinCutter(Trimmer):
    """
    Ensure that a minimum number of bases have been trimmed off each end.
    count_trimmed :: whether to consider bases cut before or during adapter trimming
        when counting the number of bases that have already been cut
    only_trimmed :: only cut read ends if they have already been adapter-trimmed.
    """
    
    display_str = "Cut conditionally"
    
    def __init__(self, lengths=[], count_trimmed=True, only_trimmed=False):
        super(MinCutter, self).__init__()
        self.front_length = sum(l for l in lengths if l > 0)
        self.back_length = sum(l for l in lengths if l < 0)
        self.count_trimmed = count_trimmed
        self.only_trimmed = only_trimmed
    
    def __call__(self, read):
        trim_front = trim_back = True
        if self.only_trimmed:
            if read.match:
                is_front = [match.is_front for match in read.match_info]
                if not any(is_front):
                    trim_front = False
                elif all(is_front):
                    trim_back = False
            else:
                return read
        
        # TODO: distinguish between adapter trimming and other trimming.
        # For things like Methyl-Seq, we want to trim additional bases
        # after adapter trimming, but we want to count other post-adapter
        # trimmed bases toward the minimum length
        def to_trim(offset, is_front):
            if self.count_trimmed:
                trimmed = read.clipped[offset] + read.clipped[offset+2]
                if read.match:
                    trimmed += sum(i.rsize_total for i in read.match_info if is_front == i.is_front)
            elif read.match:
                trimmed = read.clipped[offset+2]
            else:
                trimmed = read.clipped[offset]
            
            if is_front:
                return max(self.front_length - trimmed, 0)
            else:
                return min(trimmed + self.back_length, 0)
        
        return self.clip(
            read,
            to_trim(0, True) if trim_front else 0,
            to_trim(1, False) if trim_back else 0
        )

class LengthTagModifier(object):
    """
    Replace "length=..." strings in read names.
    """
    def __init__(self, length_tag="length="):
        self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
        self.length_tag = length_tag

    def __call__(self, read):
        read = read[:]
        if read.name.find(self.length_tag) >= 0:
            read.name = self.regex.sub(self.length_tag + str(len(read.sequence)), read.name)
        return read

class SuffixRemover(object):
    """
    Remove a given suffix from read names.
    """
    def __init__(self, suffixes=[]):
        self.suffixes = suffixes

    def __call__(self, read):
        name = read.name
        for s in self.suffixes:
            if name.endswith(s):
                name = name[:-len(s)]
        read = read[:]
        read.name = name
        return read

class PrefixSuffixAdder(object):
    """
    Add a suffix and a prefix to read names
    """
    def __init__(self, prefix="", suffix=""):
        self.prefix = prefix
        self.suffix = suffix

    def __call__(self, read):
        read = read[:]
        adapter_name = 'no_adapter' if read.match is None else read.match.adapter.name
        read.name = self.prefix.replace('{name}', adapter_name) + read.name + \
            self.suffix.replace('{name}', adapter_name)
        return read

class DoubleEncoder(object):
    """
    Double-encode colorspace reads, using characters ACGTN to represent colors.
    """
    def __init__(self):
        self.double_encode_trans = str.maketrans('0123.', 'ACGTN')

    def __call__(self, read):
        read = read[:]
        read.sequence = read.sequence.translate(self.double_encode_trans)
        return read

class ZeroCapper(object):
    """
    Change negative quality values of a read to zero
    """
    def __init__(self, quality_base=33):
        qb = quality_base
        self.zero_cap_trans = str.maketrans(''.join(map(chr, range(qb))), chr(qb) * qb)

    def __call__(self, read):
        read = read[:]
        read.qualities = read.qualities.translate(self.zero_cap_trans)
        return read

class PrimerTrimmer(Trimmer):
    display_str = "Primer-trimmed"
    
    def __call__(self, read):
        """Trim primer base from colorspace reads"""
        read = self.clip(read, 1)
        read.primer = ''
        return read

class NextseqQualityTrimmer(Trimmer):
    display_str = "Quality trimmed (NextSeq)"
    
    def __init__(self, cutoff=0, base=33):
        super(NextseqQualityTrimmer, self).__init__()
        self.cutoff = cutoff
        self.base = base
    
    def __call__(self, read):
        stop = nextseq_trim_index(read, self.cutoff, self.base)
        return self.subseq(read, end=stop)

class QualityTrimmer(Trimmer):
    display_str = "Quality-trimmed"
    
    def __init__(self, cutoff_front=0, cutoff_back=0, base=33):
        super(QualityTrimmer, self).__init__()
        self.cutoff_front = cutoff_front
        self.cutoff_back = cutoff_back
        self.base = base

    def __call__(self, read):
        start, stop = quality_trim_index(read.qualities, self.cutoff_front, self.cutoff_back, self.base)
        return self.subseq(read, start, stop)

class NEndTrimmer(Trimmer):
    """Trims Ns from the 3' and 5' end of reads"""
    
    display_str = "End Ns trimmed"
    
    def __init__(self):
        super(NEndTrimmer, self).__init__()
        self.start_trim = re.compile(r'^N+')
        self.end_trim = re.compile(r'N+$')

    def __call__(self, read):
        sequence = read.sequence
        start_cut = self.start_trim.match(sequence)
        end_cut = self.end_trim.search(sequence)
        start_cut = start_cut.end() if start_cut else 0
        end_cut = end_cut.start() if end_cut else len(read)
        return self.subseq(read, start_cut, end_cut)

class RRBSTrimmer(MinCutter):
    """
    Sequences that are adapter-trimmed are further trimmed 2 bp on the 3' end to
    remove potential methylation-biased bases from the end-repair reaction.
    """
    
    display_str = "RRBS-trimmed"
    
    def __init__(self, trim_5p=0, trim_3p=2):
        super(RRBSTrimmer, self).__init__((trim_5p, -1 * trim_3p), count_trimmed=False, only_trimmed=True)

class NonDirectionalBisulfiteTrimmer(object):
    """
    For non-directional RRBS/WGBS libraries (which implies that they were digested
    using MspI), sequences that start with either 'CAA' or 'CGA' will have 2 bp
    trimmed off the 5' end to remove potential methylation-biased bases from the
    end-repair reaction. Additionally, for RRBS reads, if CAA/CGA is not trimmed
    *and* the read has been adapter-trimmed, a minimum number of bases is trimmed
    from the 3' end.
    """
    
    display_str = "Bisulfite-trimmed (Non-directional)"
    _regex = re.compile("^C[AG]A")
    
    def __init__(self, trim_5p=2, trim_3p=2, rrbs=False):
        self._non_directional_cutter = MinCutter([trim_5p], count_trimmed=False, only_trimmed=False)
        self.rrbs = rrbs
        if rrbs:
            self._rrbs_cutter = RRBSTrimmer(trim_3p)
    
    def __call__(self, read):
        cutter = None
        if self._regex.match(read.sequence):
            cutter = self._non_directional_cutter
        elif self.rrbs:
            cutter = self._rrbs_cutter
        return cutter(read) if cutter else read
    
    @property
    def modified_bases(self):
        return self._rrbs_cutter.trimmed_bases + self._non_directional_cutter.trimmed_bases

class TruSeqBisulfiteTrimmer(MinCutter):
    """EpiGnome reads are trimmed at least 6 bp on the 5' end."""
    
    display_str = "Bisulfite-trimmed (EpiGnome/TruSeq)"
    
    def __init__(self):
        super(TruSeqBisulfiteTrimmer, self).__init__((6,), count_trimmed=True, only_trimmed=False)

class SwiftBisulfiteTrimmer(ReadPairModifier):
    """
    For WGBS libraries prepared with the Swift Accel-NGS kit, 10 bp are trimmed off the end
    of read1 and the beginning of read2.
    """
    
    display_str = "Bisulfite-trimmed (Swift)"
    
    def __init__(self, trim_5p1=0, trim_3p1=10, trim_5p2=10, trim_3p2=0):
        self._read1_cutter = MinCutter((trim_5p1, -1 * trim_3p1), count_trimmed=False, only_trimmed=False)
        self._read2_cutter = MinCutter((trim_5p2, -1 * trim_3p2), count_trimmed=False, only_trimmed=False)
    
    def __call__(self, read1, read2):
        return (self._read1_cutter(read1), self._read2_cutter(read2))

    @property
    def modified_bases(self):
        return self._read1_cutter.trimmed_bases + self._read2_cutter.trimmed_bases

# TODO: InsertAdapterCutter should save the insert match, and
# MergeOverlapping should use that rather than doing another alignment

class MergeOverlapping(ReadPairModifier):
    def __init__(self, min_overlap=0.9, error_rate=0.1):
        self.min_overlap = int(min_overlap) if min_overlap > 1 else min_overlap
        self.error_rate = error_rate
    
    def __call__(self, read1, read2):
        if read1.match and not read1.match.front and read2.match and not read2.match.front:
            min_overlap = self.min_overlap
            l1 = len(read1.sequence)
            l2 = len(read2.sequence)
            if min_overlap <= 1:
                min_overlap = max(2, round(self.min_overlap * min(l1, l2)))
                
            aligner = Aligner(read1.sequence, self.error_rate, SEMIGLOBAL)
            read2_rc = reverse_complement(read2.sequence)
            # TODO: Resolve overlap. We offer the user the choice to set mismatching bases to
            # the one with the highest quality or to 'N'. Additionally, set qualities of overlapping
            # bases to the highest value. The acutal alignment is not returned, so if there are
            # mismatches we don't know what they are, and figuring it out is complicated by indels.
            alignment = aligner.locate(read2_rc)
            if alignment:
                r1_start, r1_stop, r2_start, r2_stop, matches, errors = alignment
                if matches >= min_overlap:
                    if r2_start == 0 and r2_stop == l2:
                        # r2 is fully contained in r1
                        pass
                    elif r1_start == 0 and r1_stop == l1:
                        # r1 is fully contained in r2
                        read1.sequence = read2_rc
                        read1.qualities = "".join(reversed(read2.qualities))
                    elif r1_start > 0:
                        read1.sequence += read2_rc[r2_stop:]
                        if read1.qualities and read2.qualities:
                            read1.qualities += "".join(reversed(read2.qualities))[r2_stop:]
                    elif r2_start > 0:
                        read1.sequence = read2_rc + read1[r1_stop:]
                        if read1.qualities and read2.qualities:
                            read1.qualities = "".join(reversed(read2.qualities)) + read1.qualities[r1_stop:]
                    else:
                        raise Exception(
                            "Invalid alignment while trying to merge read {}: {}".format(
                                read1.name, ",".join(str(i) for i in alignment)))
                    
                    read1.merged = True
                    read2 = None
                
        return (read1, read2)

class Modifiers(object):
    def __init__(self, paired):
        self.modifiers = []
        self.modifier_indexes = defaultdict(lambda: [])
        self.paired = paired
    
    def add_modifier(self, mod_class, read=1|2, **kwargs):
        if issubclass(mod_class, ReadPairModifier):
            assert self.paired == "both" and read == 1|2
            mods = mod_class(**kwargs)
        else:
            mods = [None, None]
            if read & 1 > 0:
                mods[0] = mod_class(**kwargs)
            if read & 2 > 0 and self.paired == "both":
                mods[1] = mod_class(**kwargs)
            if all(m is None for m in mods):
                return None
        return self._add_modifiers(mod_class, mods)
    
    def add_modifier_pair(self, mod_class, read1_args=None, read2_args=None):
        mods = [None, None]
        if read1_args is not None:
            mods[0] = mod_class(**read1_args)
        if read2_args is not None and self.paired == "both":
            mods[1] = mod_class(**read2_args)
        if all(m is None for m in mods):
            return None
        return self._add_modifiers(mod_class, mods)
    
    def _add_modifiers(self, mod_class, mods):
        i = len(self.modifiers)
        self.modifiers.append(mods)
        self.modifier_indexes[mod_class].append(i)
        return i
    
    def get_modifiers(self, mod_class=None, read=None):
        if mod_class is None:
            mods = copy.copy(self.modifiers)
        else:
            mods = [self.modifiers[i] for i in self.modifier_indexes[mod_class]]
        if read is not None:
            mods = [m[read-1] for m in mods if m[read-1] is not None]
        return mods
    
    def get_trimmer_classes(self):
        return [
            klass for klass in self.modifier_indexes.keys()
            if issubclass(klass, Trimmer)
        ]
    
    def modify(self, record):
        bp = [0,0]
        if self.paired:
            read1, read2 = record
            bp[0] = len(read1.sequence)
            bp[1] = len(read2.sequence)
            for mods in self.modifiers:
                if isinstance(mods, ReadPairModifier):
                    read1, read2 = mods(read1, read2)
                else:
                    if mods[0] is not None:
                        read1 = mods[0](read1)
                    if mods[1] is not None:
                        read2 = mods[1](read2)
            reads = [read1, read2]
        else:
            read = record
            bp[0] = len(read.sequence)
            for mods in self.modifiers:
                read = mods[0](read)
            reads = [read]
        return (reads, bp)
