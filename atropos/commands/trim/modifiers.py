# coding: utf-8
"""This module implements all the read modifications that atropos supports.
A modifier must be callable. It is implemented as a function if no parameters
need to be stored, and as a class with a __call__ method if there are parameters
(or statistics).
"""
from collections import OrderedDict
import copy
import re
from atropos import AtroposError
from atropos.align import (
    Aligner, InsertAligner, SEMIGLOBAL, START_WITHIN_SEQ1, STOP_WITHIN_SEQ2)
from atropos.util import (
    BASE_COMPLEMENTS, reverse_complement, mean, quals2ints)
from .qualtrim import quality_trim_index, nextseq_trim_index

# Base classes

class Modifier(object):
    """Base clas for modifiers.
    """
    @property
    def name(self):
        """Modifier name.
        """
        return self.__class__.__name__

    @property
    def description(self):
        """Modifier description (for display).
        """
        return getattr(self, 'display_str', self.name)

    def summarize(self):
        """Returns a summary of the modifier's activity as a dict.
        """
        return {}

class ReadPairModifier(Modifier):
    """Base class of modifiers that edit a pair of reads simultaneously.
    """
    def __call__(self, read1, read2):
        raise NotImplementedError()

class Trimmer(Modifier):
    """Base class of modifiers that trim bases from reads.
    """
    def __init__(self):
        self.trimmed_bases = 0

    def __call__(self, read):
        raise NotImplementedError()

    def subseq(self, read, begin=0, end=None):
        """Returns a subsequence of a read.

        Args:
            read: The read to trim.
            begin: The first base of the subsequnce.
            end: The last base of the subsequence, or None for len(read).
        """
        if begin or (end is not None):
            front_bases, back_bases, new_read = read.subseq(begin, end)
            self.trimmed_bases += front_bases + back_bases
            return new_read
        else:
            return read

    def clip(self, read, front=0, back=0):
        """Returns a read with bases trimmed off the front/back.

        Args:
            read: The read to trim.
            front: The number of bases to trim from the front.
            back: The (negative) number of bases to trim from the back.
        """
        if (front or back) and len(read) > 0:
            front_bases, back_bases, new_read = read.clip(front, back)
            self.trimmed_bases += front_bases + back_bases
            return new_read
        else:
            return read

    def summarize(self):
        """Returns a summary dict.
        """
        return dict(bp_trimmed=self.trimmed_bases)

# Modifiers

class AdapterCutter(Modifier):
    """Repeatedly find one of multiple adapters in reads. The number of times
    the search is repeated is specified by the times parameter.

    Args:
        adapters: List of Adapter objects.
        times: Number of times to trim.
        action: What to do with a found adapter: None, 'trim', or 'mask'
    """
    def __init__(self, adapters=None, times=1, action='trim'):
        super(AdapterCutter, self).__init__()
        self.adapters = adapters or []
        self.times = times
        self.action = action
        self.with_adapters = 0

    def _best_match(self, read):
        """Find the best matching adapter in the given read.

        Returns:
            Either a Match instance or None if there are no matches.
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
        """Determine the adapter that best matches the given read.
        Since the best adapter is searched repeatedly, a list
        of Match instances is returned, which
        need to be applied consecutively to the read.
        The list is empty if there are no adapter matches.

        The read is converted to uppercase before it is compared to the adapter
        sequences.

        Cut found adapters from a single read. Return modified read.
        """
        if len(read) == 0:
            return read

        matches = []

        # try at most self.times times to remove an adapter
        trimmed_read = read
        for _ in range(self.times):
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
            assert len(trimmed_read) < len(read), \
                "Trimmed read isn't shorter than original"

        if self.action == 'trim':
            # read is already trimmed, nothing to do
            pass
        elif self.action == 'mask':
            # add N from last modification
            masked_sequence = trimmed_read.sequence
            for match in sorted(matches, reverse=True, key=lambda m: m.astart):
                nstr = 'N' * (
                    len(match.read.sequence) -
                    len(match.adapter.trimmed(match).sequence))
                # add N depending on match position
                if match.front:
                    masked_sequence = nstr + masked_sequence
                else:
                    masked_sequence += nstr
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

    def summarize(self):
        adapters_summary = OrderedDict()
        for adapter in self.adapters:
            adapters_summary[adapter.name] = adapter.summarize()
        return dict(
            records_with_adapters=self.with_adapters,
            adapters=adapters_summary)

# Other error correction approaches:
# https://www.ncbi.nlm.nih.gov/pubmed/25161220
# https://github.com/Malfoy/BCOOL

class ErrorCorrectorMixin(object):
    """Provides a method for error correction.

    Args:
        mismatch_action: The action to take when a mismatch between the
            overlapping portions of read1 and read2 is encountered. Valid
            values are 'liberal', 'conservative', 'N'.
        min_qual_difference: When mismatch_action=='conservative', the minimum
            difference in base quality between read1 and read2 required to
            perform the correction.
    """
    def __init__(self, mismatch_action=None, min_qual_difference=1):
        self.mismatch_action = mismatch_action
        self.r1r2_min_qual_difference = min_qual_difference
        self.r2r1_min_qual_difference = -1 * min_qual_difference
        self.corrected_pairs = 0
        self.corrected_bp = [0, 0]

    def correct_errors(self, read1, read2, insert_match, truncate_seqs=False):
        """Correct errors in overlapping reads.

        Args:
            read1: The read1 to correct.
            read2: The read2 to correct.
            insert_match: The match info telling where the reads overlap.
            truncate_seqs: Whether to truncate the sequences to equal size
                before correcting. This is necessary when the insert match is
                based on truncated sequences (e.g. when it was generated by
                InsertAligner).
        """
        # Do not attempt to correct an already corrected read
        if read1.corrected > 0 or read2.corrected > 0:
            return

        # read2 reverse-complement is the reference, read1 is the query
        r1_seq = list(read1.sequence)
        r2_seq = list(read2.sequence)
        len1 = len(r1_seq)
        len2 = len(r2_seq)

        has_quals = read1.qualities and read2.qualities
        if has_quals:
            r1_qual = list(read1.qualities)
            r2_qual = list(read2.qualities)
        elif self.mismatch_action in ('liberal', 'conservative'):
            raise ValueError(
                "Cannot perform quality-based error correction on reads "
                "lacking quality information")

        if truncate_seqs:
            if len1 > len2:
                r1_seq = r1_seq[:len2]
                if has_quals:
                    r1_qual = r1_qual[:len2]
            elif len2 > len1:
                r2_seq = r2_seq[:len1]
                if has_quals:
                    r2_qual = r2_qual[:len1]
                len2 = len1

        r1_start = insert_match[2]
        r1_end = insert_match[3]
        r1_changed = 0
        r2_start = len2 - insert_match[1]
        r2_end = len2 - insert_match[0]
        r2_changed = 0
        quals_equal = []

        for i, j in zip(
                range(r1_start, r1_end), range(r2_end - 1, r2_start - 1, -1)):
            base1 = r1_seq[i]
            base2 = BASE_COMPLEMENTS[r2_seq[j]]
            if base1 == base2:
                continue
            if self.mismatch_action == 'N':
                r1_seq[i] = 'N'
                r2_seq[j] = 'N'
                r1_changed += 1
                r2_changed += 1
            elif base1 == 'N':
                r1_seq[i] = base2
                if has_quals:
                    r1_qual[i] = r2_qual[j]
                r1_changed += 1
            elif base2 == 'N':
                r2_seq[j] = BASE_COMPLEMENTS[base1]
                if has_quals:
                    r2_qual[j] = r1_qual[i]
                r2_changed += 1
            elif has_quals:
                diff = ord(r1_qual[i]) - ord(r2_qual[j])
                if diff >= self.r1r2_min_qual_difference:
                    r2_seq[j] = BASE_COMPLEMENTS[base1]
                    r2_qual[j] = r1_qual[i]
                    r2_changed += 1
                elif diff <= self.r2r1_min_qual_difference:
                    r1_seq[i] = base2
                    r1_qual[i] = r2_qual[j]
                    r1_changed += 1
                elif self.mismatch_action == 'liberal':
                    quals_equal.append((i, j, base1, base2))

        if quals_equal:
            mean_qual1 = mean([ord(b) for b in r1_qual[r1_start:r1_end]])
            mean_qual2 = mean([ord(b) for b in r2_qual[r2_start:r2_end]])
            # Only make the corrections if one read is significantly better
            # than the other.
            # TODO: this method of determining whether one read is better
            # than the other is crude - come up with something better.
            diff = mean_qual1 - mean_qual2
            if diff > 1:
                # read1 is better than read2
                for i, j, base1, base2 in quals_equal:
                    r2_seq[j] = BASE_COMPLEMENTS[base1]
                    r2_qual[j] = r1_qual[i]
                    r2_changed += 1
            elif diff < -1:
                # read2 is better than read1
                for i, j, base1, base2 in quals_equal:
                    r1_seq[i] = base2
                    r1_qual[i] = r2_qual[j]
                    r1_changed += 1

        if r1_changed or r2_changed:
            self.corrected_pairs += 1

            def update_read(read, seq, qual, seq_len, read_num, num_changed):
                self.corrected_bp[read_num] += num_changed
                read.corrected = num_changed
                new_seq = ''.join(seq)
                partial = truncate_seqs and len(read.sequence) > seq_len
                if partial:
                    read.sequence = new_seq + read.sequence[seq_len:]
                else:
                    read.sequence = new_seq
                if has_quals:
                    new_qual = ''.join(qual)
                    if partial:
                        read.qualities = new_qual + read.qualities[seq_len:]
                    else:
                        read.qualities = new_qual

            if r1_changed:
                update_read(
                    read1, r1_seq, r1_qual if has_quals else None, len1, 0,
                    r1_changed)
            if r2_changed:
                update_read(
                    read2, r2_seq, r2_qual if has_quals else None, len2, 1,
                    r2_changed)

    def summarize(self):
        """Returns a summary dict.
        """
        return dict(
            records_corrected=self.corrected_pairs,
            bp_corrected=self.corrected_bp)

class InsertAdapterCutter(ReadPairModifier, ErrorCorrectorMixin):
    """AdapterCutter that uses InsertAligner to first try to identify
    insert overlap before falling back to semi-global adapter alignment.

    Args:
        adapter1, adapter2: Adapters.
        action: Action to take on adapter match: trim, mask (replace adapter
            with N's), lower (convert adapter bases to lower case),
            or None.
        mismatch_action: How to deal with mismatches. See
            :class:`ErrorCorrectorMixin`.
        symmetric: Whether to assume that the adapter should appear in the
            same place on overlapping reads.
        min_insert_overlap: Minimum overlap required between reads to be
            considered an insert match.
        aligner_args: Additional arguments to :class:`InsertAligner`.
    """
    def __init__(
            self, adapter1, adapter2, action='trim', mismatch_action=None,
            symmetric=True, min_insert_overlap=1, **aligner_args):
        ErrorCorrectorMixin.__init__(self, mismatch_action)
        self.adapter1 = adapter1
        self.adapter2 = adapter2
        self.aligner = InsertAligner(
            adapter1.sequence, adapter2.sequence,
            min_insert_overlap=min_insert_overlap,
            **aligner_args)
        self.min_insert_len = min_insert_overlap
        self.action = action
        self.symmetric = symmetric
        self.with_adapters = [0, 0]

    def __call__(self, read1, read2):
        read_lengths = [len(r) for r in (read1, read2)]
        if any(l < self.min_insert_len for l in read_lengths):
            return (read1, read2)

        match = self.aligner.match_insert(read1.sequence, read2.sequence)
        read1.insert_overlap = read2.insert_overlap = (match is not None)
        insert_match = None
        correct_errors = False

        if match:
            insert_match, adapter_match1, adapter_match2 = match
            correct_errors = self.mismatch_action is not None and insert_match[5] > 0
        else:
            adapter_match1 = self.adapter1.match_to(read1)
            adapter_match2 = self.adapter2.match_to(read2)
            # If the adapter matches are complementary, perform error correction
            if (
                    self.mismatch_action and adapter_match1 and
                    adapter_match2 and
                    adapter_match1.rstart == adapter_match2.rstart):
                insert_match = (
                    read_lengths[1] - adapter_match1.rstart,
                    read_lengths[1], 0, adapter_match1.rstart)
                correct_errors = True

        # If exactly one of the two alignments failed and symmetric is True,
        # duplicate the good alignment
        if self.symmetric and sum(
                bool(m) for m in (adapter_match1, adapter_match2)) == 1:

            def create_symmetric_match(match, read, read_len):
                if match.rstart > read_len:
                    return None
                match = match.copy()
                # If we're not dealing with equal-length reads, and this read
                # is shorter than the other, adjust the match end to be the
                # read length. The 'matches' and 'errors' attributes will be
                # wrong, but it shouldn't matter.
                if match.rstop < read_len:
                    match.astop -= (read_len - match.rstop)
                    match.rstop = read_len
                return match

            if adapter_match1:
                adapter_match2 = create_symmetric_match(adapter_match1, read2, read_lengths[1])
            else:
                adapter_match1 = create_symmetric_match(adapter_match2, read1, read_lengths[0])

            if self.mismatch_action and not insert_match and adapter_match1 and adapter_match2:
                # Assume that the symmetric read segments overlap and
                # perform error correction
                insert_match = (
                    read_lengths[1] - adapter_match1.rstart,
                    read_lengths[1], 0, adapter_match1.rstart)
                correct_errors = True

        if correct_errors:
            self.correct_errors(read1, read2, insert_match, truncate_seqs=True)

        return (
            self.trim(read1, self.adapter1, adapter_match1, 0),
            self.trim(read2, self.adapter2, adapter_match2, 1))

    def trim(self, read, adapter, match, read_idx):
        """Trim an adapter from a read.

        Args:
            read: The read to trim from.
            adapter: The Adapter to trim.
            match: The match details.
            read_idx: 0/1
        """
        if not match:
            read.match = None
            read.match_info = None
            return read

        match.adapter = adapter
        match.read = read
        match.front = False

        if self.action is None or match.rstart >= len(read):
            trimmed_read = read

        else:
            trimmed_read = adapter.trimmed(match)

            if self.action == 'mask':
                # add N from last modification
                masked_sequence = trimmed_read.sequence
                masked_sequence += 'N' * (len(read) - len(trimmed_read))
                # set masked sequence as sequence with original quality
                trimmed_read.sequence = masked_sequence
                trimmed_read.qualities = read.qualities
            elif self.action == 'lower':
                # TODO: offer option to mask with lower-case of trimmed base
                # This will happen as part of the refactoring to modify
                # Sequences in-place.
                pass

        trimmed_read.match = match
        trimmed_read.match_info = [match.get_info_record()]

        self.with_adapters[read_idx] += 1
        return trimmed_read

    def summarize(self):
        """Returns a summary dict.
        """
        adapters_summary = tuple(
            { adapter.name : adapter.summarize() }
            for adapter in (self.adapter1, self.adapter2))
        summary = dict(
            records_with_adapters=self.with_adapters,
            adapters=adapters_summary)
        if self.mismatch_action:
            summary.update(ErrorCorrectorMixin.summarize(self))
        return summary

class OverwriteRead(ReadPairModifier):
    """If one read is of significantly worse quality than the other, overwrite
    the poor quality read with the good-quality read.

    Computes the mean quality over the first `window_size` bases. If one read
    has mean quality < `worse_read_min_quality` and the other read has mean
    quality >= `better_read_min_quality`, replace the entire worse read with
    the better read.

    Args:
        worse_read_min_quality: The minimum quality below which a read is
            considered 'poor'
        better_read_min_quality: The quality above which a read is considered
            'good'
        window_size: The number of bases over which to compute mean quality
        base: The base for ascii-to-int quality conversion
        summary_fn: Function to summarize base quality (default=mean)
    """
    def __init__(
            self, worse_read_min_quality, better_read_min_quality,
            window_size, base=33, summary_fn=mean):
        self.worse_read_min_quality = worse_read_min_quality
        self.better_read_min_quality = better_read_min_quality
        self.window_size = window_size
        self.base = base
        self.summary_fn = summary_fn

    def __call__(self, read1, read2):
        if len(read1) < self.window_size or len(read2) < self.window_size:
            return (read1, read2)
        if not (read1.qualities and read2.qualities):
            raise ValueError(
                "OverwriteRead modifier does not work with reads "
                "lacking base qualities.")
        qual1 = list(quals2ints(read1.qualities[:self.window_size], self.base))
        summ1 = self.summary_fn(qual1)

        qual2 = list(quals2ints(read2.qualities[:self.window_size], self.base))
        summ2 = self.summary_fn(qual2)

        if (
                summ1 < self.worse_read_min_quality and
                summ2 >= self.better_read_min_quality):
            # TODO: not sure what the right value is here
            read2.corrected = 1
            read1 = read2.reverse_complement()
        elif (
                summ2 < self.worse_read_min_quality and
                summ1 >= self.better_read_min_quality):
            read1.corrected = 1
            read2 = read1.reverse_complement()

        return (read1, read2)

class UnconditionalCutter(Trimmer):
    """A modifier that unconditionally removes the first n or the last n bases
    from a read. Equivalent to MinCutter with count_trimmed=False and
    only_trimmed=False, but UnconditionalCutter is applied before adapter
    trimming.

    If the length is positive, the bases are removed from the beginning of the
    read. If the length is negative, the bases are removed from the end of the
    read.
    """
    display_str = "Cut unconditionally"

    def __init__(self, lengths=None):
        super().__init__()
        self.front_length = self.back_length = 0
        if lengths:
            self.front_length = sum(front for front in lengths if front > 0)
            self.back_length = sum(back for back in lengths if back < 0)

    def __call__(self, read):
        return self.clip(read, self.front_length, self.back_length)

class MinCutter(Trimmer):
    """Ensure that a minimum number of bases have been trimmed off each end.

    Args:
        lengths: Sequence of bases to trim. Numbers > 0 are summed to the
            total bases trimmed from the front of the read, and numbers < 0
            are summed to the total bases trimmed from the end of the read.
        count_trimmed: Whether to consider bases cut before or during adapter
            trimming when counting the number of bases that have already been
            cut.
        only_trimmed: Only cut read ends if they have already been
            adapter-trimmed.
    """
    display_str = "Cut conditionally"

    def __init__(self, lengths=None, count_trimmed=True, only_trimmed=False):
        super().__init__()
        self.front_length = self.back_length = 0
        if lengths:
            self.front_length = sum(front for front in lengths if front > 0)
            self.back_length = sum(back for back in lengths if back < 0)
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
            """Returns number of bases that need to be trimmed.
            """
            if self.count_trimmed:
                trimmed = read.clipped[offset] + read.clipped[offset+2]
                if read.match:
                    trimmed += sum(
                        i.rsize_total
                        for i in read.match_info
                        if is_front == i.is_front)
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
            to_trim(1, False) if trim_back else 0)

class LengthTagModifier(Modifier):
    """Replace "length=..." strings in read names.
    """
    def __init__(self, length_tag="length="):
        self.regex = re.compile(r"\b" + length_tag + r"[0-9]*\b")
        self.length_tag = length_tag

    def __call__(self, read):
        read = read[:]
        if read.name.find(self.length_tag) >= 0:
            read.name = self.regex.sub(
                self.length_tag + str(len(read.sequence)), read.name)
        return read

class SuffixRemover(Modifier):
    """Remove a given suffix from read names.
    """
    def __init__(self, suffixes=None):
        self.suffixes = suffixes or []

    def __call__(self, read):
        name = read.name
        for suffix in self.suffixes:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
        read = read[:]
        read.name = name
        return read

class PrefixSuffixAdder(Modifier):
    """Add a suffix and a prefix to read names.
    """
    def __init__(self, prefix="", suffix=""):
        self.prefix = prefix
        self.suffix = suffix

    def __call__(self, read):
        read = read[:]
        adapter_name = 'no_adapter'
        if read.match is not None:
            adapter_name = read.match.adapter.name
        read.name = self.prefix.replace('{name}', adapter_name) + read.name + \
            self.suffix.replace('{name}', adapter_name)
        return read

class DoubleEncoder(Modifier):
    """Double-encode colorspace reads, using characters ACGTN to represent
    colors.
    """
    def __init__(self):
        self.double_encode_trans = str.maketrans('0123.', 'ACGTN')

    def __call__(self, read):
        read = read[:]
        read.sequence = read.sequence.translate(self.double_encode_trans)
        return read

class ZeroCapper(Modifier):
    """Change negative quality values of a read to zero
    """
    def __init__(self, quality_base=33):
        qbase = quality_base
        self.zero_cap_trans = str.maketrans(
            ''.join(map(chr, range(qbase))), chr(qbase) * qbase)

    def __call__(self, read):
        read = read[:]
        read.qualities = read.qualities.translate(self.zero_cap_trans)
        return read

class PrimerTrimmer(Trimmer):
    """Trims primer base from colorspace reads.
    """
    display_str = "Primer-trimmed"

    def __call__(self, read):
        read = self.clip(read, 1)
        read.primer = ''
        return read

class NextseqQualityTrimmer(Trimmer):
    """NextSeq-specific quality trimmer.
    """
    display_str = "Quality trimmed (NextSeq)"

    def __init__(self, cutoff=0, base=33):
        super(NextseqQualityTrimmer, self).__init__()
        self.cutoff = cutoff
        self.base = base

    def __call__(self, read):
        if len(read) == 0:
            return read
        stop = nextseq_trim_index(read, self.cutoff, self.base)
        return self.subseq(read, end=stop)

class QualityTrimmer(Trimmer):
    """Trim bases from the start/end of reads based on their qualities.
    """
    display_str = "Quality-trimmed"

    def __init__(self, cutoff_front=0, cutoff_back=0, base=33):
        super(QualityTrimmer, self).__init__()
        self.cutoff_front = cutoff_front
        self.cutoff_back = cutoff_back
        self.base = base

    def __call__(self, read):
        if len(read) == 0:
            return read
        start, stop = quality_trim_index(
            read.qualities, self.cutoff_front, self.cutoff_back, self.base)
        return self.subseq(read, start, stop)

class NEndTrimmer(Trimmer):
    """Trims Ns from the 3' and 5' end of reads.
    """
    display_str = "End Ns trimmed"

    def __init__(self):
        super(NEndTrimmer, self).__init__()
        self.start_trim = re.compile(r'^N+')
        self.end_trim = re.compile(r'N+$')

    def __call__(self, read):
        if len(read) == 0:
            return read
        sequence = read.sequence
        start_cut = self.start_trim.match(sequence)
        end_cut = self.end_trim.search(sequence)
        start_cut = start_cut.end() if start_cut else 0
        end_cut = end_cut.start() if end_cut else len(read)
        return self.subseq(read, start_cut, end_cut)

class RRBSTrimmer(MinCutter):
    """Sequences that are adapter-trimmed are further trimmed 2 bp on the 3'
    end to remove potential methylation-biased bases from the end-repair
    reaction.
    """
    display_str = "RRBS-trimmed"

    def __init__(self, trim_5p=0, trim_3p=2):
        super().__init__(
            (trim_5p, -1 * trim_3p), count_trimmed=False, only_trimmed=True)

class NonDirectionalBisulfiteTrimmer(Modifier):
    """For non-directional RRBS/WGBS libraries (which implies that they were
    digested using MspI), sequences that start with either 'CAA' or 'CGA' will
    have 2 bp trimmed off the 5' end to remove potential methylation-biased
    bases from the end-repair reaction. Additionally, for RRBS reads, if CAA/CGA
    is not trimmed *and* the read has been adapter-trimmed, a minimum number of
    bases is trimmed from the 3' end.
    """
    display_str = "Bisulfite-trimmed (Non-directional)"
    _regex = re.compile(r"^C[AG]A")

    def __init__(self, trim_5p=2, trim_3p=2, rrbs=False):
        self._non_directional_cutter = MinCutter(
            [trim_5p], count_trimmed=False, only_trimmed=False)
        self.rrbs = rrbs
        if rrbs:
            self._rrbs_cutter = RRBSTrimmer(trim_3p)

    def __call__(self, read):
        if len(read) == 0:
            return read
        cutter = None
        if self._regex.match(read.sequence):
            cutter = self._non_directional_cutter
        elif self.rrbs:
            cutter = self._rrbs_cutter
        return cutter(read) if cutter else read

    def summarize(self):
        """Returns the number of trimmed bases.
        """
        bp_trimmed = (
            self._rrbs_cutter.trimmed_bases +
            self._non_directional_cutter.trimmed_bases)
        return dict(bp_trimmed=bp_trimmed)

class TruSeqBisulfiteTrimmer(MinCutter):
    """EpiGnome reads are trimmed at least 6 bp on the 5' end.
    """
    display_str = "Bisulfite-trimmed (EpiGnome/TruSeq)"

    def __init__(self):
        super().__init__((6,), count_trimmed=True, only_trimmed=False)

class SwiftBisulfiteTrimmer(ReadPairModifier):
    """For WGBS libraries prepared with the Swift Accel-NGS kit, 10 bp are
    trimmed  off the end of read1 and the beginning of read2.
    """
    display_str = "Bisulfite-trimmed (Swift)"

    def __init__(self, trim_5p1=0, trim_3p1=10, trim_5p2=10, trim_3p2=0):
        self._read1_cutter = MinCutter(
            (trim_5p1, -1 * trim_3p1), count_trimmed=False, only_trimmed=False)
        self._read2_cutter = MinCutter(
            (trim_5p2, -1 * trim_3p2), count_trimmed=False, only_trimmed=False)

    def __call__(self, read1, read2):
        return (self._read1_cutter(read1), self._read2_cutter(read2))

    def summarize(self):
        return dict(bp_trimmed=(
            self._read1_cutter.trimmed_bases,
            self._read2_cutter.trimmed_bases))

# TODO: InsertAdapterCutter should save the insert match, and
# MergeOverlapping should use that rather than doing another alignment

class MergeOverlapping(ReadPairModifier, ErrorCorrectorMixin):
    """Merge overlaping reads. The merged reads are stored in read1.
    """
    def __init__(self, min_overlap=0.9, error_rate=0.1, mismatch_action=None):
        ErrorCorrectorMixin.__init__(self, mismatch_action)
        self.min_overlap = int(min_overlap) if min_overlap > 1 else min_overlap
        self.error_rate = error_rate

    def __call__(self, read1, read2):
        len1 = len(read1.sequence)
        len2 = len(read2.sequence)
        min_overlap = self.min_overlap
        if min_overlap <= 1:
            min_overlap = max(2, round(self.min_overlap * min(len1, len2)))

        if len1 < min_overlap or len2 < min_overlap:
            return (read1, read2)

        insert_matched = read1.insert_overlap and read2.insert_overlap

        if insert_matched:
            # If we've already determined that there is an insert overlap
            # with a 3' overhang, we can constrain our alignment
            aflags = START_WITHIN_SEQ1 | STOP_WITHIN_SEQ2
        else:
            aflags = SEMIGLOBAL
        # align read1 to read2 reverse-complement to be compatible with
        # InsertAligner
        read2_rc = reverse_complement(read2.sequence)
        aligner = Aligner(read2_rc, self.error_rate, aflags)
        alignment = aligner.locate(read1.sequence)

        if alignment:
            r2_start, r2_stop, r1_start, r1_stop, matches, errors = alignment
            if matches >= min_overlap:
                # Only correct errors if we haven't already done correction in
                # the InsertAligner
                if self.mismatch_action and errors > 0 and not insert_matched:
                    self.correct_errors(read1, read2, alignment)

                if r2_start == 0 and r2_stop == len2:
                    # r2 is fully contained in r1
                    pass
                elif r1_start == 0 and r1_stop == len1:
                    # r1 is fully contained in r2
                    read1.sequence = read2_rc
                    read1.qualities = "".join(reversed(read2.qualities))
                elif r1_start > 0:
                    read1.sequence += read2_rc[r2_stop:]
                    if read1.qualities and read2.qualities:
                        read1.qualities += "".join(
                            reversed(read2.qualities))[r2_stop:]
                elif r2_start > 0:
                    read1.sequence = read2_rc + read1.sequence[r1_stop:]
                    if read1.qualities and read2.qualities:
                        read1.qualities = (
                            "".join(reversed(read2.qualities)) +
                            read1.qualities[r1_stop:])
                else:
                    raise AtroposError(
                        "Invalid alignment while trying to merge read "
                        "{}: {}".format(
                            read1.name, ",".join(str(i) for i in alignment)))

                read1.merged = True
                read2 = None

        return (read1, read2)

class Modifiers(object):
    """Base for classes that manage multiple modifiers.
    """
    def __init__(self):
        self.modifiers = []
        self.modifier_indexes = {}

    def add_modifier(self, mod_class, read=1|2, **kwargs):
        """Add a modifier of the specified type for one or both reads.

        Args:
            mod_class: The type of modifier to add.
            read: Which reads this modifier should modify; 1, 2, or 3 (1|2).
            kwargs: Additional keyword arguments to the modifier constructor.
        """
        raise NotImplementedError()

    def add_modifier_pair(self, mod_class, read1_args=None, read2_args=None):
        """Add a modifier for both reads in a pair.

        Args:
            mod_class: The type of modifier to add. Cannot be a subclass of
                ReadPairModifier.
            read1_args, read2_args: Additional keyword arguments to pass to the
                modifier constructors.
        """
        raise NotImplementedError()

    def _add_modifiers(self, mod_class, mods):
        idx = len(self.modifiers)
        self.modifiers.append(mods)
        if mod_class in self.modifier_indexes:
            self.modifier_indexes[mod_class].append(idx)
        else:
            self.modifier_indexes[mod_class] = [idx]
        return idx

    def has_modifier(self, mod_class):
        """Returns True if a modifier of the specified type has been added.
        """
        return mod_class in self.modifier_indexes

    def get_modifiers(self, mod_class=None, read=None):
        """Returns a list of modifiers that have been added.

        Args:
            mod_class: Restrict modifiers to those of a certain type.
            read: Return only the modifiers for a given read (1 or 2)
        """
        if mod_class is None:
            mods = copy.copy(self.modifiers)
        elif mod_class in self.modifier_indexes:
            mods = [self.modifiers[i] for i in self.modifier_indexes[mod_class]]
        else:
            mods = []

        if not (mods and read):
            return mods

        read_mods = []
        for mod in mods:
            if isinstance(mod, ReadPairModifier):
                read_mods.append(mod)
            elif mod[read-1] is not None:
                read_mods.append(mod[read-1])
        return read_mods

    def get_adapters(self):
        """Returns the adapters from the AdapterCutter or InsertAdapterCutter
        modifier, if any.

        Returns:
            A list [[read1_adapters], [read2_adapters]].
        """
        adapters = [[], []]
        if self.has_modifier(AdapterCutter):
            mod1, mod2 = self.get_modifiers(AdapterCutter)[0]
            if mod1:
                adapters[0] = mod1.adapters
            if mod2:
                adapters[1] = mod2.adapters
        elif self.has_modifier(InsertAdapterCutter):
            mod = self.get_modifiers(InsertAdapterCutter)[0]
            adapters[0] = [mod.adapter1]
            adapters[1] = [mod.adapter2]
        return adapters

    def modify(self, read1, read2=None):
        """Apply registered modifiers to a read/pair.

        Args:
            read1, read2: The reads to modify.

        Returns:
            A tuple of modified reads (read1, read2).
        """
        raise NotImplementedError()

    def summarize(self):
        """Returns a summary dict.
        """
        raise NotImplementedError()

class SingleEndModifiers(Modifiers):
    """Manages modifiers for single-end data.
    """
    def add_modifier(self, mod_class, read=1, **kwargs):
        if read != 1:
            raise ValueError("'read' must be 1 for single-end data")
        return self._add_modifiers(mod_class, [mod_class(**kwargs), None])

    def add_modifier_pair(self, mod_class, read1_args=None, read2_args=None):
        if read1_args is not None:
            return self.add_modifier(mod_class, **read1_args)

    def modify(self, read1, read2=None):
        for mods in self.modifiers:
            read1 = mods[0](read1)
        return (read1,)

    def summarize(self):
        summary = {}
        for mods in self.modifiers:
            mod = mods[0]
            summary[mod.name] = dict(
                (key, (value,))
                for key, value in mod.summarize().items())
            summary[mod.name]['desc'] = mod.description
        return summary

class PairedEndModifiers(Modifiers):
    """Manages modifiers for paired-end data.
    """
    def __init__(self, paired):
        super().__init__()
        self.paired = paired

    def add_modifier(self, mod_class, read=1|2, **kwargs):
        if issubclass(mod_class, ReadPairModifier):
            if self.paired != "both" and read == 1|2:
                raise ValueError(
                    "Must have paired-end reads to use modifer {}".format(
                        mod_class))
            mods = mod_class(**kwargs)
        else:
            mods = [None, None]
            if read & 1 > 0:
                mods[0] = mod_class(**kwargs)
            if read & 2 > 0 and self.paired == "both":
                mods[1] = mod_class(**kwargs)
            if not any(mods):
                return None
        return self._add_modifiers(mod_class, mods)

    def add_modifier_pair(self, mod_class, read1_args=None, read2_args=None):
        mods = [None, None]
        if read1_args is not None:
            mods[0] = mod_class(**read1_args)
        if read2_args is not None and self.paired == "both":
            mods[1] = mod_class(**read2_args)
        if any(mods):
            return self._add_modifiers(mod_class, mods)

    def modify(self, read1, read2=None):
        for mods in self.modifiers:
            if isinstance(mods, ReadPairModifier):
                read1, read2 = mods(read1, read2)
            else:
                if mods[0] is not None:
                    read1 = mods[0](read1)
                if mods[1] is not None:
                    read2 = mods[1](read2)
        return (read1, read2)

    def summarize(self):
        summary = {}
        for mods in self.modifiers:
            if isinstance(mods, ReadPairModifier):
                summary[mods.name] = mods.summarize()
                summary[mods.name]['desc'] = mods.description
            elif any(mods):
                name = desc = keys = None
                summ1 = summ2 = {}
                if mods[0]:
                    name = mods[0].name
                    desc = mods[0].description
                    summ1 = mods[0].summarize()
                    if summ1:
                        keys = summ1.keys()
                if mods[1]:
                    summ2 = mods[1].summarize()
                    if summ2:
                        if name:
                            assert name == mods[1].name
                            assert desc == mods[1].description
                            assert set(keys) == set(summ2.keys())
                        else:
                            name = mods[1].name
                            desc = mods[1].description
                            keys = summ2.keys()
                if keys:
                    summary[name] = dict(
                        (key, (summ1.get(key, None), summ2.get(key, None)))
                        for key in keys)
                    summary[name]['desc'] = desc
        return summary
