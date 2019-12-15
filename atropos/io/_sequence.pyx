# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
import copy

from atropos.errors import FormatError
from atropos.utils import truncate_string
from atropos.utils.ngs import reverse_complement


cdef class Sequence(object):
    """
    A record in a FASTQ file. Also used for FASTA (then the qualities attribute
    is None). qualities is a string and it contains the qualities encoded as
    ascii(qual+33).

    If an adapter has been matched to the sequence, the 'match' attribute is
    set to the corresponding Match instance.

    Todo:
        `name` should be `description` (reserve .name for the part before the first
        space)
    """
    cdef:
        public str name
        public str sequence
        public str qualities
        public str name2
        public int original_length
        public object match
        public object match_info
        public object clipped
        public bint insert_overlap
        public bint merged
        public int corrected
        public str umi

    def __init__(
        self, str name, str sequence, str qualities=None, str name2='',
        original_length=None, match=None, match_info=None, clipped=None,
        insert_overlap=False, merged=False, corrected=0, str umi=None, alphabet=None,
    ):
        # Validate sequence and qualities lengths are equal
        if qualities is not None:
            slen = len(sequence)
            qlen = len(qualities)

            # check that sequence and qualities are the same length
            if slen != qlen:
                rname = truncate_string(name)
                raise FormatError(
                    f"In read named {rname!r}: length of quality sequence ({qlen}) and "
                    f"length  of read ({slen}) do not match"
                )

        # If an alphabet is specified, replace any bases not in the alphabet
        # with the default character.
        if alphabet:
            sequence = alphabet.resolve_string(sequence)

        self.name = name
        self.sequence = sequence
        self.qualities = qualities
        self.name2 = name2
        self.original_length = original_length or len(sequence)
        self.match = match
        self.match_info = match_info
        self.clipped = clipped or [0,0,0,0]
        self.insert_overlap = insert_overlap
        self.merged = merged
        self.corrected = corrected
        self.umi = umi

    def subseq(self, begin=0, end=None):
        if end is None:
            new_read = self[begin:]
        else:
            new_read = self[begin:end]
        end_bases = len(self) - end
        offset = 2 if self.match else 0
        if begin:
            new_read.clipped[offset] += begin
        if end_bases:
            new_read.clipped[offset+1] += end_bases
        return begin, end_bases, new_read

    def clip(self, front=0, back=0):
        if back < 0:
            new_read = self[front:back]
            back *= -1
        else:
            new_read = self[front:]
        offset = 2 if self.match else 0
        if front:
            new_read.clipped[offset] += front
        if back:
            new_read.clipped[offset+1] += back
        return front, back, new_read

    def reverse_complement(self):
        """
        Returns a copy of this read with the sequence reverse-complemented
        and the qualities reversed.

        TODO: add option to also reverse the alignment positions
        """
        sequence = reverse_complement(self.sequence)
        qualities = clipped = match_info = None
        if self.qualities:
            qualities = ''.join(reversed(self.qualities))
        if self.match_info:
            match_info = [copy.copy(m) for m in self.match_info]

        new_read = self.__class__(
            self.name,
            sequence,
            qualities,
            self.name2,
            self.original_length,
            None,
            match_info,
            list(self.clipped),
            self.insert_overlap,
            self.merged,
            self.corrected,
            self.umi
        )

        if self.match:
            match = self.match.copy()
            match.read = new_read
            new_read.match = match

        return new_read

    @property
    def size_in_bytes(self) -> int:
        size = len(self.name) + len(self.sequence)
        if self.name2:
            size += len(self.name2)
        if self.qualities:
            size += len(self.qualities)
        return size

    def __getitem__(self, key):
        """slicing"""
        return self.__class__(
            self.name,
            self.sequence[key],
            self.qualities[key] if self.qualities is not None else None,
            self.name2,
            self.original_length,
            self.match,
            self.match_info,
            list(self.clipped),
            self.insert_overlap,
            self.merged,
            self.corrected,
            self.umi
        )

    def __repr__(self):
        qstr = ""
        if self.qualities is not None:
            qstr = f", qualities={truncate_string(self.qualities)!r}"
        return f"<Sequence(name={truncate_string(self.name)!r}, " \
               f"sequence={truncate_string(self.sequence)!r}{qstr})>"

    def __len__(self):
        return len(self.sequence)

    def __richcmp__(self, other, int op):
        if 2 <= op <= 3:
            eq = self.name == other.name and \
                self.sequence == other.sequence and \
                self.qualities == other.qualities
            if op == 2:
                return eq
            else:
                return not eq
        else:
            raise NotImplementedError()

    def __reduce__(self):
        return Sequence, (self.name, self.sequence, self.qualities, self.name2)
