# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
import copy
from atropos.io import xopen
from atropos.io.seqio import FormatError, SequenceReader
from atropos.util import reverse_complement, truncate_string

cdef class Sequence(object):
    """
    A record in a FASTQ file. Also used for FASTA (then the qualities attribute
    is None). qualities is a string and it contains the qualities encoded as
    ascii(qual+33).

    If an adapter has been matched to the sequence, the 'match' attribute is
    set to the corresponding Match instance.
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
    
    def __init__(self, str name, str sequence, str qualities=None, str name2='',
                 original_length=None, match=None, match_info=None, clipped=None,
                 insert_overlap=False, merged=False, corrected=0, alphabet=None):
        
        # Validate sequence and qualities lengths are equal
        if qualities is not None:
            slen = len(sequence)
            qlen = len(qualities)
            
            # check that sequence and qualities are the same length
            if slen != qlen:
                rname = truncate_string(name)
                raise FormatError(
                    "In read named {0!r}: length of quality sequence ({1}) and "
                    "length  of read ({2}) do not match".format(rname, qlen, slen))
        
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
        return (begin, end_bases, new_read)
    
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
        return (front, back, new_read)
    
    def reverse_complement(self):
        """Returns a copy of this read with the sequence reverse-complemented
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
            self.corrected
        )
        
        if self.match:
            match = self.match.copy()
            match.read = new_read
            new_read.match = match
        
        return new_read
    
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
            self.corrected
        )

    def __repr__(self):
        qstr = ''
        if self.qualities is not None:
            qstr = ', qualities={0!r}'.format(truncate_string(self.qualities))
        return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(
            truncate_string(self.name), truncate_string(self.sequence), qstr)

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
        return (Sequence, (self.name, self.sequence, self.qualities, self.name2))

class FastqReader(SequenceReader):
    """Reader for FASTQ files. Does not support multi-line FASTQ files.
    """
    file_format = "FASTQ"
    delivers_qualities = True
    
    def __init__(
            self, filename, quality_base=33, sequence_class=Sequence,
            alphabet=None):
        """
        file is a filename or a file-like object.
        If file is a filename, then .gz files are supported.
        """
        super().__init__(
            filename, quality_base=quality_base, alphabet=alphabet)
        self.sequence_class = sequence_class
    
    def __iter__(self):
        """
        Yield Sequence objects
        """
        cdef int i = 0
        cdef int strip
        cdef str line, name, qualities, sequence, name2
        sequence_class = self.sequence_class

        it = iter(self._file)
        line = next(it)
        if not (line and line[0] == '@'):
            raise FormatError("Line {0} in FASTQ file is expected to start "
                              "with '@', but found {1!r}".format(i+1, line[:10]))
        strip = -2 if line.endswith('\r\n') else -1
        name = line[1:strip]

        i = 1
        for line in it:
            if i == 0:
                if not (line and line[0] == '@'):
                    raise FormatError("Line {0} in FASTQ file is expected to "
                                      "start with '@', but found {1!r}".format(
                                      i+1, line[:10]))
                name = line[1:strip]
            elif i == 1:
                sequence = line[:strip]
            elif i == 2:
                if line == '+\n':  # check most common case first
                    name2 = ''
                else:
                    line = line[:strip]
                    if not (line and line[0] == '+'):
                        raise FormatError(
                            "Line {0} in FASTQ file is expected to start with "
                            "'+', but found {1!r}".format(i+1, line[:10]))
                    if len(line) > 1:
                        if not line[1:] == name:
                            raise FormatError(
                                "At line {0}: Sequence descriptions in the FASTQ file don't match "
                                "({1!r} != {2!r}).\n"
                                "The second sequence description must be either empty "
                                "or equal to the first description.".format(i+1,
                                    name, line[1:]))
                        name2 = name
                    else:
                        name2 = ''
            elif i == 3:
                if len(line) == len(sequence) - strip:
                    qualities = line[:strip]
                else:
                    qualities = line.rstrip('\r\n')
                try:
                    yield sequence_class(
                        name, sequence, qualities, name2=name2,
                        alphabet=self.alphabet)
                except Exception as err:
                    raise FormatError(
                        "Error creating sequence record at line "
                        "{}".format(i+1)) from err
            i = (i + 1) % 4
        if i != 0:
            raise FormatError("FASTQ file ended prematurely")
