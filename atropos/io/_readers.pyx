# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
from atropos.errors import FormatError
from atropos.io.readers import PrefetchSequenceReader, estimate_num_records
from atropos.io._sequence import Sequence
from atropos.utils import classproperty


class FastqReader(PrefetchSequenceReader):
    """
    Reader for FASTQ files. Does not support multi-line FASTQ files.
    """
    @classproperty
    def file_format(cls) -> str:
        return "FASTQ"

    @classproperty
    def delivers_qualities(cls) -> bool:
        return True

    @classproperty
    def has_qualfile(cls) -> bool:
        return False

    @classproperty
    def colorspace(cls) -> bool:
        return False

    @classproperty
    def interleaved(cls) -> bool:
        return False

    def __init__(
        self,
        filename,
        quality_base=33,
        sequence_class=Sequence,
        alphabet=None
    ):
        """
        file is a filename or a file-like object.
        If file is a filename, then .gz files are supported.
        """
        super().__init__(
            filename, quality_base=quality_base, alphabet=alphabet
        )
        self.sequence_class = sequence_class

    def estimate_num_records(self):
        return estimate_num_records(self.name, self._first_seq.size_in_bytes, 4, 2)

    def _iter(self):
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
