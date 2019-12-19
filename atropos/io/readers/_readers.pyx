# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
from atropos.errors import FormatError
from atropos.utils import classproperty

from ._base import PrefetchSequenceReader, estimate_num_records


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

    def estimate_num_records(self):
        return estimate_num_records(self.name, self._first_seq.size_in_bytes, 4, 2)

    def _iter(self):
        """
        Yield Sequence objects
        """
        cdef int i = 0
        cdef int strip
        cdef str line, name, qualities, sequence, name2
        sequence_factory = self._sequence_factory

        it = iter(self._file)
        line = next(it)
        if not (line and line[0] == '@'):
            raise FormatError(
                f"Line {i+1} in FASTQ file is expected to start with '@', but found "
                f"{line[:10]!r}"
            )
        strip = -2 if line.endswith('\r\n') else -1
        name = line[1:strip]

        i = 1
        for line in it:
            if i == 0:
                if not (line and line[0] == '@'):
                    raise FormatError(
                        "Line {i+1} in FASTQ file is expected to start with '@', but "
                        "found {line[:10]!r}"
                    )
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
                            f"Line {i+1} in FASTQ file is expected to start with '+', "
                            f"but found {line[:10]!r}"
                        )
                    if len(line) > 1:
                        if not line[1:] == name:
                            raise FormatError(
                                f"At line {i+1}: Sequence descriptions in the FASTQ "
                                f"file don't match ({name!r} != {line[1:]!r}).\n"
                                f"The second sequence description must be either empty "
                                f"or equal to the first description."
                            )
                        name2 = name
                    else:
                        name2 = ''
            elif i == 3:
                if len(line) == len(sequence) - strip:
                    qualities = line[:strip]
                else:
                    qualities = line.rstrip('\r\n')
                try:
                    yield sequence_factory(
                        name, sequence, qualities, name2=name2, alphabet=self._alphabet
                    )
                except Exception as err:
                    raise FormatError(
                        f"Error creating sequence record at line {i+1}"
                    ) from err
            i = (i + 1) % 4
        if i != 0:
            raise FormatError("FASTQ file ended prematurely")
