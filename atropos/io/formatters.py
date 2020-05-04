from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import List, Optional, Tuple, Union, cast

from atropos.errors import UnknownFileTypeError
from atropos.io import SequenceFileType
from atropos.io.sequence import Sequence, ColorspaceSequence


DEFAULT_SAM_HEADER = "@HD\tVN:1.6\tSO:unsorted\n"


class SequenceFileFormat(metaclass=ABCMeta):
    """
    Base class for sequence formatters.
    """

    @abstractmethod
    def format(self, read: Sequence) -> str:
        """Formats a Sequence as a string.

        Args:
            read: The Sequence object.

        Returns:
            A string representation of the sequence object in the sequence
            file format.
        """


class FastaFormat(SequenceFileFormat):
    """
    FASTA SequenceFileFormat.
    """

    def __init__(self, line_length=None):
        """
        Args:
            line_length: Max line length (in characters), or None. Determines
                whether and how lines are wrapped.
        """
        if line_length:
            from textwrap import TextWrapper

            self._text_wrapper = TextWrapper(width=line_length)
        else:
            self._text_wrapper = None

    def format(self, read: Sequence) -> str:
        return self.format_entry(read.name, read.sequence)

    def format_entry(self, name: str, sequence: str) -> str:
        """
        Converts a sequence record to a string.
        """
        if self._text_wrapper:
            sequence = self._text_wrapper.fill(sequence)

        return f">{name}\n{sequence}\n"


class ColorspaceFastaFormat(FastaFormat):
    """
    FastaFormat in which sequences are in colorspace.
    """

    def format(self, read: ColorspaceSequence) -> str:
        return self.format_entry(read.name, read.primer + read.sequence)


class FastqFormat(SequenceFileFormat):
    """
    FASTQ SequenceFileFormat.
    """

    def format(self, read: Sequence) -> str:
        return self.format_entry(read.name, read.sequence, read.qualities, read.name2)

    @classmethod
    def format_entry(
        cls, name: str, sequence: str, qualities: str, name2: Optional[str] = None
    ) -> str:
        """
        Converts a sequence record to a string.
        """
        return f"@{name}\n{sequence}\n+{name2 or ''}\n{qualities}\n"


class ColorspaceFastqFormat(FastqFormat):
    """
    FastqFormat in which sequences are in colorspace.
    """

    def format(self, read: ColorspaceSequence) -> str:
        return self.format_entry(
            read.name, read.primer + read.sequence, read.qualities, read.name2
        )


class SAMFormat(SequenceFileFormat):
    """
    SAM SequenceFileFormat.
    """

    def __init__(self, flag: int):
        self._flag = str(flag)

    def format(self, read: Sequence) -> str:
        umi = f"\tBC:Z:{read.umi}" if read.umi else ""
        return f"{read.name}\t{self._flag}\t*\t0\t0\t*\t*\t0\t0\t{read.sequence}\t" \
               f"{read.qualities}{umi}\n"


class Formatter(metaclass=ABCMeta):
    """
    Base class for Formatters.
    """

    def __init__(self):
        self.written = 0
        self.read1_bp = 0
        self.read2_bp = 0

    @abstractmethod
    def format(
        self, result: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        """
        Formats read(s) and add them to `result`.

        Args:
            result: A dict mapping file names to lists of formatted reads.
            read1: The first read to format.
            read2: The second read to format.
        """

    def _get_result_list(self, result: dict, file1: str) -> List[str]:
        if file1 in result:
            return result[file1]
        else:
            result_list = []

            header = self._get_header()

            if header:
                result_list.append(header)

            result[file1] = result_list

            return result_list

    def _get_header(self) -> Optional[str]:
        pass

    @property
    def written_bp(self) -> Tuple[int, int]:
        """
        Tuple of base-pairs written (read1_bp, read2_bp).
        """
        return self.read1_bp, self.read2_bp


class SingleEndFormatter(Formatter):
    """
    Wrapper for a SequenceFileFormat for single-end data.
    """

    def __init__(self, file1: str, seq_format: SequenceFileFormat):
        """
        Args:
            file1: The single-end file.
            seq_format: The SequenceFileFormat object.
        """
        super().__init__()
        self._file1 = file1
        self._seq_format = seq_format

    def format(
        self, result: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        """
        Formats read(s) and add them to `result`.

        Args:
            result: A dict mapping file names to lists of formatted reads.
            read1:
            read2: The reads to format.
        """
        result_list = self._get_result_list(result, self._file1)
        result_list.append(self._seq_format.format(read1))
        self.written += 1
        self.read1_bp += len(read1)


class InterleavedFormatter(Formatter):
    """
    Formats read pairs as successive reads in an interleaved file.
    """

    def __init__(
        self,
        file1: str,
        seq_format1: SequenceFileFormat,
        seq_format2: Optional[SequenceFileFormat] = None,
    ):
        super().__init__()
        self._file1 = file1
        self._seq_format1 = seq_format1
        self._seq_format2 = seq_format2 or seq_format1

    def format(
        self, result: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        result_list = self._get_result_list(result, self._file1)
        result_list.append(self._seq_format1.format(read1))
        result_list.append(self._seq_format2.format(read2))
        self.written += 1
        self.read1_bp += len(read1)
        self.read2_bp += len(read2)


class PairedEndFormatter(Formatter):
    """
    Wrapper for a SequenceFileFormat. Both reads in a pair are formatted using the
    specified format.
    """

    def __init__(
        self,
        file1: str,
        file2: str,
        seq_format1: SequenceFileFormat,
        seq_format2: Optional[SequenceFileFormat] = None,
    ):
        super().__init__()
        self._data = [(file1, seq_format1), (file2, seq_format2 or seq_format1)]

    def format(
        self, result: dict, read1: Sequence, read2: Optional[Sequence] = None
    ):
        for read, (path, seq_format) in zip((read1, read2), self._data):
            result_list = self._get_result_list(result, path)
            result_list.append(seq_format.format(read))

        self.written += 1
        self.read1_bp += len(read1)
        self.read2_bp += len(read2)


class SAMFormatterMixin:
    def __init__(self, *args, header: Optional[dict] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self._header = header

    def _get_header(self) -> Optional[str]:
        if not self._header:
            return DEFAULT_SAM_HEADER

        def create_row(_tags: Union[str, dict]):
            if isinstance(_tags, dict):
                tags_str = "\t".join(f"{key}:{val}" for key, val in _tags.items())
            else:
                tags_str = cast(str, _tags)
            return f"@{header_type}\t{tags_str}"

        rows = []

        for header_type, header_value in self._header.items():
            if isinstance(header_value, dict):
                rows.append(create_row(header_value))
            else:
                rows.extend(create_row(tags) for tags in header_value)

        return "\n".join(rows) + "\n"


class SingleEndSAMFormatter(SAMFormatterMixin, SingleEndFormatter):
    def __init__(self, file1: str, header: Optional[dict] = None):
        super().__init__(file1, SAMFormat(0), header=header)


class PairedEndSAMFormatter(SAMFormatterMixin, InterleavedFormatter):
    def __init__(self, file1: str, header: Optional[dict] = None):
        super().__init__(file1, SAMFormat(65), SAMFormat(129), header=header)


def sra_colorspace_sequence(
    name: str, sequence: str, qualities: str, name2: str
):
    """
    Factory for an SRA colorspace sequence (which has one quality value too many).
    """
    return ColorspaceSequence(name, sequence, qualities[1:], name2=name2)


def create_seq_formatter(
    file1: Union[str, Path],
    file2: Optional[Union[str, Path]] = None,
    qualities: Optional[bool] = None,
    colorspace: bool = False,
    file_format: Optional[Union[SequenceFileType, SequenceFileFormat]] = None,
    interleaved: bool = False,
    line_length: Optional[int] = None,
    bam_header: Optional[dict] = None,
):
    """
    Creates a Formatter, deriving the format name from the file extension.

    Args:
        file1: First output file.
        file2: Second output file.
        qualities: When file_format is None, this can be set to True or False to
            specify whether the written sequences will have quality values. This is
            is used in two ways:
            * If the output format cannot be determined (unrecognized extension
              etc), no exception is raised, but fasta or fastq format is chosen
              appropriately.
            * When False (no qualities available), an exception is raised when
              the auto-detected output format is FASTQ.
        colorspace: If True, instances of the Colorspace... formats are returned.
        file_format: If set to None, file format is autodetected from the file
            name extension. Set to FASTA, FASTQ, or SAM to not auto-detect.
            Colorspace is not auto-detected and must always be requested explicitly.
        interleaved: Whether the output should be interleaved (file2 must be None).
        line_length: Maximum length of a sequence line in FASTA output.
        bam_header: If writing SAM format, optional headers to add at the beginning of
            the output.

    Returns:
        A Formatter instance.
    """
    if file_format is None:
        file_format = SequenceFileType.guess_from_name(
            file1, output=True, raise_on_failure=qualities is None
        )

    if file_format is None:
        if qualities is True:
            # Format not recognized, but know we want to write reads with qualities.
            file_format = SequenceFileType.FASTQ
        elif qualities is False:
            # Same, but we know that we want to write reads without qualities.
            file_format = SequenceFileType.FASTA
        else:
            raise UnknownFileTypeError("Could not determine file type.")

    if (
        file_format in (SequenceFileType.FASTQ, SequenceFileType.SAM)
        and qualities is False
    ):
        raise ValueError(
            f"Output format cannot be {file_format} since no quality values are "
            f"available."
        )

    if file_format == SequenceFileType.SAM:
        if file2 is not None:
            raise ValueError("Only one output file allowed for SAM format")

        if interleaved:
            return PairedEndSAMFormatter(file1, bam_header)
        else:
            return SingleEndSAMFormatter(file1, bam_header)
    else:
        if file_format == SequenceFileType.FASTA:
            if colorspace:
                fmt = ColorspaceFastaFormat(line_length)
            else:
                fmt = FastaFormat(line_length)
        elif file_format == SequenceFileType.FASTQ:
            if colorspace:
                fmt = ColorspaceFastqFormat()
            else:
                fmt = FastqFormat()
        else:
            raise UnknownFileTypeError(
                f"File format {file_format!r} is unknown (expected 'fasta' or "
                f"'fastq')."
            )

        if file2 is not None:
            return PairedEndFormatter(file1, file2, fmt)
        elif interleaved:
            return InterleavedFormatter(file1, fmt)
        else:
            return SingleEndFormatter(file1, fmt)
