from abc import ABCMeta, abstractmethod
import os
from pathlib import Path
from typing import (
    Callable, Iterable, Iterator, Optional, Sequence as SequenceType, Tuple, Union
)

from xphyle import xopen
from xphyle.types import ModeArg, PathOrFile
from xphyle.utils import uncompressed_size

from atropos.io import InputRead
from atropos.io.sequence import Sequence
from atropos.utils import classproperty
from atropos.utils.ngs import Alphabet
from atropos.utils.collections import Summarizable


LINESEP_LEN = len(os.linesep)


class SequenceReaderBase(
    Summarizable, Iterable[Tuple[Sequence, ...]], metaclass=ABCMeta
):
    """
    Base class for sequence readers.
    """

    @classproperty
    @abstractmethod
    def file_format(cls) -> str:
        pass

    @classproperty
    @abstractmethod
    def delivers_qualities(cls) -> bool:
        pass

    @classproperty
    @abstractmethod
    def has_qualfile(cls) -> bool:
        pass

    @classproperty
    @abstractmethod
    def colorspace(cls) -> bool:
        pass

    @classproperty
    @abstractmethod
    def interleaved(cls) -> bool:
        pass

    @classproperty
    def input_read(self) -> Optional[InputRead]:
        return None

    @property
    @abstractmethod
    def input_names(self) -> Tuple[Union[str, SequenceType[str]], Optional[str]]:
        pass

    @property
    @abstractmethod
    def quality_base(self) -> int:
        pass

    @property
    def paired(self) -> bool:
        return self.input_read == InputRead.PAIRED

    def estimate_num_records(self) -> Optional[int]:
        """
        Estimates the total number of records in the file(s), if possible.

        The way this is typically done is to divide the size of the file
        (obtained via xphyle.utils.uncompressed_size) by the size of the first
        record.

        Returns:
            An estimate of the number of records (or record pairs), or None if
            the number cannot be estimated.
        """

    @abstractmethod
    def __iter__(self) -> Iterator[Tuple[Sequence, ...]]:
        pass

    def close(self):
        pass

    def summarize(self):
        return dict(
            input_names=self.input_names,
            input_read=self.input_read,
            file_format=self.file_format,
            delivers_qualities=self.delivers_qualities,
            quality_base=self.quality_base,
            has_qualfile=self.has_qualfile,
            colorspace=self.colorspace,
            interleaved=self.interleaved,
        )


class SequenceReader(SequenceReaderBase, metaclass=ABCMeta):
    """
    Reads possibly compressed files containing sequences.
    """

    def __init__(
        self,
        path,
        mode: Optional[ModeArg] = "r",
        quality_base: Optional[int] = 33,
        alphabet: Optional[Alphabet] = None,
        sequence_factory: Callable[..., Sequence] = Sequence,
        close_on_exit: Optional[bool] = None
    ):
        """
        Args:
            path: Path or file-like object that may be compressed (.gz, .bz2, .xz).
            mode: The file open mode.
            quality_base: The minimum quality value.
            alphabet: The alphabet to use to validate sequences. If None,
                no validation is done.
            close_on_exit: Whether to close the source file on exit.
        """
        if path is None:
            raise ValueError("'path' cannot be None")

        # TODO: if quality_base is None, detect it from the data
        self._quality_base = quality_base
        self._alphabet = alphabet
        self._sequence_factory = sequence_factory

        is_path = isinstance(path, (str, Path))

        if close_on_exit is None:
            self._close_on_exit = is_path
        else:
            self._close_on_exit = close_on_exit

        if is_path:
            self.name = str(path)
            self._file = xopen(path, mode)
            if close_on_exit is None:
                self._close_on_exit = True
        else:
            if hasattr(path, "name"):
                self.name = path.name
            else:
                # TODO: generate random unique name?
                self.name = path.__class__

            self._file = path

    @property
    def input_names(self):
        return self.name, None

    @property
    def quality_base(self) -> Optional[int]:
        return self._quality_base

    def close(self):
        """
        Closes the underlying file.
        """
        if (
            self._close_on_exit and
            self._file is not None
        ):
            if hasattr(self._file, "close"):
                self._file.close()

            self._file = None

    def __enter__(self):
        if self._file is None:
            raise ValueError("I/O operation on closed SequenceReader")

        return self

    def __exit__(self, *args):
        self.close()


class PrefetchSequenceReader(SequenceReader, metaclass=ABCMeta):
    """
    SequenceReader that prefetches the first line and stores it so it can be
    inspected before iteration begins.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._seq_iter = self._iter()
        try:
            self._first_seq = next(self._seq_iter)
        except StopIteration:
            self._first_seq = None

    def __iter__(self) -> Iterator[Tuple[Sequence, ...]]:
        if self._first_seq is None:
            return
        yield self._first_seq
        yield from self._seq_iter

    @abstractmethod
    def _iter(self) -> Iterator[Tuple[Sequence, ...]]:
        pass


def estimate_num_records(
    input_file, record_size, lineseps, format_chars=0, header_size=0
) -> Optional[int]:
    """
    Convenience method to compute the estimated number of records in the input file(s).

    Args:
        input_file: Path to input file.
        record_size: The approximate size of a single record (without
            whitespace or formatting characters).
        lineseps: Number of line separators in a record (including one at
            the end).
        format_chars: Number of formatting chars, e.g. '>' in FASTA and '@'
            and '+' in FASTQ.

    Returns:
        Integer estimate of the number of records in the file.
    """
    file_size = uncompressed_size(input_file) - header_size

    if file_size is not None:
        return int(file_size / (record_size + (lineseps * LINESEP_LEN) + format_chars))


def sequence_names_match(read1: Sequence, read2: Sequence) -> bool:
    """
    Checks whether the sequences read1 and read2 have identical names, ignoring a
    suffix of '1' or '2'. Some old paired-end reads have names that end in '/1' and
    '/2'. Also, the fastq-dump tool (used for converting SRA files to FASTQ) appends
    a .1 and .2 to paired-end reads if option -I is used.

    Args:
        read1, read2: The sequences to compare.

    Returns:
        Whether the sequences are equal.
    """
    name1 = read1.name.split(None, 1)[0]
    name2 = read2.name.split(None, 1)[0]

    if name1[-1:] in "12" and name2[-1:] in "12":
        name1 = name1[:-1]
        name2 = name2[:-1]

    return name1 == name2
