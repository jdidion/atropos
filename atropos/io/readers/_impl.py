from abc import ABCMeta, abstractmethod
from collections import namedtuple
import csv
import math
from pathlib import Path
import sys
from typing import (
    Callable, Iterator, IO, Optional, Sequence as SequenceType, Tuple, Type,
    Union, cast
)

from xphyle import STDIN, xopen

from atropos.errors import FormatError, UnknownFileTypeError
from atropos.io import InputRead, guess_format_from_name
from atropos.io.sequence import Sequence, ColorspaceSequence
from atropos.utils import classproperty, truncate_string
from atropos.utils.ngs import ALPHABETS, Alphabet

from ._base import (
    SequenceReader,
    SequenceReaderBase,
    PrefetchSequenceReader,
    estimate_num_records,
    sequence_names_match
)
from ._readers import FastqReader


class FastaReader(PrefetchSequenceReader):
    """
    Reader for FASTA files.

    Args:
        path: A path or a file-like object. In both cases, the file may
            be compressed (.gz, .bz2, .xz).
        keep_linebreaks: Whether to keep newline characters in the sequence.
        sequence_factory: The class to use when creating new sequence objects.
    """

    @classproperty
    def file_format(cls) -> str:
        return "FASTA"

    @classproperty
    def delivers_qualities(cls) -> bool:
        return False

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
        path: Union[str, Path, Iterator[str]],
        keep_linebreaks: bool = False,
        **kwargs
    ):
        self._delimiter = "\n" if keep_linebreaks else ""
        super().__init__(path, **kwargs)

    def _iter(self):
        """
        Reads the next entry from the file (single entry at a time).
        """
        name = None
        seq = []

        for i, line in enumerate(self._file):
            # strip() also removes DOS line breaks
            line = line.strip()

            if not line:
                continue

            if line and line[0] == ">":
                if name is not None:
                    yield self._sequence_factory(
                        name=name,
                        sequence=self._delimiter.join(seq),
                        alphabet=self._alphabet
                    )

                name = line[1:]
                seq = []
            elif line and line[0] == "#":
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FormatError(
                    f"At line {i + 1}: Expected '>' at beginning of FASTA record, "
                    f"but got {truncate_string(line)!r}."
                )

        if name is not None:
            yield self._sequence_factory(name, self._delimiter.join(seq), None)

    def estimate_num_records(self):
        # TODO: this will underestimate the record size (and thus overestimate the
        #  total number of records) for FASTA files with wrapped sequence lines.
        record_size = sum(seq.size_in_bytes() for seq in self._first_seq)
        return estimate_num_records([self.name], record_size, 2, 1)


class NgstreamSequenceReader(SequenceReader):
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

    def __init__(self, reader, **kwargs):
        import ngstream

        protocol = cast(ngstream.Protocol, reader)

        self._input_read = InputRead.PAIRED if protocol.paired else InputRead.SINGLE
        super().__init__(protocol, **kwargs)

    @property
    def input_read(self) -> InputRead:
        return self._input_read

    def estimate_num_records(self):
        return self._file.read_count

    def __iter__(self):
        if self._input_read == InputRead.PAIRED:
            for read in self._file:
                yield tuple(self._as_sequence(read[i]) for i in (0, 1))
        else:
            for read in self._file:
                yield self._as_sequence(read),

    def _as_sequence(self, frag):
        return self._sequence_factory(
            frag.name,
            frag.sequence,
            frag.qualities,
            alphabet=self._alphabet
        )

    def close(self):
        self._file.finish()


class FastaQualReader(SequenceReaderBase):
    """
    Reader for reads that are stored in .(CS)FASTA and .QUAL files.
    """

    @classproperty
    def file_format(cls) -> str:
        return "FastaQual"

    @classproperty
    def delivers_qualities(cls) -> bool:
        return True

    @classproperty
    def has_qualfile(cls) -> bool:
        return True

    @classproperty
    def colorspace(cls) -> bool:
        return False

    @classproperty
    def interleaved(cls) -> bool:
        return False

    @classproperty
    def input_read(self) -> InputRead:
        return InputRead.SINGLE

    def __init__(
        self,
        fastafile: Union[str, Path, Iterator[str]],
        qualfile: Union[str, Path, Iterator[str]],
        quality_base: int = 33,
        alphabet: Optional[Alphabet] = None,
        sequence_factory: Callable[..., Sequence] = Sequence,
    ):
        """
        Args:
            fastafile: filename or file-like object.
            qualfile: filename or file-like object.
            quality_base:
            sequence_factory: The class to use when creating new sequence objects.
            alphabet:
        """
        self._fastareader = FastaReader(fastafile)
        self._qualreader = FastaReader(qualfile, keep_linebreaks=True)
        self._quality_base = quality_base
        self._alphabet = alphabet
        self._sequence_factory = sequence_factory

    @property
    def input_names(self) -> Tuple[Tuple[str, str], None]:
        return (self._fastareader.name, self._qualreader.name), None

    @property
    def quality_base(self) -> Optional[int]:
        return self._quality_base

    def estimate_num_records(self):
        return self._fastareader.estimate_num_records()

    def __iter__(self):
        """
        Yields Sequence objects.
        """
        # conversion dictionary: maps strings to the appropriate ASCII-encoded
        # character
        conv = dict((str(i), chr(i + 33)) for i in range(-5, 256 - 33))

        for fastaread, qualread in zip(self._fastareader, self._qualreader):
            if fastaread.name != qualread.name:
                raise FormatError(
                    f"The read names in the FASTA and QUAL file do not match "
                    f"({fastaread.name!r} != {qualread.name!r})"
                )

            try:
                qualities = "".join(
                    [conv[value] for value in qualread.sequence.split()]
                )
            except KeyError as err:
                raise FormatError(
                    f"Within read named {fastaread.name!r}: Found invalid quality "
                    f"value {err}"
                )

            yield self._sequence_factory(
                name=fastaread.name,
                sequence=fastaread.sequence,
                qualities=qualities,
                alphabet=self._alphabet
            )

    def close(self):
        """
        Closes the underlying files.
        """
        self._fastareader.close()
        self._qualreader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class ColorspaceSequenceReaderMixin:
    """
    Base class for colorspace sequence readers.
    """

    @classproperty
    def colorspace(cls) -> bool:
        return True

    def __init__(
        self,
        reader,
        quality_base: int = 33,
        sequence_factory: Callable[..., Sequence] = ColorspaceSequence,
        **kwargs
    ):
        super().__init__(
            reader,
            quality_base=quality_base,
            sequence_factory=sequence_factory,
            **kwargs
        )


class ColorspaceNgstreamSequenceReader(
    ColorspaceSequenceReaderMixin, NgstreamSequenceReader
):
    """
    Reads colorspace sequences from an SRA accession.
    """


class ColorspaceFastaReader(ColorspaceSequenceReaderMixin, FastaReader):
    """
    Reads colorspace sequences from a FASTA.
    """


class ColorspaceFastqReader(ColorspaceSequenceReaderMixin, FastqReader):
    """
    Reads colorspace sequences from a FASTQ.
    """


class ColorspaceFastaQualReader(ColorspaceSequenceReaderMixin, FastaQualReader):
    """
    Reads sequences and qualities from separate files and returns
    :class:`ColorspaceSequence`s.
    """


class PairedSequenceReader(SequenceReaderBase):
    """
    Reads paired-end reads from two files. Wraps two SequenceReader instances,
    making sure that reads are properly paired.
    """

    @classproperty
    def has_qualfile(cls) -> bool:
        return False

    @classproperty
    def interleaved(cls) -> bool:
        return False

    @classproperty
    def input_read(self) -> InputRead:
        return InputRead.PAIRED

    def __init__(
        self,
        file1,
        file2,
        quality_base: int = 33,
        colorspace: bool = False,
        file_format: bool = None,
        alphabet: Optional[Alphabet] = None
    ):
        """
        Args:
            file1, file2: The pair of files.
            colorspace: Whether the sequences are in colorspace.
            file_format: A file_format instance.
        """
        self._reader1 = open_reader(
            file1,
            colorspace=colorspace,
            quality_base=quality_base,
            file_format=file_format,
            alphabet=alphabet,
        )
        self._reader2 = open_reader(
            file2,
            colorspace=colorspace,
            quality_base=quality_base,
            file_format=file_format,
            alphabet=alphabet,
        )

    @property
    def input_names(self):
        return self._reader1.input_names[0], self._reader2.input_names[0]

    @property
    def file_format(self) -> str:
        return self._reader1.file_format

    def colorspace(self) -> bool:
        return self._reader1.colorspace

    @property
    def quality_base(self) -> int:
        return self._reader1.quality_base

    def delivers_qualities(self) -> bool:
        return self._reader1.delivers_qualities

    def estimate_num_records(self) -> Optional[int]:
        ests = list(
            filter(
                None,
                (
                    self._reader1.estimate_num_records(),
                    self._reader2.estimate_num_records()
                ),
            )
        )

        if len(ests) > 0:
            return max(ests)

    def __getattr__(self, name):
        return getattr(self._reader1, name)

    def __iter__(self):
        """
        Iterates over the paired reads. Each item is a pair of Sequence objects.
        """
        # Avoid usage of zip() below since it will consume one item too many.
        it1 = iter(self._reader1)
        it2 = iter(self._reader2)

        while True:
            try:
                read1 = next(it1)
            except StopIteration:
                # End of file 1. Make sure that file 2 is also at end.
                try:
                    next(it2)

                    raise FormatError(
                        "Reads are improperly paired. There are more reads in file 2 "
                        "than in file 1."
                    )
                except StopIteration:
                    pass

                break

            try:
                read2 = next(it2)
            except StopIteration:
                raise FormatError(
                    "Reads are improperly paired. There are more reads in file 1 than "
                    "in file 2."
                )

            if not sequence_names_match(read1, read2):
                raise FormatError(
                    f"Reads are improperly paired. Read name '{read1.name}' in file 1 "
                    f"does not match '{read2.name}' in file 2."
                )

            yield read1, read2

    def close(self):
        """
        Closes the underlying files.
        """
        self._reader1.close()
        self._reader2.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class SequenceReaderWrapper(SequenceReaderBase):
    def __init__(self, reader: SequenceReaderBase):
        """
        Args:
            reader: The reader to wrap.
        """
        self._reader = reader

    @property
    def file_format(self) -> str:
        return self._reader.file_format

    @property
    def delivers_qualities(self) -> bool:
        return self._reader.delivers_qualities

    @property
    def colorspace(self) -> bool:
        return self._reader.colorspace

    @property
    def input_names(self) -> Tuple[Union[str, SequenceType[str]], Optional[str]]:
        return self._reader.input_names

    @property
    def quality_base(self) -> int:
        return self._reader.quality_base

    @property
    def has_qualfile(self) -> bool:
        return self._reader.has_qualfile

    @property
    def interleaved(self) -> bool:
        return self._reader.interleaved

    @property
    def input_read(self) -> InputRead:
        return self._reader.input_read

    def estimate_num_records(self) -> Optional[int]:
        return self._reader.estimate_num_records()

    def __iter__(self):
        return iter(self._reader)

    def close(self):
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class InterleavedSequenceReader(SequenceReaderWrapper):
    """
    Reads paired-end reads from an interleaved FASTQ file.
    """

    @classproperty
    def has_qualfile(cls) -> bool:
        return False

    @classproperty
    def interleaved(cls) -> bool:
        return True

    @classproperty
    def input_read(self) -> InputRead:
        return InputRead.PAIRED

    def __init__(
        self,
        path,
        quality_base: int = 33,
        colorspace: bool = False,
        file_format: Optional[str] = None,
        alphabet: Optional[Alphabet] = None
    ):
        """
        Args:
            path: The interleaved FASTQ file.
            colorspace: Whether the sequences are in colorspace.
            file_format: A file_format instance.
        """
        super().__init__(open_reader(
            path,
            quality_base=quality_base,
            colorspace=colorspace,
            file_format=file_format,
            alphabet=alphabet,
        ))

    def estimate_num_records(self) -> int:
        return math.ceil(self._reader.estimate_num_records() / 2)

    def __iter__(self) -> Iterator[Tuple[Sequence, Sequence]]:
        # Avoid usage of zip() below since it will consume one item too many.
        itr = iter(self._reader)

        for read1 in itr:
            try:
                read2 = next(itr)
            except StopIteration:
                raise FormatError(
                    "Interleaved input file incomplete: Last record has no " "partner."
                )

            if not sequence_names_match(read1, read2):
                raise FormatError(
                    f"Reads are improperly paired. Name {read1.name!r} (first) does "
                    f"not match {read2.name!r} (second)."
                )

            yield read1, read2


class PairedToSingleEndReader(SequenceReaderWrapper):
    """
    Wrapper that yields either first or second reads from a paired-end reader.
    """

    def __init__(self, reader: SequenceReaderBase, input_read: InputRead):
        super().__init__(reader)
        self._input_read = input_read

    def __iter__(self):
        idx = 1 if self._input_read == InputRead.READ2 else 0
        for record in self._reader:
            yield record[idx]


class SAMReader(SequenceReaderBase, metaclass=ABCMeta):
    """
    Reader for SAM/BAM files. Paired-end files must be name-sorted. Does not support
    secondary/supplementary reads.

    TODO: either drop support for BAM files (require user to pipe output from
     `samtools view`) and implement pure-python SAM parsing, or replace pysam with a
     different BAM-parsing library (pysam is not Python 3.8 compatible).
    """

    @classproperty
    def file_format(cls) -> str:
        return "SAM"

    @classproperty
    def delivers_qualities(cls) -> bool:
        return True

    @classproperty
    def has_qualfile(cls) -> bool:
        return False

    @classproperty
    def colorspace(cls) -> bool:
        return False

    def __init__(
        self,
        path: Union[str, Path, Iterator[str], Iterator[bytes]],
        quality_base: int = 33,
        sequence_factory: Type[Sequence] = Sequence,
        alphabet: Alphabet = None,
        pysam_kwargs: Optional[dict] = None,
    ):
        """
        Args:
            path: A filename or a file-like object. If a filename, then .gz files
                are supported.
            sequence_factory: The class to use when creating new sequence objects.
        """
        self._close_on_exit = is_path = isinstance(path, (str, Path))
        self._quality_base = quality_base
        self._sequence_factory = sequence_factory
        self._alphabet = alphabet

        if is_path:
            self.name = str(path)
        elif hasattr(path, "name"):
            self.name = getattr(path, "name")
        else:
            self.name = path.__class__

        if self.name.endswith(".sam"):
            if is_path:
                self._file = xopen(path, "rt")
            else:
                self._file = cast(IO, path)

            self._sam_iter = SAMParser(self._file)
        elif self.name.endswith(".bam"):
            if is_path:
                self._file = xopen(path, "rb")
            else:
                self._file = cast(IO, path)

            if pysam_kwargs is None:
                pysam_kwargs = {"check_sq": False}

            self._sam_iter = BAMParser(self._file, **pysam_kwargs)
        else:
            # Todo: open file and check for BAM magic number
            raise ValueError(f"Cannot detect type of file {path}")

    @property
    def input_names(self):
        return self.name, None

    @property
    def quality_base(self) -> int:
        return self._quality_base

    def estimate_num_records(self):
        return self._sam_iter.estimate_num_records()

    def __iter__(self):
        # pysam raises an error of the SAM/BAM header is not complete,
        # even though it's unnecessary for our puroses. In that case,
        # we use our own simple SAM parser.
        return self._iter(self._sam_iter)

    @abstractmethod
    def _iter(self, sam):
        """
        Creates an iterator over records in the SAM/BAM file.
        """

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """
        Closes the underling AlignmentFile.
        """
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def _as_sequence(self, read):
        return self._sequence_factory(
            read.query_name,
            read.query_sequence,
            read.query_qualities,
            alphabet=self._alphabet
        )


SAMRead = namedtuple(
    "SAMRead",
    ("query_name", "query_sequence", "query_qualities", "is_read1", "is_read2"),
)


class SAMParser:
    def __init__(self, sam_file: IO):
        self._sam_file = sam_file
        self._reader = csv.reader(sam_file, delimiter="\t")
        self._header_size = 0

        line = None

        for line in self._reader:
            # skip header lines
            if not line[0].startswith("@"):
                self._header_size += len(line)
                break

        self._next_line = line

    def estimate_num_records(self):
        record_len = len("\t".join(self._next_line))
        return estimate_num_records(
            self._sam_file, record_len, 1, header_size=self._header_size
        )

    def __iter__(self):
        return self

    def __next__(self):
        if self._next_line is None:
            raise StopIteration()

        is_read1 = (int(self._next_line[1]) & 64) > 0

        read = SAMRead(
            self._next_line[0],
            self._next_line[9],
            self._next_line[10],
            is_read1,
            not is_read1,
        )

        try:
            self._next_line = next(self._reader)
        except StopIteration:
            self._next_line = None

        return read


class BAMParser:
    def __init__(self, bam_file: IO, **kwargs):
        import pysam

        self._reader = pysam.AlignmentFile(str(bam_file), **kwargs)

    def estimate_num_records(self):
        """
        Todo: Not sure how to estimate the number of records in a name-sorted BAM file.
         There is no index, and bgzip doesn't have an option to get the uncompressed
         size. Not sure if gzip -l will work correctly for a bgzip file.
        """

    def __iter__(self):
        return self

    def __next__(self):
        read = next(self._reader)

        return SAMRead(
            read.query_name,
            read.query_sequence,
            "".join(chr(33 + q) for q in read.query_qualities),
            read.is_read1,
            read.is_read2,
        )


class SingleEndSAMReader(SAMReader):
    """
    Reader for single-end SAM/BAM files.
    """
    @classproperty
    def interleaved(cls) -> bool:
        return False

    @property
    def input_read(self) -> InputRead:
        return InputRead.SINGLE

    def _iter(self, sam):
        for read in sam:
            yield self._as_sequence(read)


class PairedEndEstimatorMixin:
    def estimate_num_records(self):
        est = cast(SequenceReader, super()).estimate_num_records()
        if est is None:
            return None
        else:
            return math.ceil(est / 2)


class Read1SingleEndSAMReader(SAMReader, PairedEndEstimatorMixin):
    """
    Reads a paired-end SAM/BAM file as if it were single-end, yielding only the
    first read from each pair.
    """
    @classproperty
    def interleaved(cls) -> bool:
        return False

    @property
    def input_read(self) -> InputRead:
        return InputRead.READ1

    def _iter(self, sam):
        for read in sam:
            if read.is_read1:
                yield self._as_sequence(read)


class Read2SingleEndSAMReader(SAMReader, PairedEndEstimatorMixin):
    """
    Reads a paired-end SAM/BAM file as if it were single-end, yielding only the
    second read from each pair.
    """

    @classproperty
    def interleaved(cls) -> bool:
        return False

    @property
    def input_read(self) -> InputRead:
        return InputRead.READ2

    def _iter(self, sam):
        for read in sam:
            if read.is_read2:
                yield self._as_sequence(read)


class PairedEndSAMReader(SAMReader, PairedEndEstimatorMixin):
    """
    Reads pairs of reads from a SAM/BAM file. The file must be name-sorted.
    """

    @classproperty
    def interleaved(cls) -> bool:
        return True

    @property
    def input_read(self) -> InputRead:
        return InputRead.PAIRED

    def _iter(self, sam):
        for reads in zip(sam, sam):
            if reads[0].query_name != reads[1].query_name:
                raise FormatError(
                    f"Consecutive reads {reads[0].query_name}, {reads[1].query_name} "
                    f"in paired-end SAM/BAM file do not have the same name; make sure "
                    f"your file is name-sorted and does not contain any "
                    f"secondary/supplementary alignments."
                )

            if not reads[0].is_read1 ^ reads[1].is_read1:
                raise ValueError(
                    f"There must be exactly one read with read1 flag set in pair "
                    f"{reads}"
                )

            if not reads[0].is_read2 ^ reads[1].is_read2:
                raise ValueError(
                    f"There must be exactly one read with read2 flag set in pair "
                    f"{reads}"
                )

            if reads[1].is_read1:
                # Swap reads so the first read is always read1
                reads = (reads[1], reads[0])

            yield tuple(self._as_sequence(r) for r in reads)


class FileWithPrependedLine:
    """
    A file-like object that allows to "prepend" a single line to an already
    opened file. That is, further reads on the file will return the provided
    line and only then the actual content. This is needed to solve the problem
    of autodetecting input from a stream: As soon as the first line has been
    read, we know the file type, but also that line is "gone" and unavailable
    for further processing.

    Args:
        file: An already opened file-like object.
        line: A single string (newline will be appended if not included).
    """

    def __init__(self, file: IO, line: str):
        if not line.endswith("\n"):
            line += "\n"
        self.first_line = line
        self._file = file

    @property
    def name(self) -> str:
        return self._file.name

    def __iter__(self):
        yield self.first_line
        yield from self._file

    def close(self):
        """
        Closes the underlying file.
        """
        self._file.close()


def open_reader(
    file1=None,
    file2=None,
    ngstream_reader=None,
    qualfile: Union[str, Path, Iterator[str]] = None,
    quality_base: Optional[int] = None,
    colorspace: bool = False,
    file_format: Optional[str] = None,
    interleaved: bool = False,
    input_read: Optional[InputRead] = None,
    alphabet: Optional[Alphabet] = None,
) -> SequenceReaderBase:
    """
    Open sequence files in FASTA or FASTQ format for reading. This is a factory that
    returns an instance of one of the ...Reader classes also defined in this module.

    Args:
        file1: Path to read1 regular or compressed file or file-like object.
        file2: Path to read2 regular or compressed file or file-like object.
        ngstream_reader: An ngstream.Protocol object.
        qualfile: Path to qualfile regular or compressed file or file-like object. If
            specified, then file1 must be a FASTA file and sequences are single-end.
            One of file2 and qualfile must always be None (no paired-end data is
            supported when reading qualfiles).
        quality_base: Base for quality values.
        colorspace: If True, instances of the Colorspace... classes
            are returned.
        file_format: If set to None, file format is autodetected from the file
            name extension. Set to 'fasta', 'fastq', 'sra-fastq', 'sam', or
            'bam' to not auto-detect. Colorspace is not auto-detected and must
            always be requested explicitly.
        interleaved: If True, then file1 contains interleaved paired-end data.
            file2 and qualfile must be None in this case.
        input_read: When file1 is a paired-end interleaved or SAM/BAM
            file, this specifies whether to only use the first or second read
            (1 or 2) or to use both reads (3).
        alphabet: An Alphabet instance - the alphabet to use to validate sequences.
    """
    if interleaved and (file2 is not None or qualfile is not None):
        raise ValueError("When interleaved is set, file2 and qualfile must be None")

    if file2 is not None and qualfile is not None:
        raise ValueError("Setting both file2 and qualfile is not supported")

    if alphabet and isinstance(alphabet, str):
        if alphabet not in ALPHABETS:
            raise ValueError(f"Invalid alphabet {alphabet}")

        alphabet = ALPHABETS[alphabet]

    if file2 is not None:
        return PairedSequenceReader(
            file1,
            file2,
            quality_base=quality_base,
            colorspace=colorspace,
            file_format=file_format,
            alphabet=alphabet,
        )

    if qualfile is not None:
        # read from .(CS)FASTA/.QUAL
        reader_class = ColorspaceFastaQualReader if colorspace else FastaQualReader
        return reader_class(
            file1,
            qualfile=qualfile,
            quality_base=quality_base,
            alphabet=alphabet
        )

    reader = None

    if file1:
        if file_format is None:
            if file1 not in (STDIN, sys.stdin):
                file_format = guess_format_from_name(file1)

            if file_format is None:
                if file1 == STDIN:
                    file1 = sys.stdin

                for line in file1:
                    if line.startswith("#"):
                        # Skip comment lines (needed for csfasta)
                        continue

                    if line.startswith(">"):
                        file_format = "fasta"
                    elif line.startswith("@"):
                        file_format = "fastq"

                    # TODO: guess SAM/BAM from data
                    file1 = FileWithPrependedLine(file1, line)

                    break

        if file_format is not None:
            file_format = file_format.lower()

            if file_format in ("sam", "bam"):
                if colorspace:
                    raise ValueError(
                        "SAM/BAM format is not currently supported for colorspace reads"
                    )

                if interleaved:
                    return PairedEndSAMReader(
                        file1, quality_base=quality_base, alphabet=alphabet
                    )
                elif input_read == InputRead.READ1:
                    return Read1SingleEndSAMReader(
                        file1, quality_base=quality_base, alphabet=alphabet
                    )
                elif input_read == InputRead.READ2:
                    return Read2SingleEndSAMReader(
                        file1, quality_base=quality_base, alphabet=alphabet
                    )
                else:
                    return SingleEndSAMReader(
                        file1, quality_base=quality_base, alphabet=alphabet
                    )

            if interleaved:
                reader = InterleavedSequenceReader(
                    file1,
                    quality_base=quality_base,
                    colorspace=colorspace,
                    file_format=file_format,
                    alphabet=alphabet
                )
            elif file_format == "fasta":
                reader_class = ColorspaceFastaReader if colorspace else FastaReader
                reader = reader_class(file1, alphabet=alphabet)
            elif file_format == "fastq":
                reader_class = ColorspaceFastqReader if colorspace else FastqReader
                reader = reader_class(
                    file1, quality_base=quality_base, alphabet=alphabet
                )
            elif file_format == "sra-fastq" and colorspace:
                def sra_colorspace_sequence_factory(
                    name: str,
                    sequence: str,
                    qualities: str,
                    *args,
                    **kwargs
                ):
                    """
                    Factory for an SRA colorspace sequence (which has one quality
                    value too many).
                    """
                    return ColorspaceSequence(
                        name, sequence, qualities[1:], *args, **kwargs
                    )

                reader = ColorspaceFastqReader(
                    file1,
                    quality_base=quality_base,
                    alphabet=alphabet,
                    sequence_factory=sra_colorspace_sequence_factory
                )
    elif ngstream_reader:
        if colorspace:
            reader_class = ColorspaceNgstreamSequenceReader
        else:
            reader_class = NgstreamSequenceReader

        reader = reader_class(
            ngstream_reader,
            quality_base=quality_base,
            alphabet=alphabet
        )

    if reader:
        if input_read == InputRead.PAIRED or not reader.paired:
            return reader
        else:
            return PairedToSingleEndReader(reader, input_read)

    raise UnknownFileTypeError(
        f"File format {file_format!r} is unknown (expected 'sra-fastq' (only for "
        f"colorspace), 'fasta', 'fastq', 'sam', or 'bam')."
    )
