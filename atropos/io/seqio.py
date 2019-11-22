# coding: utf-8
"""Sequence I/O classes: Reading and writing of FASTA and FASTQ files.

TODO
- Sequence.name should be Sequence.description or so (reserve .name for the part
  before the first space)
"""
import sys
from atropos import AtroposError
from atropos.io import STDOUT, xopen
from atropos.io.compression import splitext_compressed
from atropos.util import Summarizable, truncate_string, ALPHABETS

SINGLE = 0
READ1 = 1
READ2 = 2
PAIRED = 1|2

class FormatError(AtroposError):
    """Raised when an input file (FASTA or FASTQ) is malformatted."""
    pass

class UnknownFileType(AtroposError):
    """Raised when open could not autodetect the file type."""
    pass

## Reading sequences from files ##

class SequenceReaderBase(Summarizable): # pylint: disable=no-member
    """Sequence readers must provide the following properties:
    - input_names: (read1 file, read2 file)
    - input_read: 1, 2, or 3
    - file_format: string, e.g. FASTA, FASTQ, etc.
    - delivers_qualities: bool
    - has_qualfile: bool
    - quality_base: int
    - colorspace: bool
    - interleaved: bool
    """
    def summarize(self):
        return dict(
            input_names=self.input_names,
            input_read=self.input_read,
            file_format=self.file_format,
            delivers_qualities=self.delivers_qualities,
            quality_base=self.quality_base,
            has_qualfile=self.has_qualfile,
            colorspace=self.colorspace,
            interleaved=self.interleaved)

class SequenceReader(SequenceReaderBase):
    """Read possibly compressed files containing sequences.

    Args:
        file is a path or a file-like object. In both cases, the file may
            be compressed (.gz, .bz2, .xz).
        mode: The file open mode.
        quality_base: The minimum quality value.
        alphabet: The alphabet to use to validate sequences. If None, no
            validation is done.
    """
    delivers_qualities = False
    has_qualfile = False
    colorspace = False
    interleaved = False
    input_read = SINGLE
    _close_on_exit = False

    def __init__(self, path, mode='r', quality_base=None, alphabet=None):
        self.quality_base = quality_base
        self.alphabet = alphabet
        if isinstance(path, str):
            self.name = path
            self._file = xopen(path, mode)
            self._close_on_exit = True
        else:
            if hasattr(path, 'name'):
                self.name = path.name
            else:
                # TODO: generate random unique name?
                self.name = path.__class__
            self._file = path

    @property
    def input_names(self):
        return (self.name, None)

    def close(self):
        """Close the underlying file.
        """
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        if self._file is None:
            raise ValueError("I/O operation on closed SequenceReader")
        return self

    def __exit__(self, *args):
        self.close()

try:
    from ._seqio import Sequence, FastqReader
except ImportError:
    pass

class ColorspaceSequence(Sequence):
    """Sequence object for colorspace reads.
    """
    def __init__(
            self, name, sequence, qualities, primer=None, name2='',
            original_length=None, match=None, match_info=None, clipped=None,
            insert_overlap=False, merged=False, corrected=0, alphabet=None):
        # In colorspace, the first character is the last nucleotide of the
        # primer base and the second character encodes the transition from the
        # primer base to the first real base of the read.
        if primer is None:
            self.primer = sequence[0:1]
            sequence = sequence[1:]
        else:
            self.primer = primer
        if qualities is not None and len(sequence) != len(qualities):
            rname = truncate_string(name)
            raise FormatError(
                "In read named {0!r}: length of colorspace quality "
                "sequence ({1}) and length of read ({2}) do not match (primer "
                "is: {3!r})".format(
                    rname, len(qualities), len(sequence), self.primer))
        super().__init__(
            name, sequence, qualities, name2, original_length, match,
            match_info, clipped, insert_overlap, merged, corrected,
            alphabet=alphabet)
        # TODO: use 'alphabet' here
        if not self.primer in ('A', 'C', 'G', 'T'):
            raise FormatError(
                "Primer base is {0!r} in read {1!r}, but it should be one of "
                "A, C, G, T.".format(self.primer, truncate_string(name)))

    def __repr__(self):
        fmt_str = \
            '<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>'
        qstr = ''
        if self.qualities is not None:
            qstr = ', qualities={0!r}'.format(truncate_string(self.qualities))
        return fmt_str.format(
            truncate_string(self.name), self.primer,
            truncate_string(self.sequence), qstr)

    def __getitem__(self, key):
        return self.__class__(
            self.name,
            self.sequence[key],
            self.qualities[key] if self.qualities is not None else None,
            self.primer,
            self.name2,
            self.original_length,
            self.match,
            self.match_info,
            self.clipped,
            self.insert_overlap,
            self.merged,
            self.corrected)

class SraSequenceReader(SequenceReader):
    delivers_qualities = True
    file_format = "fastq"

    def __init__(
            self, reader, quality_base=None, sequence_class=Sequence,
            alphabet=None):
        super().__init__(reader, quality_base=quality_base, alphabet=alphabet)
        self.input_read = PAIRED if reader.paired else SINGLE
        self.sequence_class = sequence_class

    def __iter__(self):
        if self.input_read == PAIRED:
            for read in self._file:
                yield tuple(self._as_sequence(frag) for frag in read[:2])
        else:
            for read in self._file:
                yield self._as_sequence(read[0])

    def _as_sequence(self, frag):
        return self.sequence_class(*frag, alphabet=self.alphabet)

    def close(self):
        self._file.finish()

class SraColorspaceSequenceReader(SraSequenceReader):
    """Reads colorspace sequences from an SRA accession.
    """
    colorspace = True

    def __init__(self, reader, quality_base=33, alphabet=None):
        super().__init__(
            reader, quality_base=quality_base,
            sequence_class=ColorspaceSequence,
            alphabet=alphabet)

class FileWithPrependedLine(object):
    """A file-like object that allows to "prepend" a single line to an already
    opened file. That is, further reads on the file will return the provided
    line and only then the actual content. This is needed to solve the problem
    of autodetecting input from a stream: As soon as the first line has been
    read, we know the file type, but also that line is "gone" and unavailable
    for further processing.

    Args:
        file: An already opened file-like object.
        line: A single string (newline will be appended if not included).
    """
    def __init__(self, file, line):
        if not line.endswith('\n'):
            line += '\n'
        self.first_line = line
        self._file = file

    @property
    def name(self):
        return self._file.name

    def __iter__(self):
        yield self.first_line
        for line in self._file:
            yield line

    def close(self):
        """Close the underlying file.
        """
        self._file.close()

class FastaReader(SequenceReader):
    """Reader for FASTA files.

    Args:
        path: A path or a file-like object. In both cases, the file may
            be compressed (.gz, .bz2, .xz).
        keep_linebreaks: Whether to keep newline characters in the sequence.
        sequence_class: The class to use when creating new sequence objects.
    """
    file_format = "FASTA"

    def __init__(
            self, path, keep_linebreaks=False, sequence_class=Sequence,
            alphabet=None):
        super().__init__(path, alphabet=alphabet)
        self.sequence_class = sequence_class
        self._delimiter = '\n' if keep_linebreaks else ''

    def __iter__(self):
        """Read next entry from the file (single entry at a time).
        """
        name = None
        seq = []
        for i, line in enumerate(self._file):
            # strip() also removes DOS line breaks
            line = line.strip()
            if not line:
                continue
            if line and line[0] == '>':
                if name is not None:
                    yield self.sequence_class(
                        name, self._delimiter.join(seq), None,
                        alphabet=self.alphabet)
                name = line[1:]
                seq = []
            elif line and line[0] == '#':
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FormatError(
                    "At line {0}: Expected '>' at beginning of FASTA record, "
                    "but got {1!r}.".format(i+1, truncate_string(line)))

        if name is not None:
            yield self.sequence_class(
                name, self._delimiter.join(seq), None,
                alphabet=self.alphabet)

class ColorspaceFastaReader(FastaReader):
    """Reads colorspace sequences from a FASTA.

    Args:
        path: The file to read.
        keep_linebreaks: Whether to keep linebreaks in wrapped sequences.
    """
    colorspace = True

    def __init__(self, path, keep_linebreaks=False, alphabet=None):
        super().__init__(
            path, keep_linebreaks, sequence_class=ColorspaceSequence,
            alphabet=alphabet)

class ColorspaceFastqReader(FastqReader):
    """Reads colorspace sequences from a FASTQ.
    """
    colorspace = True

    def __init__(self, path, quality_base=33, alphabet=None):
        super().__init__(
            path, quality_base=quality_base, sequence_class=ColorspaceSequence,
            alphabet=alphabet)

class SRAColorspaceFastqReader(FastqReader):
    """Reads SRA-formatted colorspace sequences from a FASTQ.
    """
    colorspace = True

    def __init__(self, path, quality_base=33, alphabet=None):
        super().__init__(
            path, quality_base=quality_base,
            sequence_class=sra_colorspace_sequence,
            alphabet=alphabet)

class FastaQualReader(SequenceReaderBase):
    """Reader for reads that are stored in .(CS)FASTA and .QUAL files.

    Args:
        fastafile and qualfile are filenames or file-like objects.
            If a filename is used, then .gz files are recognized.
        sequence_class: The class to use when creating new sequence objects.
    """
    file_format = "FastaQual"
    delivers_qualities = True
    has_qualfile = True
    colorspace = False
    interleaved = False
    input_read = SINGLE

    def __init__(
            self, fastafile, qualfile, quality_base=33,
            sequence_class=Sequence, alphabet=None):
        self.fastareader = FastaReader(fastafile)
        self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
        self.quality_base = quality_base
        self.sequence_class = sequence_class
        self.alphabet = alphabet

    @property
    def input_names(self):
        return ((self.fastareader.name, self.qualreader.name), None)

    def __iter__(self):
        """Yield Sequence objects.
        """
        # conversion dictionary: maps strings to the appropriate ASCII-encoded
        # character
        conv = dict()
        for i in range(-5, 256 - 33):
            conv[str(i)] = chr(i + 33)
        for fastaread, qualread in zip(self.fastareader, self.qualreader):
            if fastaread.name != qualread.name:
                raise FormatError(
                    "The read names in the FASTA and QUAL file do not match "
                    "({0!r} != {1!r})".format(fastaread.name, qualread.name))
            try:
                qualities = ''.join(
                    [conv[value] for value in qualread.sequence.split()])
            except KeyError as err:
                raise FormatError(
                    "Within read named {0!r}: Found invalid quality "
                    "value {1}".format(fastaread.name, err))
            assert fastaread.name == qualread.name
            yield self.sequence_class(
                fastaread.name, fastaread.sequence, qualities,
                alphabet=self.alphabet)

    def close(self):
        """Close the underlying files.
        """
        self.fastareader.close()
        self.qualreader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

class ColorspaceFastaQualReader(FastaQualReader):
    """Reads sequences and qualities from separate files and returns
    :class:`ColorspaceSequence`s.

    Args:
        fastafile, qualfile: FASTA files that contain the sequences and
            qualities, respectively.
    """
    colorspace = True

    def __init__(self, fastafile, qualfile, quality_base=33, alphabet=None):
        super().__init__(
            fastafile, qualfile, quality_base=quality_base,
            sequence_class=ColorspaceSequence, alphabet=alphabet)

class PairedSequenceReader(SequenceReaderBase):
    """Read paired-end reads from two files. Wraps two SequenceReader
    instances, making sure that reads are properly paired.

    Args:
        file1, file2: The pair of files.
        colorspace: Whether the sequences are in colorspace.
        file_format: A file_format instance.
    """
    input_read = PAIRED
    interleaved = False

    def __init__(
            self, file1, file2, quality_base=33, colorspace=False,
            file_format=None, alphabet=None):
        self.reader1 = open_reader(
            file1, colorspace=colorspace, quality_base=quality_base,
            file_format=file_format, alphabet=alphabet)
        self.reader2 = open_reader(
            file2, colorspace=colorspace, quality_base=quality_base,
            file_format=file_format, alphabet=alphabet)

    @property
    def input_names(self):
        return (
            self.reader1.input_names[0],
            self.reader2.input_names[0])

    def __getattr__(self, name):
        return getattr(self.reader1, name)

    def __iter__(self):
        """Iterate over the paired reads. Each item is a pair of Sequence
        objects.
        """
        # Avoid usage of zip() below since it will consume one item too many.
        it1, it2 = iter(self.reader1), iter(self.reader2)
        while True:
            try:
                read1 = next(it1)
            except StopIteration:
                # End of file 1. Make sure that file 2 is also at end.
                try:
                    next(it2)
                    raise FormatError(
                        "Reads are improperly paired. There are more reads in "
                        "file 2 than in file 1.")
                except StopIteration:
                    pass
                break
            try:
                read2 = next(it2)
            except StopIteration:
                raise FormatError(
                    "Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.")
            if not sequence_names_match(read1, read2):
                raise FormatError(
                    "Reads are improperly paired. Read name '{0}' in file 1 "
                    "does not match '{1}' in file 2.".format(
                        read1.name, read2.name))
            yield (read1, read2)

    def close(self):
        """Close the underlying files.
        """
        self.reader1.close()
        self.reader2.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

class InterleavedSequenceReader(SequenceReaderBase):
    """Read paired-end reads from an interleaved FASTQ file.

    Args:
        path: The interleaved FASTQ file.
        colorspace: Whether the sequences are in colorspace.
        file_format: A file_format instance.
    """
    input_read = PAIRED
    interleaved = True

    def __init__(
            self, path, quality_base=33, colorspace=False, file_format=None,
            alphabet=None):
        self.reader = open_reader(
            path, quality_base=quality_base, colorspace=colorspace,
            file_format=file_format, alphabet=alphabet)

    def __getattr__(self, name):
        return getattr(self.reader, name)

    def __iter__(self):
        # Avoid usage of zip() below since it will consume one item too many.
        itr = iter(self.reader)
        for read1 in itr:
            try:
                read2 = next(itr)
            except StopIteration:
                raise FormatError(
                    "Interleaved input file incomplete: Last record has no "
                    "partner.")
            if not sequence_names_match(read1, read2):
                raise FormatError(
                    "Reads are improperly paired. Name {0!r} (first) does not "
                    "match {1!r} (second).".format(read1.name, read2.name))
            yield (read1, read2)

    def close(self):
        """Close the underlying reader.
        """
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

# TODO: SAM/BAM classes need unit tests

class SAMReader(SequenceReaderBase):
    """Reader for SAM/BAM files. Paired-end files must be name-sorted. Does
    not support secondary/supplementary reads. This is an abstract class.

    Args:
        path: A filename or a file-like object. If a filename, then .gz files
            are supported.
        sequence_class: The class to use when creating new sequence objects.
    """
    file_format = "SAM"
    delivers_qualities = True
    interleaved = False
    has_qualfile = False
    colorspace = False

    def __init__(
            self, path, quality_base=33, sequence_class=Sequence,
            alphabet=None, pysam_kwargs=None):
        self._close_on_exit = False
        if isinstance(path, str):
            path = xopen(path, 'rb')
            self._close_on_exit = True
        if isinstance(path, str):
            self.name = path
        else:
            self.name = path.name
        self._file = path
        self.quality_base = quality_base
        self.sequence_class = sequence_class
        self.alphabet = alphabet
        self.pysam_kwargs = pysam_kwargs or dict(check_sq=False)

    @property
    def input_names(self):
        return (self.name, None)

    def __iter__(self):
        import pysam
        return self._iter(pysam.AlignmentFile(self._file, **self.pysam_kwargs))

    def _iter(self, sam):
        """Create an iterator over records in the SAM/BAM file.
        """
        raise NotImplementedError()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """Close the underling AlignmentFile.
        """
        if self._close_on_exit and self._file is not None:
            self._file.close()
            self._file = None

    def _as_sequence(self, read):
        return self.sequence_class(
            read.query_name,
            read.query_sequence,
            ''.join(chr(33 + q) for q in read.query_qualities),
            alphabet=self.alphabet)

class SingleEndSAMReader(SAMReader):
    """Reader for single-end SAM/BAM files.
    """
    input_read = SINGLE

    def _iter(self, sam):
        for read in sam:
            yield self._as_sequence(read)

class Read1SingleEndSAMReader(SAMReader):
    """Reads a paired-end SAM/BAM file as if it were single-end, yielding
    only the first read from each pair.
    """
    input_read = READ1

    def _iter(self, sam):
        for read in sam:
            if read.is_read1:
                yield self._as_sequence(read)

class Read2SingleEndSAMReader(SAMReader):
    """Reads a paired-end SAM/BAM file as if it were single-end, yielding
    only the second read from each pair.
    """
    input_read = READ2

    def _iter(self, sam):
        for read in sam:
            if read.is_read2:
                yield self._as_sequence(read)

class PairedEndSAMReader(SAMReader):
    """Reads pairs of reads from a SAM/BAM file. The file must be name-sorted.
    """
    input_read = PAIRED
    interleaved = True

    def _iter(self, sam):
        for reads in zip(sam, sam):
            if reads[0].query_name != reads[1].query_name:
                raise AtroposError(
                    "Consecutive reads {}, {} in paired-end SAM/BAM file do "
                    "not have the same name; make sure your file is "
                    "name-sorted and does not contain any "
                    "secondary/supplementary alignments.",
                    reads[0].query_name, reads[1].query_name)

            if reads[0].is_read1:
                assert reads[1].is_read2
            else:
                assert reads[1].is_read1
                reads = (reads[1], reads[0])

            yield tuple(self._as_sequence(r) for r in reads)

class SequenceFileFormat():
    """Base class for sequence formatters.
    """
    def format(self, read):
        """Format a Sequence as a string.

        Args:
            read: The Sequence object.

        Returns:
            A string representation of the sequence object in the sequence
            file format.
        """
        raise NotImplementedError()

class FastaFormat(SequenceFileFormat):
    """FASTA SequenceFileFormat.

    Args:
        line_length: Max line length (in characters), or None. Determines
            whether and how lines are wrapped.
    """
    def __init__(self, line_length=None):
        self.text_wrapper = None
        if line_length:
            from textwrap import TextWrapper
            self.text_wrapper = TextWrapper(width=line_length)

    def format(self, read):
        return self.format_entry(read.name, read.sequence)

    def format_entry(self, name, sequence):
        """Convert a sequence record to a string.
        """
        if self.text_wrapper:
            sequence = self.text_wrapper.fill(sequence)
        return "".join((">", name, "\n", sequence, "\n"))

class ColorspaceFastaFormat(FastaFormat):
    """FastaFormat in which sequences are in colorspace.
    """
    def format(self, read):
        return self.format_entry(read.name, read.primer + read.sequence)

class FastqFormat(SequenceFileFormat):
    """FASTQ SequenceFileFormat.
    """
    def format(self, read):
        return self.format_entry(
            read.name, read.sequence, read.qualities, read.name2)

    def format_entry(self, name, sequence, qualities, name2=""):
        """Convert a sequence record to a string.
        """
        return "".join((
            '@', name, '\n',
            sequence, '\n+',
            name2, '\n',
            qualities, '\n'))

class ColorspaceFastqFormat(FastqFormat):
    """FastqFormat in which sequences are in colorspace.
    """
    def format(self, read):
        return self.format_entry(
            read.name, read.primer + read.sequence, read.qualities)

class SingleEndFormatter():
    """Wrapper for a SequenceFileFormat for single-end data.

    Args:
        seq_format: The SequenceFileFormat object.
        file1: The single-end file.
    """
    def __init__(self, seq_format, file1):
        self.seq_format = seq_format
        self.file1 = file1
        self.written = 0
        self.read1_bp = 0
        self.read2_bp = 0

    def format(self, result, read1, read2=None):
        """Format read(s) and add them to `result`.

        Args:
            result: A dict mapping file names to lists of formatted reads.
            read1, read2: The reads to format.
        """
        result[self.file1].append(self.seq_format.format(read1))
        self.written += 1
        self.read1_bp += len(read1)

    @property
    def written_bp(self):
        """Tuple of base-pairs written (read1_bp, read2_bp).
        """
        return (self.read1_bp, self.read2_bp)

class InterleavedFormatter(SingleEndFormatter):
    """Format read pairs as successive reads in an interleaved file.
    """
    def format(self, result, read1, read2=None):
        result[self.file1].extend((
            self.seq_format.format(read1),
            self.seq_format.format(read2)))
        self.written += 1
        self.read1_bp += len(read1)
        self.read2_bp += len(read2)

class PairedEndFormatter(SingleEndFormatter):
    """Wrapper for a SequenceFileFormat. Both reads in a pair are formatted
    using the specified format.
    """
    def __init__(self, seq_format, file1, file2):
        super(PairedEndFormatter, self).__init__(seq_format, file1)
        self.file2 = file2

    def format(self, result, read1, read2):
        result[self.file1].append(self.seq_format.format(read1))
        result[self.file2].append(self.seq_format.format(read2))
        self.written += 1
        self.read1_bp += len(read1)
        self.read2_bp += len(read2)

def sra_colorspace_sequence(name, sequence, qualities, name2, alphabet=None):
    """Factory for an SRA colorspace sequence (which has one quality value
    too many).
    """
    return ColorspaceSequence(
        name, sequence, qualities[1:], name2=name2, alphabet=alphabet)

def sequence_names_match(read1, read2):
    """Check whether the sequences read1 and read2 have identical names,
    ignoring a suffix of '1' or '2'. Some old paired-end reads have names that
    end in '/1' and '/2'. Also, the fastq-dump tool (used for converting SRA
    files to FASTQ) appends a .1 and .2 to paired-end reads if option -I is
    used.

    Args:
        read1, read2: The sequences to compare.

    Returns:
        Whether the sequences are equal.
    """
    name1 = read1.name.split(None, 1)[0]
    name2 = read2.name.split(None, 1)[0]
    if name1[-1:] in '12' and name2[-1:] in '12':
        name1 = name1[:-1]
        name2 = name2[:-1]
    return name1 == name2

def paired_to_read1(reader):
    """Generator that yields the first read from an iterator over read pairs.
    """
    for read1, _ in reader:
        yield read1

def paired_to_read2(reader):
    """Generator that yields the second read from an iterator over read pairs.
    """
    for _, read2 in reader:
        yield read2

def open_reader(
        file1=None, file2=None, qualfile=None, quality_base=None,
        colorspace=False, file_format=None, interleaved=False,
        input_read=None, alphabet=None):
    """Open sequence files in FASTA or FASTQ format for reading. This is
    a factory that returns an instance of one of the ...Reader
    classes also defined in this module.

    Args:
        file1, file2, qualfile: Paths to regular or compressed files or
            file-like objects. Use file1 if data is single-end. If file2 is also
            provided, sequences are paired. If qualfile is given, then file1
            must be a FASTA file and sequences are single-end. One of file2 and
            qualfile must always be None (no paired-end data is supported when
            reading qualfiles).
        quality_base: Base for quality values.
        interleaved:If True, then file1 contains interleaved paired-end data.
            file2 and qualfile must be None in this case.
        colorspace: If True, instances of the Colorspace... classes
            are returned.
        file_format: If set to None, file format is autodetected from the file
            name extension. Set to 'fasta', 'fastq', 'sra-fastq', 'sam', or
            'bam' to not auto-detect. Colorspace is not auto-detected and must
            always be requested explicitly.
        input_read: When file1 is a paired-end interleaved or SAM/BAM
            file, this specifies whether to only use the first or second read
            (1 or 2) or to use both reads (None).
        alphabet: An Alphabet instance - the alphabet to use to validate
            sequences.
    """
    if interleaved and (file2 is not None or qualfile is not None):
        raise ValueError(
            "When interleaved is set, file2 and qualfile must be None")
    if file2 is not None and qualfile is not None:
        raise ValueError("Setting both file2 and qualfile is not supported")

    if alphabet and isinstance(alphabet, str):
        if not alphabet in ALPHABETS:
            raise ValueError("Invalid alphabet {}".format(alphabet))
        alphabet = ALPHABETS[alphabet]

    if file2 is not None:
        return PairedSequenceReader(
            file1, file2, quality_base=quality_base, colorspace=colorspace,
            file_format=file_format, alphabet=alphabet)

    if qualfile is not None:
        if colorspace:
            # read from .(CS)FASTA/.QUAL
            return ColorspaceFastaQualReader(
                file1, qualfile, quality_base=quality_base, alphabet=alphabet)
        else:
            return FastaQualReader(
                file1, qualfile, quality_base=quality_base, alphabet=alphabet)

    if file_format is None and file1 != STDOUT:
        file_format = guess_format_from_name(file1)

    if file_format is None:
        if file1 == STDOUT:
            file1 = sys.stdin
        for line in file1:
            if line.startswith('#'):
                # Skip comment lines (needed for csfasta)
                continue
            if line.startswith('>'):
                file_format = 'fasta'
            elif line.startswith('@'):
                file_format = 'fastq'
            # TODO: guess SAM/BAM from data
            file1 = FileWithPrependedLine(file1, line)
            break

    if file_format is not None:
        file_format = file_format.lower()
        if file_format in ("sam", "bam"):
            if colorspace:
                raise ValueError(
                    "SAM/BAM format is not currently supported for colorspace "
                    "reads")
            if interleaved:
                return PairedEndSAMReader(
                    file1, quality_base=quality_base, alphabet=alphabet)
            elif input_read == READ1:
                return Read1SingleEndSAMReader(
                    file1, quality_base=quality_base, alphabet=alphabet)
            elif input_read == READ2:
                return Read2SingleEndSAMReader(
                    file1, quality_base=quality_base, alphabet=alphabet)
            else:
                return SingleEndSAMReader(
                    file1, quality_base=quality_base, alphabet=alphabet)
        elif interleaved:
            reader = InterleavedSequenceReader(
                file1, quality_base=quality_base, colorspace=colorspace,
                file_format=file_format, alphabet=alphabet)
            if input_read == READ1:
                return paired_to_read1(reader)
            elif input_read == READ2:
                return paired_to_read2(reader)
            else:
                return reader
        elif file_format == 'fasta':
            fasta_handler = ColorspaceFastaReader if colorspace else FastaReader
            return fasta_handler(file1, alphabet=alphabet)
        elif file_format == 'fastq':
            fastq_handler = ColorspaceFastqReader if colorspace else FastqReader
            return fastq_handler(
                file1, quality_base=quality_base, alphabet=alphabet)
        elif file_format == 'sra-fastq' and colorspace:
            return SRAColorspaceFastqReader(
                file1, quality_base=quality_base, alphabet=alphabet)

    raise UnknownFileType(
        "File format {0!r} is unknown (expected 'sra-fastq' (only for "
        "colorspace), 'fasta', 'fastq', 'sam', or 'bam').".format(
            file_format or '<Undetected>'
        ))

def sra_reader(
        reader, quality_base=None, colorspace=False, input_read=None,
        alphabet=None):
    """Wrap an existing SraReader. The reader must 1) have a 'paired' property,
    and 2) be iterable. Furthermore, each value yielded by the iterator must
    be a list of N tuples, where N = 2 if paired else 1, and where each tuple
    is (name, sequence, qualities).

    Args:
        reader: An existing reader.
        quality_base: Base for quality values.
        colorspace: If True, instances of the Colorspace... classes
            are returned.
        input_read: When file1 is a paired-end interleaved or SAM/BAM
            file, this specifies whether to only use the first or second read
            (1 or 2) or to use both reads (None).
        alphabet: An Alphabet instance - the alphabet to use to validate
            sequences.
    """
    if colorspace:
        wrapped = SraColorspaceSequenceReader(
            reader, quality_base=quality_base, alphabet=alphabet)
    else:
        wrapped = SraSequenceReader(
            reader, quality_base=quality_base, alphabet=alphabet)

    if not reader.paired or input_read == PAIRED:
        return wrapped

    if input_read == READ1:
        return paired_to_read1(wrapped)
    else:
        return paired_to_read2(wrapped)

def guess_format_from_name(path, raise_on_failure=False):
    """Detect file format based on the file name.

    Args:
        path: The filename to guess.
        raise_on_failure: Whether to raise an exception if the filename cannot
            be detected.

    Returns:
        The format name.
    """
    name = None
    if isinstance(path, str):
        name = path
    elif hasattr(path, "name"):	 # seems to be an open file-like object
        name = path.name

    if name:
        name, ext1, _ = splitext_compressed(name)
        ext = ext1.lower()
        if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
            return 'fasta'
        elif ext in ['.fastq', '.fq'] or (
                ext == '.txt' and name.endswith('_sequence')):
            return 'fastq'
        elif ext in ('.sam', '.bam'):
            return ext[1:]

    if raise_on_failure:
        raise UnknownFileType(
            "Could not determine whether file {0!r} is FASTA or FASTQ: file "
            "name extension {1!r} not recognized".format(path, ext))

def create_seq_formatter(file1, file2=None, interleaved=False, **kwargs):
    """Create a formatter, deriving the format name from the file extension.

    Args:
        file1, file2: Output files.
        interleaved: Whether the output should be interleaved (file2 must be
            None).
        kwargs: Additional arguments to pass to :method:`get_format`.
    """
    seq_format = get_format(file1, **kwargs)
    if file2 is not None:
        return PairedEndFormatter(seq_format, file1, file2)
    elif interleaved:
        return InterleavedFormatter(seq_format, file1)
    else:
        return SingleEndFormatter(seq_format, file1)

def get_format(
        path, file_format=None, colorspace=False, qualities=None,
        line_length=None):
    """Create a SequenceFileFormat instance.

    Args:
        path: The filename.

        file_format: If set to None, file format is autodetected from the file
        name extension. Set to 'fasta', 'fastq', or 'sra-fastq' to not
        auto-detect. Colorspace is not auto-detected and must always be
        requested explicitly.

        colorspace: If True, instances of the Colorspace... formats are
            returned.

        qualities: When file_format is None, this can be set to True or False to
            specify whether the written sequences will have quality values.
            This is is used in two ways:
            * If the output format cannot be determined (unrecognized extension
              etc), no exception is raised, but fasta or fastq format is chosen
              appropriately.
            * When False (no qualities available), an exception is raised when
              the auto-detected output format is FASTQ.

    Returns:
        A SequenceFileFormat object.
    """
    if file_format is None:
        file_format = guess_format_from_name(
            path, raise_on_failure=qualities is None)

    if file_format is None:
        if qualities is True:
            # Format not recognized, but know we want to write reads with
            # qualities.
            file_format = 'fastq'
        elif qualities is False:
            # Same, but we know that we want to write reads without qualities.
            file_format = 'fasta'
        else:
            raise UnknownFileType("Could not determine file type.")

    file_format = file_format.lower()

    if file_format == 'fastq' and qualities is False:
        raise ValueError(
            "Output format cannot be FASTQ since no quality values are "
            "available.")

    if file_format == 'fasta':
        if colorspace:
            return ColorspaceFastaFormat(line_length)
        else:
            return FastaFormat(line_length)
    elif file_format == 'fastq':
        if colorspace:
            return ColorspaceFastqFormat()
        else:
            return FastqFormat()
    else:
        raise UnknownFileType(
            "File format {0!r} is unknown (expected 'fasta' or "
            "'fastq').".format(file_format))
