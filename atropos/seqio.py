# coding: utf-8
"""
Sequence I/O classes: Reading and writing of FASTA and FASTQ files.

TODO

- Sequence.name should be Sequence.description or so (reserve .name for the part
  before the first space)
"""
import re
import sys
import io
from os.path import splitext
from .filters import NoFilter
from .xopen import xopen, open_output

def _shorten(s, n=100):
    """Shorten string s to at most n characters, appending "..." if necessary."""
    if s is None:
        return None
    if len(s) > n:
        s = s[:n-3] + '...'
    return s

class FormatError(Exception):
    """
    Raised when an input file (FASTA or FASTQ) is malformatted.
    """

class UnknownFileType(Exception):
    """
    Raised when open could not autodetect the file type.
    """

## Reading sequences from files ##

class SequenceReader(object):
    """Read possibly compressed files containing sequences"""
    _close_on_exit = False

    def __init__(self, file):
        """
        file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).
        """
        if isinstance(file, str):
            file = xopen(file)
            self._close_on_exit = True
        self._file = file

    def close(self):
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
    def __init__(self, name, sequence, qualities, primer=None, name2='', original_length=None, match=None, match_info=None, clipped=[0,0,0]):
        # In colorspace, the first character is the last nucleotide of the primer base
        # and the second character encodes the transition from the primer base to the
        # first real base of the read.
        if primer is None:
            self.primer = sequence[0:1]
            sequence = sequence[1:]
        else:
            self.primer = primer
        if qualities is not None and len(sequence) != len(qualities):
            rname = _shorten(name)
            raise FormatError("In read named {0!r}: length of colorspace quality "
                "sequence ({1}) and length of read ({2}) do not match (primer "
                "is: {3!r})".format(rname, len(qualities), len(sequence), self.primer))
        super(ColorspaceSequence, self).__init__(name, sequence, qualities, name2, original_length, match, match_info, clipped)
        if not self.primer in ('A', 'C', 'G', 'T'):
            raise FormatError("Primer base is {0!r} in read {1!r}, but it "
                "should be one of A, C, G, T.".format(
                    self.primer, _shorten(name)))

    def __repr__(self):
        qstr = ''
        if self.qualities is not None:
            qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
        return '<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>'.format(
            _shorten(self.name), self.primer, _shorten(self.sequence), qstr)

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
            self.clipped)

def sra_colorspace_sequence(name, sequence, qualities, name2):
    """Factory for an SRA colorspace sequence (which has one quality value too many)"""
    return ColorspaceSequence(name, sequence, qualities[1:], name2=name2)

class FileWithPrependedLine(object):
    """
    A file-like object that allows to "prepend" a single
    line to an already opened file. That is, further
    reads on the file will return the provided line and
    only then the actual content. This is needed to solve
    the problem of autodetecting input from a stream:
    As soon as the first line has been read, we know
    the file type, but also that line is "gone" and
    unavailable for further processing.
    """
    def __init__(self, file, line):
        """
        file is an already opened file-like object.
        line is a single string (newline will be appended if not included)
        """
        if not line.endswith('\n'):
            line += '\n'
        self.first_line = line
        self._file = file

    def __iter__(self):
        yield self.first_line
        for line in self._file:
            yield line

    def close(self):
        self._file.close()

class FastaReader(SequenceReader):
    """
    Reader for FASTA files.
    """
    def __init__(self, file, keep_linebreaks=False, sequence_class=Sequence):
        """
        file is a path or a file-like object. In both cases, the file may
        be compressed (.gz, .bz2, .xz).

        keep_linebreaks -- whether to keep newline characters in the sequence
        """
        super(FastaReader, self).__init__(file)
        self.sequence_class = sequence_class
        self.delivers_qualities = False
        self._delimiter = '\n' if keep_linebreaks else ''

    def __iter__(self):
        """
        Read next entry from the file (single entry at a time).
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
                    yield self.sequence_class(name, self._delimiter.join(seq), None)
                name = line[1:]
                seq = []
            elif line and line[0] == '#':
                continue
            elif name is not None:
                seq.append(line)
            else:
                raise FormatError("At line {0}: Expected '>' at beginning of "
                    "FASTA record, but got {1!r}.".format(i+1, _shorten(line)))

        if name is not None:
            yield self.sequence_class(name, self._delimiter.join(seq), None)

class ColorspaceFastaReader(FastaReader):
    def __init__(self, file, keep_linebreaks=False):
        super(ColorspaceFastaReader, self).__init__(file, keep_linebreaks,
            sequence_class=ColorspaceSequence)

class ColorspaceFastqReader(FastqReader):
    def __init__(self, file):
        super(ColorspaceFastqReader, self).__init__(file, sequence_class=ColorspaceSequence)

class SRAColorspaceFastqReader(FastqReader):
    def __init__(self, file):
        super(SRAColorspaceFastqReader, self).__init__(file, sequence_class=sra_colorspace_sequence)

class FastaQualReader(object):
    """
    Reader for reads that are stored in .(CS)FASTA and .QUAL files.
    """
    delivers_qualities = True

    def __init__(self, fastafile, qualfile, sequence_class=Sequence):
        """
        fastafile and qualfile are filenames or file-like objects.
        If a filename is used, then .gz files are recognized.

        The objects returned when iteritng over this file are instances of the
        given sequence_class.
        """
        self.fastareader = FastaReader(fastafile)
        self.qualreader = FastaReader(qualfile, keep_linebreaks=True)
        self.sequence_class = sequence_class

    def __iter__(self):
        """
        Yield Sequence objects.
        """
        # conversion dictionary: maps strings to the appropriate ASCII-encoded character
        conv = dict()
        for i in range(-5, 256 - 33):
            conv[str(i)] = chr(i + 33)
        for fastaread, qualread in zip(self.fastareader, self.qualreader):
            if fastaread.name != qualread.name:
                raise FormatError("The read names in the FASTA and QUAL file "
                    "do not match ({0!r} != {1!r})".format(fastaread.name, qualread.name))
            try:
                qualities = ''.join([conv[value] for value in qualread.sequence.split()])
            except KeyError as e:
                raise FormatError("Within read named {0!r}: Found invalid quality "
                    "value {1}".format(fastaread.name, e))
            assert fastaread.name == qualread.name
            yield self.sequence_class(fastaread.name, fastaread.sequence, qualities)

    def close(self):
        self.fastareader.close()
        self.qualreader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

class ColorspaceFastaQualReader(FastaQualReader):
    def __init__(self, fastafile, qualfile):
        super(ColorspaceFastaQualReader, self).__init__(fastafile, qualfile,
            sequence_class=ColorspaceSequence)

def sequence_names_match(r1, r2):
    """
    Check whether the sequences r1 and r2 have identical names, ignoring a
    suffix of '1' or '2'. Some old paired-end reads have names that end in '/1'
    and '/2'. Also, the fastq-dump tool (used for converting SRA files to FASTQ)
    appends a .1 and .2 to paired-end reads if option -I is used.
    """
    name1 = r1.name.split(None, 1)[0]
    name2 = r2.name.split(None, 1)[0]
    if name1[-1:] in '12' and name2[-1:] in '12':
        name1 = name1[:-1]
        name2 = name2[:-1]
    return name1 == name2

class PairedSequenceReader(object):
    """
    Read paired-end reads from two files.

    Wraps two SequenceReader instances, making sure that reads are properly
    paired.
    """
    def __init__(self, file1, file2, colorspace=False, fileformat=None):
        self.reader1 = open_reader(file1, colorspace=colorspace, fileformat=fileformat)
        self.reader2 = open_reader(file2, colorspace=colorspace, fileformat=fileformat)
        self.delivers_qualities = self.reader1.delivers_qualities

    def __iter__(self):
        """
        Iterate over the paired reads. Each item is a pair of Sequence objects.
        """
        # Avoid usage of zip() below since it will consume one item too many.
        it1, it2 = iter(self.reader1), iter(self.reader2)
        while True:
            try:
                r1 = next(it1)
            except StopIteration:
                # End of file 1. Make sure that file 2 is also at end.
                try:
                    next(it2)
                    raise FormatError("Reads are improperly paired. There are more reads in "
                        "file 2 than in file 1.")
                except StopIteration:
                    pass
                break
            try:
                r2 = next(it2)
            except StopIteration:
                raise FormatError("Reads are improperly paired. There are more reads in "
                    "file 1 than in file 2.")
            if not sequence_names_match(r1, r2):
                raise FormatError("Reads are improperly paired. Read name '{0}' "
                    "in file 1 does not match '{1}' in file 2.".format(r1.name, r2.name))
            yield (r1, r2)

    def close(self):
        self.reader1.close()
        self.reader2.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

class InterleavedSequenceReader(object):
    """
    Read paired-end reads from an interleaved FASTQ file.
    """
    def __init__(self, file, colorspace=False, fileformat=None):
        self.reader = open_reader(file, colorspace=colorspace, fileformat=fileformat)
        self.delivers_qualities = self.reader.delivers_qualities

    def __iter__(self):
        # Avoid usage of zip() below since it will consume one item too many.
        it = iter(self.reader)
        for r1 in it:
            try:
                r2 = next(it)
            except StopIteration:
                raise FormatError("Interleaved input file incomplete: Last record has no partner.")
            if not sequence_names_match(r1, r2):
                raise FormatError("Reads are improperly paired. Name {0!r} "
                    "(first) does not match {1!r} (second).".format(r1.name, r2.name))
            yield (r1, r2)

    def close(self):
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

def open_reader(file1, file2=None, qualfile=None, colorspace=False, fileformat=None,
                interleaved=False, qualities=None):
    """
    Open sequence files in FASTA or FASTQ format for reading. This is
    a factory that returns an instance of one of the ...Reader
    classes also defined in this module.

    file1, file2, qualfile -- Paths to regular or compressed files or file-like
        objects. Use file1 if data is single-end. If also file2 is provided,
        sequences are paired. If qualfile is given, then file1 must be a FASTA
        file and sequences are single-end. One of file2 and qualfile must always
        be None (no paired-end data is supported when reading qualfiles).

    interleaved -- If True, then file1 contains interleaved paired-end data.
        file2 and qualfile must be None in this case.

    colorspace -- If True, instances of the Colorspace... classes
        are returned.

    fileformat -- If set to None, file format is autodetected from the file name
        extension. Set to 'fasta', 'fastq', or 'sra-fastq' to not auto-detect.
        Colorspace is not auto-detected and must always be requested explicitly.

    qualities -- When mode is 'w' and fileformat is None, this can be set to
        True or False to specify whether the written sequences will have quality
        values. This is is used in two ways:
        * If the output format cannot be determined (unrecognized extension
          etc), no exception is raised, but fasta or fastq format is chosen
          appropriately.
        * When False (no qualities available), an exception is raised when the
          auto-detected output format is FASTQ.

    """
    if interleaved and (file2 is not None or qualfile is not None):
        raise ValueError("When interleaved is set, file2 and qualfile must be None")
    if file2 is not None and qualfile is not None:
        raise ValueError("Setting both file2 and qualfile is not supported")
    
    if file2 is not None:
        return PairedSequenceReader(file1, file2, colorspace, fileformat)

    if interleaved:
        return InterleavedSequenceReader(file1, colorspace, fileformat)
        
    if qualfile is not None:
        if colorspace:
            # read from .(CS)FASTA/.QUAL
            return ColorspaceFastaQualReader(file1, qualfile)
        else:
            return FastaQualReader(file1, qualfile)
    
    fastq_handler = ColorspaceFastqReader if colorspace else FastqReader
    fasta_handler = ColorspaceFastaReader if colorspace else FastaReader

    if fileformat is None and file1 != "-":
        fileformat = guess_format_from_name(file1)
    
    if fileformat is None:
        if file1 == "-":
            file1 = sys.stdin if mode == 'r' else sys.stdout
        for line in file1:
            if line.startswith('#'):
                # Skip comment lines (needed for csfasta)
                continue
            if line.startswith('>'):
                fileformat = 'fasta'
            elif line.startswith('@'):
                fileformat = 'fastq'
            file1 = FileWithPrependedLine(file1, line)
            break
        
    fileformat = fileformat.lower()
    if fileformat == 'fastq' and qualities is False:
        raise ValueError("Output format cannot be FASTQ since no quality values are available.")
    if fileformat == 'fasta':
        return fasta_handler(file1)
    if fileformat == 'fastq':
        return fastq_handler(file1)
    if fileformat == 'sra-fastq' and colorspace:
        return SRAColorspaceFastqReader(file1)
    
    raise UnknownFileType(
        "File format {0!r} is unknown (expected 'sra-fastq' (only for colorspace), "
        "'fasta' or 'fastq').".format(fileformat))

def guess_format_from_name(file, raise_on_failure=False):
    # Detect file format
    name = None
    if isinstance(file, str):
        name = file
    elif hasattr(file, "name"):	 # seems to be an open file-like object
        name = file.name
    
    if name:
        name, ext1, ext2 = _splitext(name)
        ext = ext1.lower()
        if ext in ['.fasta', '.fa', '.fna', '.csfasta', '.csfa']:
            return 'fasta'
        elif ext in ['.fastq', '.fq'] or (ext == '.txt' and name.endswith('_sequence')):
            return 'fastq'
    
    if raise_on_failure:
        raise UnknownFileType("Could not determine whether file {0!r} is FASTA "
            "or FASTQ: file name extension {1!r} not recognized".format(file, ext))

def add_suffix_to_path(path, suffix):
    """
    Add the suffix (str or int) after the file name but
    before the extension.
    """
    name, ext1, ext2 = _splitext(path)
    return "{}{}{}{}".format(name, suffix, ext1, ext2 or "")
    
def _splitext(name):
    ext1 = ext2 = None
    for ext in ('.gz', '.xz', '.bz2'):
        if name.endswith(ext):
            ext2 = ext
            name = name[:-len(ext)]
            break
    name, ext1 = splitext(name)
    return (name, ext1, ext2)

## Converting reads to strings ##

class FastaFormat(object):
    def __init__(self, line_length=None):
        self.text_wrapper = None
        if line_length:
            from textwrap import TextWrapper
            self.text_wrapper = TextWrapper(width=line_length)
    
    def format(self, read):
        return self.format_entry(read.name, read.sequence)
    
    def format_entry(self, name, sequence):
        if self.text_wrapper:
            sequence = self.text_wrapper.fill(sequence)
        return "".join((">", name, "\n", sequence, "\n"))

class ColorspaceFastaFormat(FastaFormat):
    def format(self, read):
        return self.format_entry(read.name, read.primer + read.sequence)

class FastqFormat(object):
    def format(self, read):
        return self.format_entry(read.name, read.sequence, read.qualities, read.name2)
    
    def format_entry(self, name, sequence, qualities, name2=""):
        return "".join((
            '@', name, '\n',
            sequence, '\n+',
            name2, '\n',
            qualities, '\n'
        ))

class ColorspaceFastqFormat(FastqFormat):
    def format(self, read):
        return self.format_entry(read.name, read.primer + read.sequence, read.qualities)

class SingleEndFormatter(object):
    def __init__(self, seq_format, file1):
        self.seq_format = seq_format
        self.file1 = file1
        self.written = 0
        self.read1_bp = 0
        self.read2_bp = 0
    
    def format(self, result, read1, read2=None):
        result[self.file1].append(self.seq_format.format(read1))
        self.written += 1
        self.read1_bp += len(read1)
    
    @property
    def written_bp(self):
        return (self.read1_bp, self.read2_bp)

class InterleavedFormatter(SingleEndFormatter):
    def format(self, result, read1, read2=None):
        result[self.file1].extend((
            self.seq_format.format(read1),
            self.seq_format.format(read2)
        ))
        self.written += 1
        self.read1_bp += len(read1)
        self.read2_bp += len(read2)

class PairedEndFormatter(SingleEndFormatter):
    def __init__(self, seq_format, file1, file2):
        super(PairedEndFormatter, self).__init__(seq_format, file1)
        self.file2 = file2
    
    def format(self, result, read1, read2):
        result[self.file1].append(self.seq_format.format(read1))
        result[self.file2].append(self.seq_format.format(read2))
        self.written += 1
        self.read1_bp += len(read1)
        self.read2_bp += len(read2)

class DelimFormatter(object):
    def __init__(self, path, delim=' '):
        self.path = path
        self.delim = delim
    
    def _format(self, result, fields):
        result[self.path].append("".join((
            self.delim.join(str(f) for f in fields),
            "\n"
        )))

class RestFormatter(DelimFormatter):
    def format(self, result, read):
        if read.match:
            rest = read.match.rest()
            if len(rest) > 0:
                self._format(result, (rest, read.name))

class InfoFormatter(DelimFormatter):
    def __init__(self, path):
        super(InfoFormatter, self).__init__(path, delim='\t')
    
    def format(self, result, read):
        if read.match:
            for m in read.match_info:
                self._format(result, m[0:11])
        else:
            seq = read.sequence
            qualities = read.qualities if read.qualities is not None else ''
            self._format(result, (read.name, -1, seq, qualities))

class WildcardFormatter(DelimFormatter):
    def format(self, result, read):
        if read.match:
            self._format(result, (read.match.wildcards(), read.name))

def create_seq_formatter(file1, file2=None, interleaved=False, **kwargs):
    seq_format = get_format(file1, **kwargs)
    if file2 is not None:
        return PairedEndFormatter(seq_format, file1, file2)
    elif interleaved:
        return InterleavedFormatter(seq_format, file1)
    else:
        return SingleEndFormatter(seq_format, file1)

def get_format(file, fileformat=None, colorspace=False, qualities=None, line_length=None):
    """
    fileformat -- If set to None, file format is autodetected from the file name
        extension. Set to 'fasta', 'fastq', or 'sra-fastq' to not auto-detect.
        Colorspace is not auto-detected and must always be requested explicitly.

    colorspace -- If True, instances of the Colorspace... formats
        are returned.
    
    qualities -- When fileformat is None, this can be set to
        True or False to specify whether the written sequences will have quality
        values. This is is used in two ways:
        * If the output format cannot be determined (unrecognized extension
          etc), no exception is raised, but fasta or fastq format is chosen
          appropriately.
        * When False (no qualities available), an exception is raised when the
          auto-detected output format is FASTQ.
    """
    if fileformat is None:
        fileformat = guess_format_from_name(file, raise_on_failure=qualities is None)
    
    if fileformat is None:
        if qualities is True:
            # Format not recognized, but know we want to write reads with qualities
            fileformat = 'fastq'
        elif qualities is False:
            # Same, but we know that we want to write reads without qualities
            fileformat = 'fasta'
    
    if fileformat is None:
        raise UnknownFileType("Could not determine file type.")
    
    if fileformat == 'fastq' and qualities is False:
        raise ValueError("Output format cannot be FASTQ since no quality values are available.")
    
    fileformat = fileformat.lower()
    if fileformat == 'fasta':
        if colorspace:
            return ColorspaceFastaFormat(line_length)
        else:
            return FastaFormat(line_length)
    elif fileformat == 'fastq':
        if colorspace:
            return ColorspaceFastqFormat()
        else:
            return FastqFormat()
    else:
        raise UnknownFileType("File format {0!r} is unknown (expected 'fasta' or 'fastq').".format(fileformat))

class Formatters(object):
    def __init__(self, output, seq_formatter_args):
        self.output = output
        self.multiplexed = output is not None and '{name}' in output
        self.seq_formatter_args = seq_formatter_args
        self.seq_formatters = {}
        self.mux_formatters = {}
        self.info_formatters = []
        self.discarded = 0
    
    def add_seq_formatter(self, filter_type, file1, file2=None):
        self.seq_formatters[filter_type] = create_seq_formatter(
            file1, file2, **self.seq_formatter_args)
    
    def add_info_formatter(self, formatter):
        self.info_formatters.append(formatter)
    
    def get_mux_formatter(self, name):
        assert self.multiplexed
        if name not in self.mux_formatters:
            path = self.output.format(name=name)
            self.mux_formatters[name] = create_seq_formatter(path, **self.seq_formatter_args)
        return self.mux_formatters[name]
    
    def get_seq_formatters(self):
        return (
            set(f for f in self.seq_formatters.values() if f.written > 0) |
            set(f for f in self.mux_formatters.values() if f.written > 0)
        )
    
    def summary(self):
        seq_formatters = self.get_seq_formatters()
        written = sum(f.written for f in seq_formatters)
        read_bp = [
            sum(f.read1_bp for f in seq_formatters),
            sum(f.read2_bp for f in seq_formatters)
        ]
        return (written, read_bp)
    
    def format(self, result, dest, read1, read2=None):
        if self.multiplexed and (dest == NoFilter) and read1.match:
            name = read1.match.adapter.name
            formatter = self.get_mux_formatter(name)
            formatter.format(result, read1, read2)
        elif dest in self.seq_formatters:
            self.seq_formatters[dest].format(result, read1, read2)
        else:
            self.discarded += 1
        
        for f in self.info_formatters:
            f.format(result, read1)
            if read2:
                f.format(result, read2)

## Writing sequences to files ##

class Writers(object):
    def __init__(self, force_create):
        self.writers = {}
        self.force_create = force_create
        self.suffix = None
    
    def write_result(self, result, compressed=False):
        """
        Write results to file.
        
        result :: dict with keys being file descriptors
        and values being data (either bytes or strings).
        Strings are expected to already have appropriate
        line-endings.
        compressed :: whether data has already been compressed.
        """
        for file_desc, data in result.items():
            self.write(file_desc, data, compressed)
    
    def write(self, file_desc, data, compressed=False):
        """
        Write data to file.
        """
        if compressed:
            path, mode = file_desc
        else:
            path = file_desc
        if path not in self.writers:
            if self.suffix:
                real_path = add_suffix_to_path(path, self.suffix)
            else:
                real_path = path
            # TODO: test whether O_NONBLOCK allows non-blocking write to NFS
            if compressed:
                self.writers[path] = open_output(real_path, mode)
            else:
                self.writers[path] = xopen(real_path, "w")
        self.writers[path].write(data)
    
    def close(self):
        for path in self.force_create:
            if path not in self.writers:
                with open_output(path, "w"):
                    pass
        for writer in self.writers.values():
            writer.close()
