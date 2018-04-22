# coding: utf-8
from pytest import raises
from collections import defaultdict
import random
import os
from io import StringIO
import shutil
from textwrap import dedent
from tempfile import mkdtemp
from unittest import skipIf, TestCase
from atropos.io import xopen, open_output
from atropos.io.seqio import (
    Sequence,
    ColorspaceSequence,
    FormatError,
    FastaReader,
    FastqReader,
    FastaQualReader,
    InterleavedSequenceReader,
    FastaFormat,
    FastqFormat,
    InterleavedFormatter,
    SingleEndFormatter,
    SingleEndSAMFormatter,
    PairedEndSAMFormatter,
    create_seq_formatter,
    open_reader as openseq,
    sequence_names_match,
)
from atropos.util import ALPHABETS
from .utils import temporary_path

# files tests/data/simple.fast{q,a}
simple_fastq = [
    Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
    Sequence("second_sequence", "SEQUENCE2", "83<??:(61"),
]
simple_fasta = [Sequence(x.name, x.sequence, None) for x in simple_fastq]


class TestAlphabet(TestCase):

    def test_alphabet(self):
        alphabet = ALPHABETS['dna']
        for base in ('A', 'C', 'G', 'T', 'N'):
            assert base in alphabet
        assert 'X' not in alphabet
        assert alphabet.resolve('X') == 'N'


class TestSequence(TestCase):

    def test_too_many_qualities(self):
        with raises(FormatError):
            Sequence(name="name", sequence="ACGT", qualities="#####")

    def test_too_many_qualities_colorspace(self):
        with raises(FormatError):
            ColorspaceSequence(name="name", sequence="T0123", qualities="#####")

    def test_invalid_primer(self):
        with raises(FormatError):
            ColorspaceSequence(name="name", sequence="K0123", qualities="####")


class TestFastaReader(TestCase):

    def test(self):
        with self.assertRaises(ValueError):
            with FastaReader(None) as _:
                pass

        with FastaReader("tests/data/simple.fasta") as f:
            reads = list(f)
        assert reads == simple_fasta

        fasta = StringIO(">first_sequence\nSEQUENCE1\n>second_sequence\nSEQUENCE2\n")
        reads = list(FastaReader(fasta))
        assert reads == simple_fasta

    def test_with_comments(self):
        fasta = StringIO(
            dedent("""
            # a comment
            # another one
            >first_sequence
            SEQUENCE1
            >second_sequence
            SEQUENCE2
            """)
        )
        reads = list(FastaReader(fasta))
        assert reads == simple_fasta

    def test_wrong_format(self):
        with raises(FormatError):
            fasta = StringIO(
                dedent("""
                # a comment
                # another one
                unexpected
                >first_sequence
                SEQUENCE1
                >second_sequence
                SEQUENCE2
                """)
            )
            _ = list(FastaReader(fasta))

    def test_fastareader_keeplinebreaks(self):
        with FastaReader("tests/data/simple.fasta", keep_linebreaks=True) as f:
            reads = list(f)
        assert reads[0] == simple_fasta[0]
        assert reads[1].sequence == 'SEQUEN\nCE2'

    def test_context_manager(self):
        filename = "tests/data/simple.fasta"
        with open(filename) as f:
            assert not f.closed
            _ = list(openseq(f))
            assert not f.closed
        assert f.closed
        with FastaReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            _ = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None
        # Open it a second time
        with FastaReader(filename) as _:
            pass


class TestFastqReader(TestCase):
    def test_fastqreader(self):
        with FastqReader("tests/data/simple.fastq") as f:
            reads = list(f)
        assert reads == simple_fastq

    def test_fastqreader_dos(self):
        with FastqReader("tests/data/dos.fastq") as f:
            dos_reads = list(f)
        with FastqReader("tests/data/small.fastq") as f:
            unix_reads = list(f)
        assert dos_reads == unix_reads

    def test_fastq_wrongformat(self):
        with raises(FormatError), FastqReader("tests/data/withplus.fastq") as f:
            _ = list(f)

    def test_fastq_incomplete(self):
        fastq = StringIO("@name\nACGT+\n")
        with raises(FormatError), FastqReader(fastq) as fq:
            list(fq)

    def test_context_manager(self):
        filename = "tests/data/simple.fastq"
        with open(filename) as f:
            assert not f.closed
            _ = list(openseq(f))
            assert not f.closed
        assert f.closed
        with FastqReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            _ = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None

    def test_alphabet(self):
        filename = "tests/data/bad_bases.fq"
        with FastqReader(filename, alphabet=ALPHABETS['dna']) as f:
            reads = list(f)
            assert reads[0].sequence == 'ACGNGGACT'
            assert reads[1].sequence == 'CGGACNNNC'


class TestFastaQualReader(TestCase):
    def test_mismatching_read_names(self):
        with raises(FormatError):
            fasta = StringIO(">name\nACG")
            qual = StringIO(">nome\n3 5 7")
            list(FastaQualReader(fasta, qual))

    def test_invalid_quality_value(self):
        with raises(FormatError):
            fasta = StringIO(">name\nACG")
            qual = StringIO(">name\n3 xx 7")
            list(FastaQualReader(fasta, qual))


class TestSeqioOpen(TestCase):
    _tmpdir = None

    def setUp(self):
        self._tmpdir = mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tmpdir)

    def test_sequence_reader(self):
        # test the autodetection
        with openseq("tests/data/simple.fastq") as f:
            reads = list(f)
        assert reads == simple_fastq
        with openseq("tests/data/simple.fasta") as f:
            reads = list(f)
        assert reads == simple_fasta
        with open("tests/data/simple.fastq") as f:
            reads = list(openseq(f))
        assert reads == simple_fastq
        # make the name attribute unavailable
        f = StringIO(open("tests/data/simple.fastq").read())
        reads = list(openseq(f))
        assert reads == simple_fastq
        f = StringIO(open("tests/data/simple.fasta").read())
        reads = list(openseq(f))
        assert reads == simple_fasta

    def test_autodetect_fasta_format(self):
        path = os.path.join(self._tmpdir, 'tmp.fasta')
        fmt = create_seq_formatter(path)
        assert isinstance(fmt, SingleEndFormatter)
        assert isinstance(fmt.seq_format, FastaFormat)
        write_seq_output(simple_fasta, fmt)
        assert list(openseq(path)) == simple_fasta

    def test_write_qualities_to_fasta(self):
        path = os.path.join(self._tmpdir, 'tmp.fasta')
        fmt = create_seq_formatter(path, qualities=True)
        assert isinstance(fmt, SingleEndFormatter)
        assert isinstance(fmt.seq_format, FastaFormat)
        write_seq_output(simple_fasta, fmt)
        assert list(openseq(path)) == simple_fasta

    def test_autodetect_fastq_format(self):
        path = os.path.join(self._tmpdir, 'tmp.fastq')
        fmt = create_seq_formatter(path)
        assert isinstance(fmt, SingleEndFormatter)
        assert isinstance(fmt.seq_format, FastqFormat)
        write_seq_output(simple_fastq, fmt)
        assert list(openseq(path)) == simple_fastq

    def test_fastq_qualities_missing(self):
        with raises(ValueError):
            path = os.path.join(self._tmpdir, 'tmp.fastq')
            create_seq_formatter(path, qualities=False)


def write_seq_output(reads, fmt):
    result = defaultdict(list)
    for read in reads:
        fmt.format(result, read)
    for path, seqs in result.items():
        with xopen(path, "w") as f:
            f.write("".join(seqs))


class TestInterleavedReader(TestCase):
    def test(self):
        expected = [
            (
                Sequence('read1/1 some text', 'TTATTTGTCTCCAGC', '##HHHHHHHHHHHHH'),
                Sequence('read1/2 other text', 'GCTGGAGACAAATAA', 'HHHHHHHHHHHHHHH'),
            ),
            (
                Sequence('read3/1', 'CCAACTTGATATTAATAACA', 'HHHHHHHHHHHHHHHHHHHH'),
                Sequence('read3/2', 'TGTTATTAATATCAAGTTGG', '#HHHHHHHHHHHHHHHHHHH'),
            ),
        ]
        reads = list(InterleavedSequenceReader("tests/cut/interleaved.fastq"))
        for (r1, r2), (e1, e2) in zip(reads, expected):
            print(r1, r2, e1, e2)
        assert reads == expected
        with openseq("tests/cut/interleaved.fastq", interleaved=True) as f:
            reads = list(f)
        assert reads == expected

    def test_missing_partner(self):
        with raises(FormatError):
            s = StringIO('@r1\nACG\n+\nHHH')
            list(InterleavedSequenceReader(s))

    def test_incorrectly_paired(self):
        with raises(FormatError):
            s = StringIO('@r1/1\nACG\n+\nHHH\n@wrong_name\nTTT\n+\nHHH')
            list(InterleavedSequenceReader(s))


class TestFastaWriter(TestCase):
    _tmpdir = None
    path = None

    def setUp(self):
        self._tmpdir = mkdtemp()
        self.path = os.path.join(self._tmpdir, 'tmp.fasta')

    def tearDown(self):
        shutil.rmtree(self._tmpdir)

    def test(self):
        fmt = FastaFormat()
        with open_output(self.path, "w") as fw:
            fw.write(fmt.format_entry("name", "CCATA"))
            fw.write(fmt.format_entry("name2", "HELLO"))
        with open(self.path) as t:
            assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'

    def test_linelength(self):
        fmt = FastaFormat(line_length=3)
        with open_output(self.path, "w") as fw:
            fw.write(fmt.format_entry("r1", "ACG"))
            fw.write(fmt.format_entry("r2", "CCAT"))
            fw.write(fmt.format_entry("r3", "TACCAG"))
        with open(self.path) as t:
            x = t.read()
            print(x)
            assert x == '>r1\nACG\n>r2\nCCA\nT\n>r3\nTAC\nCAG\n'

    def test_write_sequence_object(self):
        fmt = FastaFormat()
        with open_output(self.path, "w") as fw:
            fw.write(fmt.format(Sequence("name", "CCATA")))
            fw.write(fmt.format(Sequence("name2", "HELLO")))
        with open(self.path) as t:
            assert t.read() == '>name\nCCATA\n>name2\nHELLO\n'

    def test_write_zero_length_sequence(self):
        fmt = FastaFormat()
        s = fmt.format_entry("name", "")
        assert s == '>name\n\n', '{0!r}'.format(s)


class TestFastqWriter(TestCase):
    _tmpdir = None
    path = None

    def setUp(self):
        self._tmpdir = mkdtemp()
        self.path = os.path.join(self._tmpdir, 'tmp.fastq')

    def tearDown(self):
        shutil.rmtree(self._tmpdir)

    def test(self):
        fmt = FastqFormat()
        with open_output(self.path, "w") as fw:
            fw.write(fmt.format_entry("name", "CCATA", "!#!#!"))
            fw.write(fmt.format_entry("name2", "HELLO", "&&&!&&"))
        with open(self.path) as t:
            assert t.read() == '@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n'

    def test_twoheaders(self):
        fmt = FastqFormat()
        with open_output(self.path, "w") as fw:
            fw.write(fmt.format(Sequence("name", "CCATA", "!#!#!", name2="name")))
            fw.write(fmt.format(Sequence("name2", "HELLO", "&&&!&", name2="name2")))
        with open(self.path) as t:
            assert t.read(
            ) == '@name\nCCATA\n+name\n!#!#!\n@name2\nHELLO\n+name2\n&&&!&\n'


class TestInterleavedWriter(TestCase):
    def test(self):
        reads = [
            (
                Sequence('A/1 comment', 'TTA', '##H'),
                Sequence('A/2 comment', 'GCT', 'HH#'),
            ),
            (Sequence('B/1', 'CC', 'HH'), Sequence('B/2', 'TG', '#H')),
        ]
        fmt = InterleavedFormatter("foo", FastqFormat())
        result = defaultdict(lambda: [])
        for read1, read2 in reads:
            fmt.format(result, read1, read2)
        assert fmt.written == 2
        assert fmt.read1_bp == 5
        assert fmt.read2_bp == 5
        assert "foo" in result
        assert "".join(
            result["foo"]
        ) == (
            '@A/1 comment\nTTA\n+\n##H\n@A/2 comment\nGCT\n+\nHH#\n@B/1\nCC\n+\nHH\n'
            '@B/2\nTG\n+\n#H\n'
        )


class TestSAMWriter(TestCase):
    def test_single_end(self):
        reads = [Sequence('A/1', 'TTA', '##H'), Sequence('B/1', 'CC', 'HH')]
        fmt = SingleEndSAMFormatter("foo")
        result = defaultdict(lambda: [])
        for read in reads:
            fmt.format(result, read)
        assert fmt.written == 2
        assert fmt.read1_bp == 5
        assert fmt.read2_bp == 0
        assert "foo" in result
        result_str = "".join(result["foo"])
        expected = (
            '@HD\tVN:1.5\tSO:unsorted\nA/1\t0\t*\t0\t0\t*\t*\t0\t0\tTTA\t##H\nB/1\t0'
            '\t*\t0\t0\t*\t*\t0\t0\tCC\tHH\n'
        )
        print(result_str)
        print(expected)
        assert result_str == expected

    def test_paired_end(self):
        reads = [
            (Sequence('A/1', 'TTA', '##H'), Sequence('A/2', 'GCT', 'HH#')),
            (Sequence('B/1', 'CC', 'HH'), Sequence('B/2', 'TG', '#H')),
        ]
        fmt = PairedEndSAMFormatter("foo")
        result = defaultdict(lambda: [])
        for read1, read2 in reads:
            fmt.format(result, read1, read2)
        assert fmt.written == 2
        assert fmt.read1_bp == 5
        assert fmt.read2_bp == 5
        assert "foo" in result
        result_str = "".join(result["foo"])
        expected = (
            '@HD\tVN:1.5\tSO:unsorted\nA/1\t65\t*\t0\t0\t*\t*\t0\t0\tTTA\t##H\nA/2\t129'
            '\t*\t0\t0\t*\t*\t0\t0\tGCT\tHH#\nB/1\t65\t*\t0\t0\t*\t*\t0\t0\tCC\tHH\nB/2'
            '\t129\t*\t0\t0\t*\t*\t0\t0\tTG\t#H\n'
        )
        assert result_str == expected


class TestPairedSequenceReader(TestCase):
    def test_sequence_names_match(self):

        def match(name1, name2):
            seq1 = Sequence(name1, 'ACGT')
            seq2 = Sequence(name2, 'AACC')
            return sequence_names_match(seq1, seq2)

        assert match('abc', 'abc')
        assert match('abc/1', 'abc/2')
        assert match('abc.1', 'abc.2')
        assert match('abc1', 'abc2')
        assert not match('abc', 'xyz')


try:
    import ngstream
    ngs_available = True
except ModuleNotFoundError:
    ngs_available = False


@skipIf(not ngs_available, "ngstream library not available")
class TestSraReader(TestCase):
    def test_sra_reader(self):
        with ngstream.sra.SraReader(item_limit=10) as _:
            # atropos_reader = sra_reader(reader)
            pass


def create_truncated_file(path):
    # Random text
    text = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZ') for _ in range(200))
    f = xopen(path, 'w')
    f.write(text)
    f.close()
    f = open(path, 'a')
    f.truncate(os.stat(path).st_size - 10)
    f.close()


def test_truncated_gz():
    with raises(EOFError), temporary_path('truncated.gz') as path:
        create_truncated_file(path)
        f = xopen(path, 'r')
        f.read()
        f.close()


def test_truncated_gz_iter():
    with raises(EOFError), temporary_path('truncated.gz') as path:
        create_truncated_file(path)
        f = xopen(path, 'r', use_system=False)  # work around bug in py3.4
        for _ in f:
            pass
        f.close()
