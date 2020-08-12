from collections import defaultdict
from io import StringIO
import os
import random
from textwrap import dedent

import pytest
from xphyle import open_, xopen

from atropos import __version__
from atropos.errors import FormatError
from atropos.io import InputRead
from atropos.io.sequence import Sequence, ColorspaceSequence
from atropos.io.formatters import (
    FastaFormat,
    FastqFormat,
    InterleavedFormatter,
    SingleEndFormatter,
    SingleEndSAMFormatter,
    PairedEndSAMFormatter,
    create_seq_formatter,
)
from atropos.io.readers import (
    FastaReader,
    FastqReader,
    FastaQualReader,
    InterleavedSequenceReader,
    open_reader,
    sequence_names_match
)
from atropos.utils import no_import
from atropos.utils.ngs import ALPHABETS

from .utils import no_internet

# files tests/data/simple.fast{q,a}
simple_fastq = [
    Sequence("first_sequence", "SEQUENCE1", ":6;;8<=:<"),
    Sequence("second_sequence", "SEQUENCE2", "83<??:(61"),
]
simple_fasta = [Sequence(x.name, x.sequence, None) for x in simple_fastq]


class TestAlphabet:
    def test_alphabet(self):
        alphabet = ALPHABETS["dna"]
        for base in ("A", "C", "G", "T", "N"):
            assert base in alphabet
        assert "X" not in alphabet
        assert alphabet.resolve("X") == "N"


class TestSequence:
    def test_too_many_qualities(self):
        with pytest.raises(FormatError):
            Sequence(name="name", sequence="ACGT", qualities="#####")

    def test_too_many_qualities_colorspace(self):
        with pytest.raises(FormatError):
            ColorspaceSequence(name="name", sequence="T0123", qualities="#####")

    def test_invalid_primer(self):
        with pytest.raises(FormatError):
            ColorspaceSequence(name="name", sequence="K0123", qualities="####")


class TestFastaReader:
    def test(self, input_data):
        with pytest.raises(ValueError):
            with FastaReader(None) as _:
                pass

        with FastaReader(input_data("simple.fasta")) as f:
            reads = list(f)
        assert reads == simple_fasta

        fasta = StringIO(">first_sequence\nSEQUENCE1\n>second_sequence\nSEQUENCE2\n")
        reads = list(FastaReader(fasta))
        assert reads == simple_fasta

    def test_with_comments(self):
        fasta = StringIO(
            dedent(
                """
            # a comment
            # another one
            >first_sequence
            SEQUENCE1
            >second_sequence
            SEQUENCE2
            """
            )
        )
        reads = list(FastaReader(fasta))
        assert reads == simple_fasta

    def test_wrong_format(self):
        with pytest.raises(FormatError):
            fasta = StringIO(
                dedent(
                    """
                # a comment
                # another one
                unexpected
                >first_sequence
                SEQUENCE1
                >second_sequence
                SEQUENCE2
                """
                )
            )
            _ = list(FastaReader(fasta))

    def test_fastareader_keeplinebreaks(self, input_data):
        with FastaReader(input_data("simple.fasta"), keep_linebreaks=True) as f:
            reads = list(f)
        assert reads[0] == simple_fasta[0]
        assert reads[1].sequence == "SEQUEN\nCE2"

    def test_context_manager(self, input_data):
        filename = input_data("simple.fasta")
        with open(filename) as f:
            assert not f.closed
            _ = list(open_reader(f))
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


class TestFastqReader:
    def test_fastqreader(self, input_data):
        with FastqReader(input_data("simple.fastq")) as f:
            reads = list(f)
        assert reads == simple_fastq

    def test_fastqreader_dos(self, input_data):
        with FastqReader(input_data("dos.fastq")) as f:
            dos_reads = list(f)
        with FastqReader(input_data("small.fastq")) as f:
            unix_reads = list(f)
        assert dos_reads == unix_reads

    def test_fastq_wrongformat(self, input_data):
        with pytest.raises(FormatError), FastqReader(input_data("withplus.fastq")) as f:
            _ = list(f)

    def test_fastq_incomplete(self):
        fastq = StringIO("@name\nACGT+\n")
        with pytest.raises(FormatError), FastqReader(fastq) as fq:
            list(fq)

    def test_context_manager(self, input_data):
        filename = input_data("simple.fastq")
        with open(filename) as f:
            assert not f.closed
            _ = list(open_reader(f))
            assert not f.closed
        assert f.closed
        with FastqReader(filename) as sr:
            tmp_sr = sr
            assert not sr._file.closed
            _ = list(sr)
            assert not sr._file.closed
        assert tmp_sr._file is None

    def test_alphabet(self, input_data):
        filename = input_data("bad_bases.fq")
        with FastqReader(filename, alphabet=ALPHABETS["dna"]) as f:
            reads = list(f)
            assert reads[0].sequence == "ACGNGGACT"
            assert reads[1].sequence == "CGGACNNNC"


class TestFastaQualReader:
    def test_mismatching_read_names(self):
        with pytest.raises(FormatError):
            fasta = StringIO(">name\nACG")
            qual = StringIO(">nome\n3 5 7")
            list(FastaQualReader(fasta, qual))

    def test_invalid_quality_value(self):
        with pytest.raises(FormatError):
            fasta = StringIO(">name\nACG")
            qual = StringIO(">name\n3 xx 7")
            list(FastaQualReader(fasta, qual))


class TestSeqioOpen:
    def test_sequence_reader(self, input_data, tmp_path):
        # test the autodetection
        with open_reader(input_data("simple.fastq")) as f:
            reads = list(f)
        assert reads == simple_fastq
        with open_reader(input_data("simple.fasta")) as f:
            reads = list(f)
        assert reads == simple_fasta
        with open(input_data("simple.fastq")) as f:
            reads = list(open_reader(f))
        assert reads == simple_fastq
        # make the name attribute unavailable
        with open(input_data("simple.fastq")) as inp:
            f = StringIO(inp.read())
        reads = list(open_reader(f))
        assert reads == simple_fastq
        with open(input_data("simple.fasta")) as inp:
            f = StringIO(inp.read())
        reads = list(open_reader(f))
        assert reads == simple_fasta

    def test_autodetect_fasta_format(self, tmp_path):
        path = tmp_path / "tmp.fasta"
        fmt = create_seq_formatter(path)
        assert isinstance(fmt, SingleEndFormatter)
        assert isinstance(fmt._seq_format, FastaFormat)
        write_seq_output(simple_fasta, fmt)
        assert list(open_reader(path)) == simple_fasta

    def test_write_qualities_to_fasta(self, tmp_path):
        path = tmp_path / "tmp.fasta"
        fmt = create_seq_formatter(path, qualities=True)
        assert isinstance(fmt, SingleEndFormatter)
        assert isinstance(fmt._seq_format, FastaFormat)
        write_seq_output(simple_fasta, fmt)
        assert list(open_reader(path)) == simple_fasta

    def test_autodetect_fastq_format(self, tmp_path):
        path = tmp_path / "tmp.fastq"
        fmt = create_seq_formatter(path)
        assert isinstance(fmt, SingleEndFormatter)
        assert isinstance(fmt._seq_format, FastqFormat)
        write_seq_output(simple_fastq, fmt)
        assert list(open_reader(path)) == simple_fastq

    def test_fastq_qualities_missing(self, tmp_path):
        with pytest.raises(ValueError):
            path = tmp_path / "tmp.fastq"
            create_seq_formatter(path, qualities=False)


def write_seq_output(reads, fmt):
    result = defaultdict(list)
    for read in reads:
        fmt.format(result, read)
    for path, seqs in result.items():
        with xopen(path, "w") as f:
            f.write("".join(seqs))


class TestInterleavedReader:
    def test(self, expected_data):
        expected = [
            (
                Sequence("read1/1 some text", "TTATTTGTCTCCAGC", "##HHHHHHHHHHHHH"),
                Sequence("read1/2 other text", "GCTGGAGACAAATAA", "HHHHHHHHHHHHHHH"),
            ),
            (
                Sequence("read3/1", "CCAACTTGATATTAATAACA", "HHHHHHHHHHHHHHHHHHHH"),
                Sequence("read3/2", "TGTTATTAATATCAAGTTGG", "#HHHHHHHHHHHHHHHHHHH"),
            ),
        ]
        reads = list(InterleavedSequenceReader(expected_data("interleaved.fastq")))
        assert reads == expected

        with open_reader(
            expected_data("interleaved.fastq"),
            interleaved=True,
            input_read=InputRead.PAIRED
        ) as f:
            reads = list(f)
        assert reads == expected

    def test_missing_partner(self):
        with pytest.raises(FormatError):
            s = StringIO("@r1\nACG\n+\nHHH")
            list(InterleavedSequenceReader(s))

    def test_incorrectly_paired(self):
        with pytest.raises(FormatError):
            s = StringIO("@r1/1\nACG\n+\nHHH\n@wrong_name\nTTT\n+\nHHH")
            list(InterleavedSequenceReader(s))


class TestFastaWriter:
    def test(self, tmp_path):
        path = tmp_path / "tmp"
        fmt = FastaFormat()
        with open_(path, "w") as fw:
            fw.write(fmt.format_entry("name", "CCATA"))
            fw.write(fmt.format_entry("name2", "HELLO"))
        with open_(path) as t:
            assert t.read() == ">name\nCCATA\n>name2\nHELLO\n"

    def test_linelength(self, tmp_path):
        path = tmp_path / "tmp"
        fmt = FastaFormat(line_length=3)
        with open_(path, "w") as fw:
            fw.write(fmt.format_entry("r1", "ACG"))
            fw.write(fmt.format_entry("r2", "CCAT"))
            fw.write(fmt.format_entry("r3", "TACCAG"))
        with open_(path) as t:
            x = t.read()
            print(x)
            assert x == ">r1\nACG\n>r2\nCCA\nT\n>r3\nTAC\nCAG\n"

    def test_write_sequence_object(self, tmp_path):
        path = tmp_path / "tmp"
        fmt = FastaFormat()
        with open_(path, "w") as fw:
            fw.write(fmt.format(Sequence("name", "CCATA")))
            fw.write(fmt.format(Sequence("name2", "HELLO")))
        with open_(path) as t:
            assert t.read() == ">name\nCCATA\n>name2\nHELLO\n"

    def test_write_zero_length_sequence(self):
        fmt = FastaFormat()
        s = fmt.format_entry("name", "")
        assert s == ">name\n\n", "{0!r}".format(s)


class TestFastqWriter:
    def test(self, tmp_path):
        path = tmp_path / "tmp.fastq"
        fmt = FastqFormat()
        with open_(path, "w") as fw:
            fw.write(fmt.format_entry("name", "CCATA", "!#!#!"))
            fw.write(fmt.format_entry("name2", "HELLO", "&&&!&&"))
        with open_(path) as t:
            assert t.read() == "@name\nCCATA\n+\n!#!#!\n@name2\nHELLO\n+\n&&&!&&\n"

    def test_twoheaders(self, tmp_path):
        path = tmp_path / "tmp.fastq"
        fmt = FastqFormat()
        with open_(path, "w") as fw:
            fw.write(fmt.format(Sequence("name", "CCATA", "!#!#!", name2="name")))
            fw.write(fmt.format(Sequence("name2", "HELLO", "&&&!&", name2="name2")))
        with open_(path) as t:
            assert (
                t.read() == "@name\nCCATA\n+name\n!#!#!\n@name2\nHELLO\n+name2\n&&&!&\n"
            )


class TestInterleavedWriter:
    def test(self):
        reads = [
            (
                Sequence("A/1 comment", "TTA", "##H"),
                Sequence("A/2 comment", "GCT", "HH#"),
            ),
            (Sequence("B/1", "CC", "HH"), Sequence("B/2", "TG", "#H")),
        ]
        fmt = InterleavedFormatter("foo", FastqFormat())
        result = defaultdict(lambda: [])
        for read1, read2 in reads:
            fmt.format(result, read1, read2)
        assert fmt.written == 2
        assert fmt.read1_bp == 5
        assert fmt.read2_bp == 5
        assert "foo" in result
        assert "".join(result["foo"]) == (
            "@A/1 comment\nTTA\n+\n##H\n@A/2 comment\nGCT\n+\nHH#\n@B/1\nCC\n+\nHH\n"
            "@B/2\nTG\n+\n#H\n"
        )


class TestSAMWriter:
    def test_single_end(self):
        reads = [Sequence("A/1", "TTA", "##H"), Sequence("B/1", "CC", "HH")]
        expected_header = "@HD\tVN:1.6\tSO:unsorted\n"
        expected_pg = f"@PG\tID:Atropos\tPN:Atropos\tVN:{__version__}\tCL:test\n"
        expected_reads = (
            "A/1\t4\t*\t0\t0\t*\t*\t0\t0\tTTA\t##H\n"
            "B/1\t4\t*\t0\t0\t*\t*\t0\t0\tCC\tHH\n"
        )
        fmt = SingleEndSAMFormatter("foo", command="test")

        # simulate writing multiple batches - should only have one header
        # batch 1
        batch1 = defaultdict(lambda: [])
        for read in reads:
            fmt.format(batch1, read)
        assert fmt.written == 2
        assert fmt.read1_bp == 5
        assert fmt.read2_bp == 0
        assert "foo" in batch1
        result_str1 = "".join(batch1["foo"])
        assert result_str1 == expected_header + expected_pg + expected_reads

        # batch 2
        batch2 = defaultdict(lambda: [])
        for read in reads:
            fmt.format(batch2, read)
        assert fmt.written == 4
        assert fmt.read1_bp == 10
        assert fmt.read2_bp == 0
        assert "foo" in batch2
        result_str2 = "".join(batch2["foo"])
        assert result_str2 == expected_reads

    def test_paired_end(self):
        reads = [
            (Sequence("A/1", "TTA", "##H"), Sequence("A/2", "GCT", "HH#")),
            (Sequence("B/1", "CC", "HH"), Sequence("B/2", "TG", "#H")),
        ]
        expected_header = "@HD\tVN:1.6\tSO:unsorted\n"
        expected_pg = f"@PG\tID:Atropos\tPN:Atropos\tVN:{__version__}\tCL:test\n"
        expected_reads = (
            "A/1\t77\t*\t0\t0\t*\t*\t0\t0\tTTA\t##H\n"
            "A/2\t141\t*\t0\t0\t*\t*\t0\t0\tGCT\tHH#\n"
            "B/1\t77\t*\t0\t0\t*\t*\t0\t0\tCC\tHH\n"
            "B/2\t141\t*\t0\t0\t*\t*\t0\t0\tTG\t#H\n"
        )
        fmt = PairedEndSAMFormatter("foo", command="test")

        # simulate writing multiple batches - should only have one header
        # batch 1
        batch1 = defaultdict(lambda: [])
        for read1, read2 in reads:
            fmt.format(batch1, read1, read2)
        assert fmt.written == 2
        assert fmt.read1_bp == 5
        assert fmt.read2_bp == 5
        assert "foo" in batch1
        result_str1 = "".join(batch1["foo"])
        assert result_str1 == expected_header + expected_pg + expected_reads

        # batch2
        batch2 = defaultdict(lambda: [])
        for read1, read2 in reads:
            fmt.format(batch2, read1, read2)
        assert fmt.written == 4
        assert fmt.read1_bp == 10
        assert fmt.read2_bp == 10
        assert "foo" in batch2
        result_str2 = "".join(batch2["foo"])
        assert result_str2 == expected_reads

    def test_tags(self):
        reads = [
            (
                Sequence("A/1", "TTA", "##H", annotations={"XX": "XX:i:0"}),
                Sequence("A/2", "GCT", "HH#", annotations={"XX": "XX:i:1"})
            )
        ]
        expected_header = "@HD\tVN:1.6\tSO:unsorted\n"
        expected_pg = f"@PG\tID:Atropos\tPN:Atropos\tVN:{__version__}\tCL:test\n"
        expected_reads = (
            "A/1\t77\t*\t0\t0\t*\t*\t0\t0\tTTA\t##H\tXX:i:0\n"
            "A/2\t141\t*\t0\t0\t*\t*\t0\t0\tGCT\tHH#\tXX:i:1\n"
        )
        fmt = PairedEndSAMFormatter("foo", command="test")
        # simulate writing multiple batches - should only have one header
        # batch 1
        batch1 = defaultdict(lambda: [])
        for read1, read2 in reads:
            fmt.format(batch1, read1, read2)
        assert fmt.written == 1
        assert fmt.read1_bp == 3
        assert fmt.read2_bp == 3
        assert "foo" in batch1
        result_str1 = "".join(batch1["foo"])
        assert result_str1 == expected_header + expected_pg + expected_reads


class TestPairedSequenceReader:
    def test_sequence_names_match(self):
        def match(name1, name2):
            seq1 = Sequence(name1, "ACGT")
            seq2 = Sequence(name2, "AACC")
            return sequence_names_match(seq1, seq2)

        assert match("abc", "abc")
        assert match("abc/1", "abc/2")
        assert match("abc.1", "abc.2")
        assert match("abc1", "abc2")
        assert not match("abc", "xyz")


SRA_ACCESSION = "ERR2009169"
SRA_SEQ1 = Sequence(
    "D00442:178:C87AVANXX:2:2203:1465:2211",
    "CAGCTTCTTCATCATGTCCTCTACTTTCTTGGCCCGCTCGGCAGGCCCCAAC",
    "CCCCCFGGGGGGGGGGGGGGGFGGGGGGGGGCGGGGGGGG@BEDDGGGGGGG"
)


@pytest.mark.skipif(
    no_internet("https://ncbi.nlm.nih.gov") or
    no_import("ngstream") or
    no_import("ngs"),
    reason="ngstream library not available"
)
class TestSraReader:
    def test_sra_reader(self):
        import ngstream
        with ngstream.open(SRA_ACCESSION, "sra") as stream:
            reader = open_reader(ngstream_reader=stream)
            assert next(iter(reader)) == SRA_SEQ1


def create_truncated_file(path):
    # Random text
    text = "".join(random.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZ") for _ in range(200))
    f = xopen(path, "w")
    f.write(text)
    f.close()
    f = open(path, "a")
    f.truncate(os.stat(path).st_size - 10)
    f.close()


def test_truncated_gz(tmp_path):
    with pytest.raises(EOFError):
        path = tmp_path / "truncated.gz"
        create_truncated_file(path)
        f = xopen(path, "r")
        f.read()
        f.close()


def test_truncated_gz_iter(tmp_path):
    with pytest.raises(EOFError):
        path = tmp_path / "truncated.gz"
        create_truncated_file(path)
        f = xopen(path, "r", use_system=False)  # work around bug in py3.4
        for _ in f:
            pass
        f.close()
