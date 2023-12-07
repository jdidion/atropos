# coding: utf-8
# TODO
# test with the --output option
# test reading from standard input
from io import StringIO
import os
from pytest import raises
import sys
from atropos.commands import execute_cli, get_command
from unittest import skipIf
from .utils import (
    run,
    files_equal,
    datapath,
    cutpath,
    redirect_stderr,
    no_internet,
    no_import,
)


def test_example(tmp_path):
    run(tmp_path, "-N -b ADAPTER", "example.fa", "example.fa")


def test_small(tmp_path):
    run(tmp_path, "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq")


def test_empty(tmp_path):
    """empty input"""
    run(tmp_path, "-a TTAGACATATCTCCGTCG", "empty.fastq", "empty.fastq")


def test_newlines(tmp_path):
    """DOS/Windows newlines"""
    run(tmp_path, "-e 0.12 -b TTAGACATATCTCCGTCG", "dos.fastq", "dos.fastq")


def test_lowercase(tmp_path):
    """lowercase adapter"""
    run(tmp_path, "-b ttagacatatctccgtcg", "lowercase.fastq", "small.fastq")


def test_rest(tmp_path):
    """-r/--rest-file"""
    path = tmp_path / "rest.tmp"
    run(tmp_path, ["-b", "ADAPTER", "-N", "-r", str(path)], "rest.fa", "rest.fa")
    assert files_equal(datapath("rest.txt"), str(path))


def test_restfront(tmp_path):
    path = tmp_path / "rest.tmp"
    run(tmp_path, ["-g", "ADAPTER", "-N", "-r", str(path)], "restfront.fa", "rest.fa")
    assert files_equal(datapath("restfront.txt"), str(path))


def test_discard(tmp_path):
    """--discard"""
    run(tmp_path, "-b TTAGACATATCTCCGTCG --discard", "discard.fastq", "small.fastq")


def test_discard_untrimmed(tmp_path):
    """--discard-untrimmed"""
    run(
        tmp_path,
        "-b CAAGAT --discard-untrimmed",
        "discard-untrimmed.fastq",
        "small.fastq",
    )


def test_plus(tmp_path):
    """test if sequence name after the "+" is retained"""
    run(tmp_path, "-e 0.12 -b TTAGACATATCTCCGTCG", "plus.fastq", "plus.fastq")


def test_extensiontxtgz(tmp_path):
    """automatic recognition of "_sequence.txt.gz" extension"""
    run(tmp_path, "-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz")


def test_format(tmp_path):
    """the -f/--format parameter"""
    run(
        tmp_path,
        "-f fastq -b TTAGACATATCTCCGTCG",
        "small.fastq",
        "small.myownextension",
    )


def test_minimum_length(tmp_path):
    """-m/--minimum-length"""
    run(tmp_path, "-c -m 5 -a 330201030313112312", "minlen.fa", "lengths.fa")


def test_too_short(tmp_path):
    """--too-short-output"""
    run(
        tmp_path,
        "-c -m 5 -a 330201030313112312 --too-short-output tooshort.tmp.fa",
        "minlen.fa",
        "lengths.fa",
    )
    assert files_equal(datapath("tooshort.fa"), "tooshort.tmp.fa")
    os.remove("tooshort.tmp.fa")


def test_too_short_no_primer(tmp_path):
    """--too-short-output and --trim-primer"""
    run(
        tmp_path,
        "-c -m 5 -a 330201030313112312 --trim-primer --too-short-output tooshort.tmp.fa",
        "minlen.noprimer.fa",
        "lengths.fa",
    )
    assert files_equal(datapath("tooshort.noprimer.fa"), "tooshort.tmp.fa")
    os.remove("tooshort.tmp.fa")


def test_maximum_length(tmp_path):
    """-M/--maximum-length"""
    run(tmp_path, "-c -M 5 -a 330201030313112312", "maxlen.fa", "lengths.fa")


def test_too_long(tmp_path):
    """--too-long-output"""
    run(
        tmp_path,
        "-c -M 5 --too-long-output toolong.tmp.fa -a 330201030313112312",
        "maxlen.fa",
        "lengths.fa",
    )
    assert files_equal(datapath("toolong.fa"), "toolong.tmp.fa")
    os.remove("toolong.tmp.fa")


def test_length_tag(tmp_path):
    """454 data; -n and --length-tag"""
    run(
        tmp_path,
        "-n 3 -e 0.1 --length-tag length= "
        "-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "
        "-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
        "454.fa",
        "454.fa",
    )


def test_overlap_a(tmp_path):
    """-O/--overlap with -a (-c omitted on purpose)"""
    run(tmp_path, "-O 10 -a 330201030313112312 -e 0.0 -N", "overlapa.fa", "overlapa.fa")


def test_overlap_b(tmp_path):
    """-O/--overlap with -b"""
    run(tmp_path, "-O 10 -b TTAGACATATCTCCGTCG -N", "overlapb.fa", "overlapb.fa")


def test_qualtrim(tmp_path):
    """-q with low qualities"""
    run(tmp_path, "-q 10 -a XXXXXX", "lowqual.fastq", "lowqual.fastq")


def test_qualbase(tmp_path):
    """-q with low qualities, using ascii(quality+64) encoding"""
    run(
        tmp_path,
        "-q 10 --quality-base 64 -a XXXXXX",
        "illumina64.fastq",
        "illumina64.fastq",
    )


def test_quality_trim_only(tmp_path):
    """only trim qualities, do not remove adapters"""
    run(tmp_path, "-q 10 --quality-base 64", "illumina64.fastq", "illumina64.fastq")


def test_twoadapters(tmp_path):
    """two adapters"""
    run(
        tmp_path,
        "-a AATTTCAGGAATT -a GTTCTCTAGTTCT",
        "twoadapters.fasta",
        "twoadapters.fasta",
    )


def test_polya(tmp_path):
    """poly-A tails"""
    run(
        tmp_path,
        "-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "polya.fasta",
        "polya.fasta",
    )


def test_polya_brace_notation(tmp_path):
    """poly-A tails"""
    run(tmp_path, "-m 24 -O 10 -a A{35}", "polya.fasta", "polya.fasta")


def test_mask_adapter(tmp_path):
    """mask adapter with N (reads maintain the same length)"""
    run(
        tmp_path,
        "-b CAAG -n 3 --mask-adapter",
        "anywhere_repeat.fastq",
        "anywhere_repeat.fastq",
    )


def test_gz_multiblock(tmp_path):
    """compressed gz file with multiple blocks (created by concatenating two .gz files)"""
    run(tmp_path, "-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz")


def test_suffix(tmp_path):
    """-y/--suffix parameter, combined with _F3"""
    run(
        tmp_path,
        "-c -e 0.12 -a 1=330201030313112312 -y _my_suffix_{name} --strip-f3",
        "suffix.fastq",
        "solid.csfasta",
        qualfile="solid.qual",
    )


def test_read_wildcard(tmp_path):
    """test wildcards in reads"""
    run(tmp_path, "--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa")


def test_adapter_wildcard(tmp_path):
    """wildcards in adapter"""
    for adapter_type, expected in (
        ("-a", "wildcard_adapter.fa"),
        ("-b", "wildcard_adapter_anywhere.fa"),
    ):
        path = tmp_path / "wildcardtmp.txt"
        run(
            tmp_path,
            "--wildcard-file {0} {1} ACGTNNNACGT".format(path, adapter_type),
            expected,
            "wildcard_adapter.fa",
        )
        with open(path) as wct:
            lines = wct.readlines()
        lines = [line.strip() for line in lines]
        assert lines == ["AAA 1", "GGG 2", "CCC 3b", "TTT 4b"]


def test_wildcard_N(tmp_path):
    """test 'N' wildcard matching with no allowed errors"""
    run(
        tmp_path,
        "-e 0 -a GGGGGGG --match-read-wildcards",
        "wildcardN.fa",
        "wildcardN.fa",
    )


def test_illumina_adapter_wildcard(tmp_path):
    run(
        tmp_path,
        "-a VCCGAMCYUCKHRKDCUBBCNUWNSGHCGU",
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_adapter_front(tmp_path):
    """test adapter in front"""
    run(tmp_path, "--front ADAPTER -N", "examplefront.fa", "example.fa")


def test_literal_N(tmp_path):
    """test matching literal 'N's"""
    run(tmp_path, "-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta")


def test_literal_N2(tmp_path):
    run(tmp_path, "-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta")


def test_literal_N_brace_notation(tmp_path):
    """test matching literal 'N's"""
    run(tmp_path, "-N -e 0.2 -a N{14}", "trimN3.fasta", "trimN3.fasta")


def test_literal_N2_brace_notation(tmp_path):
    run(tmp_path, "-N -O 1 -g N{14}", "trimN5.fasta", "trimN5.fasta")


def test_anchored_front(tmp_path):
    run(tmp_path, "-g ^FRONTADAPT -N", "anchored.fasta", "anchored.fasta")


def test_anchored_front_ellipsis_notation(tmp_path):
    run(tmp_path, "-a FRONTADAPT... -N", "anchored.fasta", "anchored.fasta")


def test_anchored_back(tmp_path):
    run(tmp_path, "-a BACKADAPTER$ -N", "anchored-back.fasta", "anchored-back.fasta")


def test_anchored_back_no_indels(tmp_path):
    run(
        tmp_path,
        "-a BACKADAPTER$ -N --no-indels",
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_no_indels(tmp_path):
    run(
        tmp_path,
        "-a TTAGACATAT -g GAGATTGCCA --no-indels",
        "no_indels.fasta",
        "no_indels.fasta",
    )


def test_issue_46(tmp_path):
    """issue 46 - IndexError with --wildcard-file"""
    path = tmp_path / "wildcardtmp.txt"
    run(
        tmp_path,
        "--anywhere=AACGTN --wildcard-file={0}".format(path),
        "issue46.fasta",
        "issue46.fasta",
    )


def test_strip_suffix(tmp_path):
    run(
        tmp_path,
        "--strip-suffix _sequence -a XXXXXXX",
        "stripped.fasta",
        "simple.fasta",
    )


def test_info_file(tmp_path):
    # The true adapter sequence in the illumina.fastq.gz data set is
    # GCCTAACTTCTTAGACTGCCTTAAGGACGT (fourth base is different)
    #
    path = str(tmp_path / "infotmp.txt")
    run(
        tmp_path,
        ["--info-file", path, "-a", "adapt=GCCGAACTTCTTAGACTGCCTTAAGGACGT"],
        "illumina.fastq",
        "illumina.fastq.gz",
    )
    assert files_equal(cutpath("illumina.info.txt"), path)


def test_info_file_times(tmp_path):
    path = tmp_path / "infotmp.txt"
    run(
        tmp_path,
        [
            "--info-file",
            str(path),
            "--times",
            "2",
            "-a",
            "adapt=GCCGAACTTCTTA",
            "-a",
            "adapt2=GACTGCCTTAAGGACGT",
        ],
        "illumina5.fastq",
        "illumina5.fastq",
    )
    assert files_equal(cutpath("illumina5.info.txt"), str(path))


def test_info_file_fasta(tmp_path):
    path = tmp_path / "infotmp.txt"
    # Just make sure that it runs
    run(
        tmp_path,
        [
            "--info-file",
            str(path),
            "-a",
            "TTAGACATAT",
            "-g",
            "GAGATTGCCA",
            "--no-indels",
        ],
        "no_indels.fasta",
        "no_indels.fasta",
    )


def test_named_adapter(tmp_path):
    run(
        tmp_path,
        "-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_adapter_with_U(tmp_path):
    run(
        tmp_path,
        "-a GCCGAACUUCUUAGACUGCCUUAAGGACGU",
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_no_trim(tmp_path):
    """--no-trim"""
    run(
        tmp_path,
        "--no-trim --discard-untrimmed -a CCCTAGTTAAAC",
        "no-trim.fastq",
        "small.fastq",
    )


def test_bzip2(tmp_path):
    """test bzip2 support"""
    run(tmp_path, "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.bz2")


try:
    import lzma

    def test_xz(tmp_path):
        """test xz support"""
        run(tmp_path, "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.xz")

except ImportError:
    pass


def test_qualfile_only(tmp_path):
    with raises(SystemExit), redirect_stderr():
        execute_cli(["-sq", datapath("E3M.qual")])


def test_no_args(tmp_path):
    with redirect_stderr():
        assert execute_cli() != 0


def test_two_fastqs(tmp_path):
    with raises(SystemExit), redirect_stderr():
        execute_cli(
            ["-pe1", datapath("paired.1.fastq"), "-pe2", datapath("paired.2.fastq")]
        )


def test_anchored_no_indels(tmp_path):
    """anchored 5' adapter, mismatches only (no indels)"""
    run(
        tmp_path,
        "-g ^TTAGACATAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
    )


def test_anchored_no_indels_wildcard_read(tmp_path):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the read count as matches"""
    run(
        tmp_path,
        "-g ^TTAGACATAT --match-read-wildcards --no-indels -e 0.1",
        "anchored_no_indels_wildcard.fasta",
        "anchored_no_indels.fasta",
    )


def test_anchored_no_indels_wildcard_adapt(tmp_path):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the adapter count as matches"""
    run(
        tmp_path,
        "-g ^TTAGACANAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
    )


def test_unconditional_cut_front(tmp_path):
    run(tmp_path, "-u 5", "unconditional-front.fastq", "small.fastq")


def test_unconditional_cut_back(tmp_path):
    run(tmp_path, "-u -5", "unconditional-back.fastq", "small.fastq")


def test_unconditional_cut_both(tmp_path):
    run(tmp_path, "-u -5 -u 5", "unconditional-both.fastq", "small.fastq")


def test_untrimmed_output(tmp_path):
    path = tmp_path / "untrimmed.tmp.fastq"
    run(
        tmp_path,
        ["-a", "TTAGACATATCTCCGTCG", "--untrimmed-output", str(path)],
        "small.trimmed.fastq",
        "small.fastq",
    )
    assert files_equal(cutpath("small.untrimmed.fastq"), str(path))


def test_adapter_file(tmp_path):
    run(
        tmp_path,
        "-a file:" + datapath("adapter.fasta"),
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_adapter_file_5p_anchored(tmp_path):
    run(
        tmp_path,
        "-N -g file:" + datapath("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
    )


def test_adapter_file_3p_anchored(tmp_path):
    run(
        tmp_path,
        "-N -a file:" + datapath("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_adapter_file_5p_anchored_no_indels(tmp_path):
    run(
        tmp_path,
        "-N --no-indels -g file:" + datapath("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
    )


def test_adapter_file_3p_anchored_no_indels(tmp_path):
    run(
        tmp_path,
        "-N --no-indels -a file:" + datapath("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_demultiplex(tmp_path):
    multiout = os.path.join(
        os.path.dirname(__file__), "data", "tmp-demulti.{name}.fasta"
    )
    params = [
        "-a",
        "first=AATTTCAGGAATT",
        "-a",
        "second=GTTCTCTAGTTCT",
        "-o",
        multiout,
        "-se",
        datapath("twoadapters.fasta"),
    ]
    command = get_command("trim")
    result = command.execute(params)
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result[0] == 0
    assert files_equal(
        cutpath("twoadapters.first.fasta"), multiout.format(name="first")
    )
    assert files_equal(
        cutpath("twoadapters.second.fasta"), multiout.format(name="second")
    )
    assert files_equal(
        cutpath("twoadapters.unknown.fasta"), multiout.format(name="unknown")
    )
    os.remove(multiout.format(name="first"))
    os.remove(multiout.format(name="second"))
    os.remove(multiout.format(name="unknown"))


def test_max_n(tmp_path):
    run(tmp_path, "--max-n 0", "maxn0.fasta", "maxn.fasta")
    run(tmp_path, "--max-n 1", "maxn1.fasta", "maxn.fasta")
    run(tmp_path, "--max-n 2", "maxn2.fasta", "maxn.fasta")
    run(tmp_path, "--max-n 0.2", "maxn0.2.fasta", "maxn.fasta")
    run(tmp_path, "--max-n 0.4", "maxn0.4.fasta", "maxn.fasta")


def test_quiet_is_quiet(tmp_path):
    captured_standard_output = StringIO()
    captured_standard_error = StringIO()
    try:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = captured_standard_output
        sys.stderr = captured_standard_error
        execute_cli(
            [
                "-o",
                "/dev/null",
                "--quiet",
                "-a",
                "XXXX",
                "-se",
                datapath("illumina.fastq.gz"),
            ]
        )
    finally:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
    print(captured_standard_output.getvalue())
    assert captured_standard_output.getvalue() == ""
    assert captured_standard_error.getvalue() == ""


def test_nextseq(tmp_path):
    run(tmp_path, "--nextseq-trim 22", "nextseq.fastq", "nextseq.fastq")


def test_linked(tmp_path):
    run(tmp_path, "-a AAAAAAAAAA...TTTTTTTTTT", "linked.fasta", "linked.fasta")


def test_fasta(tmp_path):
    run(tmp_path, "-a TTAGACATATCTCCGTCG", "small.fasta", "small.fastq")


def test_custom_bisulfite_1(tmp_path):
    run(
        tmp_path,
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,1",
        "small.fastq",
        "small.fastq",
    )


def test_custom_bisulfite_2(tmp_path):
    run(
        tmp_path,
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 15,15,1,1",
        "small_mincut1.fastq",
        "small.fastq",
    )


def test_custom_bisulfite_3(tmp_path):
    run(
        tmp_path,
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,0",
        "small_mincut2.fastq",
        "small.fastq",
    )


def test_custom_bisulfite_4(tmp_path):
    run(
        tmp_path,
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,0,0",
        "small_mincut3.fastq",
        "small.fastq",
    )


@skipIf(
    no_internet("https://ncbi.nlm.nih.gov") or no_import("srastream"),
    "No internet connection or srastream not importable",
)
def test_sra(tmp_path):
    run(
        tmp_path,
        "-b CTGGAGTTCAGACGTGTGCTCT --max-reads 100",
        "SRR2040662_trimmed.fq",
        sra_accn="SRR2040662",
    )
