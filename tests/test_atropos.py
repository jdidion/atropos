# TODO:
#  test with the --output option
#  test reading from standard input
import os

import pytest
from xphyle import open_

from atropos.console import execute_cli
from atropos.commands.trim.console import TrimCommandConsole

from .utils import (
    run,
    files_equal,
    datapath,
    cutpath,
    redirect_stderr,
    no_internet,
    no_import,
    intercept_stdout,
    intercept_stderr,
)

try:
    import lzma
except ImportError:
    lzma = None


def test_example(tmp_path_factory):
    run("-N -b ADAPTER", "example.fa", "example.fa", tmp_path_factory=tmp_path_factory)


def test_example_stdout():
    run('-N -b ADAPTER', 'example.fa', 'example.fa', stdout=True)


def test_small(tmp_path_factory):
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_empty(tmp_path_factory):
    """empty input"""
    run(
        "-a TTAGACATATCTCCGTCG", "empty.fastq", "empty.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_newlines(tmp_path_factory):
    """DOS/Windows newlines"""
    run(
        "-e 0.12 -b TTAGACATATCTCCGTCG", "dos.fastq", "dos.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_lowercase(tmp_path_factory):
    """lowercase adapter"""
    run(
        "-b ttagacatatctccgtcg", "lowercase.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_rest(tmp_path_factory):
    """-r/--rest-file"""
    rest_tmp = tmp_path_factory.mktemp("rest.tmp")
    run(
        ["-b", "ADAPTER", "-N", "-r", rest_tmp], "rest.fa", "rest.fa",
        tmp_path_factory=tmp_path_factory
    )
    assert files_equal(datapath("rest.txt"), rest_tmp)


def test_restfront(tmp_path_factory):
    path = tmp_path_factory.mktemp("rest.txt")
    run(
        ["-g", "ADAPTER", "-N", "-r", path], "restfront.fa", "rest.fa",
        tmp_path_factory=tmp_path_factory
    )
    assert files_equal(datapath("restfront.txt"), path)


def test_discard(tmp_path_factory):
    """--discard"""
    run(
        "-b TTAGACATATCTCCGTCG --discard", "discard.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_discard_untrimmed(tmp_path_factory):
    """--discard-untrimmed"""
    run(
        "-b CAAGAT --discard-untrimmed", "discard-untrimmed.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_plus(tmp_path_factory):
    """test if sequence name after the "+" is retained"""
    run(
        "-e 0.12 -b TTAGACATATCTCCGTCG", "plus.fastq", "plus.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_extensiontxtgz(tmp_path_factory):
    """automatic recognition of "_sequence.txt.gz" extension"""
    run(
        "-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz",
        tmp_path_factory=tmp_path_factory
    )


def test_format(tmp_path_factory):
    """the -f/--format parameter"""
    run(
        "-f fastq -b TTAGACATATCTCCGTCG", "small.fastq", "small.myownextension",
        tmp_path_factory=tmp_path_factory

    )


def test_minimum_length(tmp_path_factory):
    """-m/--minimum-length"""
    run(
        "-c -m 5 -a 330201030313112312", "minlen.fa", "lengths.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_too_short(tmp_path_factory):
    """--too-short-output"""
    run(
        "-c -m 5 -a 330201030313112312 --too-short-output tooshort.tmp.fa",
        "minlen.fa",
        "lengths.fa",
        tmp_path_factory=tmp_path_factory,
    )
    assert files_equal(datapath("tooshort.fa"), "tooshort.tmp.fa")
    os.remove("tooshort.tmp.fa")


def test_too_short_no_primer(tmp_path_factory):
    """--too-short-output and --trim-primer"""
    run(
        "-c -m 5 -a 330201030313112312 --trim-primer --too-short-output tooshort.tmp.fa",
        "minlen.noprimer.fa",
        "lengths.fa",
        tmp_path_factory=tmp_path_factory,
    )
    assert files_equal(datapath("tooshort.noprimer.fa"), "tooshort.tmp.fa")
    os.remove("tooshort.tmp.fa")


def test_maximum_length(tmp_path_factory):
    """-M/--maximum-length"""
    run(
        "-c -M 5 -a 330201030313112312", "maxlen.fa", "lengths.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_too_long(tmp_path_factory):
    """--too-long-output"""
    run(
        "-c -M 5 --too-long-output toolong.tmp.fa -a 330201030313112312",
        "maxlen.fa",
        "lengths.fa",
        tmp_path_factory=tmp_path_factory,
    )
    assert files_equal(datapath("toolong.fa"), "toolong.tmp.fa")
    os.remove("toolong.tmp.fa")


def test_length_tag(tmp_path_factory):
    """454 data; -n and --length-tag"""
    run(
        "-n 3 -e 0.1 --length-tag length= "
        "-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "
        "-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
        "454.fa",
        "454.fa",
        tmp_path_factory=tmp_path_factory,
    )


def test_overlap_a(tmp_path_factory):
    """-O/--overlap with -a (-c omitted on purpose)"""
    run(
        "-O 10 -a 330201030313112312 -e 0.0 -N", "overlapa.fa", "overlapa.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_overlap_b(tmp_path_factory):
    """-O/--overlap with -b"""
    run(
        "-O 10 -b TTAGACATATCTCCGTCG -N", "overlapb.fa", "overlapb.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_qualtrim(tmp_path_factory):
    """-q with low qualities"""
    run(
        "-q 10 -a XXXXXX", "lowqual.fastq", "lowqual.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_qualbase(tmp_path_factory):
    """-q with low qualities, using ascii(quality+64) encoding"""
    run(
        "-q 10 --quality-base 64 -a XXXXXX", "illumina64.fastq", "illumina64.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_quality_trim_only(tmp_path_factory):
    """only trim qualities, do not remove adapters"""
    run(
        "-q 10 --quality-base 64", "illumina64.fastq", "illumina64.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_twoadapters(tmp_path_factory):
    """two adapters"""
    run(
        "-a AATTTCAGGAATT -a GTTCTCTAGTTCT", "twoadapters.fasta", "twoadapters.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_polya(tmp_path_factory):
    """poly-A tails"""
    run(
        "-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "polya.fasta",
        "polya.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_polya_brace_notation(tmp_path_factory):
    """poly-A tails"""
    run(
        "-m 24 -O 10 -a A{35}", "polya.fasta", "polya.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_mask_adapter(tmp_path_factory):
    """mask adapter with N (reads maintain the same length)"""
    run(
        "-b CAAG -n 3 --mask-adapter", "anywhere_repeat.fastq", "anywhere_repeat.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_gz_multiblock(tmp_path_factory):
    """compressed gz file with multiple blocks (created by concatenating two .gz files)
    """
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz",
        tmp_path_factory=tmp_path_factory
    )


def test_suffix(tmp_path_factory):
    """-y/--suffix parameter, combined with _F3"""
    run(
        "-c -e 0.12 -a 1=330201030313112312 -y _my_suffix_{name} --strip-f3",
        "suffix.fastq",
        "solid.csfasta",
        qualfile="solid.qual",
        tmp_path_factory=tmp_path_factory,
    )


def test_read_wildcard(tmp_path_factory):
    """test wildcards in reads"""
    run(
        "--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_adapter_wildcard(tmp_path_factory):
    """wildcards in adapter"""
    for adapter_type, expected in (
        ("-a", "wildcard_adapter.fa"),
        ("-b", "wildcard_adapter_anywhere.fa"),
    ):
        wildcardtmp = tmp_path_factory.mktemp("wildcardtmp.txt")
        run(
            "--wildcard-file {0} {1} ACGTNNNACGT".format(wildcardtmp, adapter_type),
            expected,
            "wildcard_adapter.fa",
            tmp_path_factory=tmp_path_factory
        )
        with open_(wildcardtmp) as wct:
            lines = wct.readlines()
        lines = [line.strip() for line in lines]
        assert lines == ["AAA 1", "GGG 2", "CCC 3b", "TTT 4b"]


def test_wildcard_n(tmp_path_factory):
    """test 'N' wildcard matching with no allowed errors"""
    run(
        "-e 0 -a GGGGGGG --match-read-wildcards", "wildcardN.fa", "wildcardN.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_illumina_adapter_wildcard(tmp_path_factory):
    run(
        "-a VCCGAMCYUCKHRKDCUBBCNUWNSGHCGU", "illumina.fastq", "illumina.fastq.gz",
        tmp_path_factory=tmp_path_factory
    )


def test_adapter_front(tmp_path_factory):
    """test adapter in front"""
    run(
        "--front ADAPTER -N", "examplefront.fa", "example.fa",
        tmp_path_factory=tmp_path_factory
    )


def test_literal_n(tmp_path_factory):
    """test matching literal 'N's"""
    run(
        "-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_literal_n2(tmp_path_factory):
    run(
        "-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_literal_n_brace_notation(tmp_path_factory):
    """test matching literal 'N's"""
    run(
        "-N -e 0.2 -a N{14}", "trimN3.fasta", "trimN3.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_literal_n2_brace_notation(tmp_path_factory):
    run(
        "-N -O 1 -g N{14}", "trimN5.fasta", "trimN5.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_anchored_front(tmp_path_factory):
    run(
        "-g ^FRONTADAPT -N", "anchored.fasta", "anchored.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_anchored_front_ellipsis_notation(tmp_path_factory):
    run(
        "-a FRONTADAPT... -N", "anchored.fasta", "anchored.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_anchored_back(tmp_path_factory):
    run(
        "-a BACKADAPTER$ -N", "anchored-back.fasta", "anchored-back.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_anchored_back_no_indels(tmp_path_factory):
    run(
        "-a BACKADAPTER$ -N --no-indels", "anchored-back.fasta", "anchored-back.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_no_indels(tmp_path_factory):
    run(
        "-a TTAGACATAT -g GAGATTGCCA --no-indels", "no_indels.fasta", "no_indels.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_issue_46(tmp_path_factory):
    """issue 46 - IndexError with --wildcard-file"""
    wildcardtmp = tmp_path_factory.mktemp("wildcardtmp.txt")
    run(
        "--anywhere=AACGTN --wildcard-file={0}".format(wildcardtmp),
        "issue46.fasta",
        "issue46.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_strip_suffix(tmp_path_factory):
    run(
        "--strip-suffix _sequence -a XXXXXXX", "stripped.fasta", "simple.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_info_file(tmp_path_factory):
    # The true adapter sequence in the illumina.fastq.gz data set is
    # GCCTAACTTCTTAGACTGCCTTAAGGACGT (fourth base is different)
    #
    infotmp = tmp_path_factory.mktemp("infotmp.txt")
    run(
        ["--info-file", infotmp, "-a", "adapt=GCCGAACTTCTTAGACTGCCTTAAGGACGT"],
        "illumina.fastq",
        "illumina.fastq.gz",
        tmp_path_factory=tmp_path_factory
    )
    assert files_equal(cutpath("illumina.info.txt"), infotmp)


def test_info_file_times(tmp_path_factory):
    infotmp = tmp_path_factory.mktemp("infotmp.txt")
    run(
        [
            "--info-file",
            infotmp,
            "--times",
            "2",
            "-a",
            "adapt=GCCGAACTTCTTA",
            "-a",
            "adapt2=GACTGCCTTAAGGACGT",
        ],
        "illumina5.fastq",
        "illumina5.fastq",
        tmp_path_factory=tmp_path_factory
    )
    assert files_equal(cutpath("illumina5.info.txt"), infotmp)


def test_info_file_fasta(tmp_path_factory):
    infotmp = tmp_path_factory.mktemp("infotmp.txt")
    # Just make sure that it runs
    run(
        [
            "--info-file",
            infotmp,
            "-a",
            "TTAGACATAT",
            "-g",
            "GAGATTGCCA",
            "--no-indels",
        ],
        "no_indels.fasta",
        "no_indels.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_named_adapter(tmp_path_factory):
    run(
        "-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz",
        tmp_path_factory=tmp_path_factory
    )


def test_adapter_with_u(tmp_path_factory):
    run(
        "-a GCCGAACUUCUUAGACUGCCUUAAGGACGU", "illumina.fastq", "illumina.fastq.gz",
        tmp_path_factory=tmp_path_factory
    )


def test_no_trim(tmp_path_factory):
    """ --no-trim """
    run(
        "--no-trim --discard-untrimmed -a CCCTAGTTAAAC", "no-trim.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_bzip2(tmp_path_factory):
    """test bzip2 support"""
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.bz2",
        tmp_path_factory=tmp_path_factory
    )


@pytest.mark.skipif(lzma is None)
def test_xz(tmp_path_factory):
    """test xz support"""
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.xz",
        tmp_path_factory=tmp_path_factory
    )


def test_qualfile_only():
    with pytest.raises(SystemExit), redirect_stderr():
        execute_cli(["-sq", datapath("E3M.qual")])


def test_no_args():
    with redirect_stderr():
        assert execute_cli() != 0


def test_anchored_no_indels(tmp_path_factory):
    """anchored 5' adapter, mismatches only (no indels)"""
    run(
        "-g ^TTAGACATAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_anchored_no_indels_wildcard_read(tmp_path_factory):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the read
    count as matches.
    """
    run(
        "-g ^TTAGACATAT --match-read-wildcards --no-indels -e 0.1",
        "anchored_no_indels_wildcard.fasta",
        "anchored_no_indels.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_anchored_no_indels_wildcard_adapt(tmp_path_factory):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the adapter
    count as matches.
    """
    run(
        "-g ^TTAGACANAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_unconditional_cut_front(tmp_path_factory):
    run(
        "-u 5", "unconditional-front.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_unconditional_cut_back(tmp_path_factory):
    run(
        "-u -5", "unconditional-back.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_unconditional_cut_both(tmp_path_factory):
    run(
        "-u -5 -u 5", "unconditional-both.fastq", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_untrimmed_output(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("untrimmed.tmp.fastq")
    run(
        ["-a", "TTAGACATATCTCCGTCG", "--untrimmed-output", tmp],
        "small.trimmed.fastq",
        "small.fastq",
        tmp_path_factory=tmp_path_factory
    )
    assert files_equal(cutpath("small.untrimmed.fastq"), tmp)


def test_adapter_file(tmp_path_factory):
    run(
        "-a file:" + datapath("adapter.fasta"), "illumina.fastq", "illumina.fastq.gz",
        tmp_path_factory=tmp_path_factory
    )


def test_adapter_file_5p_anchored(tmp_path_factory):
    run(
        "-N -g file:" + datapath("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_adapter_file_3p_anchored(tmp_path_factory):
    run(
        "-N -a file:" + datapath("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_adapter_file_5p_anchored_no_indels(tmp_path_factory):
    run(
        "-N --no-indels -g file:" + datapath("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_adapter_file_3p_anchored_no_indels(tmp_path_factory):
    run(
        "-N --no-indels -a file:" + datapath("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
        tmp_path_factory=tmp_path_factory,
    )


def test_demultiplex():
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
    result = TrimCommandConsole.execute(params)
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


def test_max_n(tmp_path_factory):
    run("--max-n 0", "maxn0.fasta", "maxn.fasta", tmp_path_factory=tmp_path_factory)
    run("--max-n 1", "maxn1.fasta", "maxn.fasta", tmp_path_factory=tmp_path_factory)
    run("--max-n 2", "maxn2.fasta", "maxn.fasta", tmp_path_factory=tmp_path_factory)
    run("--max-n 0.2", "maxn0.2.fasta", "maxn.fasta", tmp_path_factory=tmp_path_factory)
    run("--max-n 0.4", "maxn0.4.fasta", "maxn.fasta", tmp_path_factory=tmp_path_factory)


def test_quiet_is_quiet():
    with intercept_stdout() as captured_stdout, intercept_stderr() as captured_stderr:
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
        assert captured_stdout.getvalue() == ""
        assert captured_stderr.getvalue() == ""


def test_nextseq(tmp_path_factory):
    run(
        "--nextseq-trim 22", "nextseq.fastq", "nextseq.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_linked(tmp_path_factory):
    run(
        "-a AAAAAAAAAA...TTTTTTTTTT", "linked.fasta", "linked.fasta",
        tmp_path_factory=tmp_path_factory
    )


def test_fasta(tmp_path_factory):
    run(
        "-a TTAGACATATCTCCGTCG", "small.fasta", "small.fastq",
        tmp_path_factory=tmp_path_factory
    )


def test_custom_bisulfite_1(tmp_path_factory):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,1", "small.fastq",
        "small.fastq", tmp_path_factory=tmp_path_factory
    )


def test_custom_bisulfite_2(tmp_path_factory):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 15,15,1,1",
        "small_mincut1.fastq",
        "small.fastq",
        tmp_path_factory=tmp_path_factory,
    )


def test_custom_bisulfite_3(tmp_path_factory):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,0",
        "small_mincut2.fastq",
        "small.fastq",
        tmp_path_factory=tmp_path_factory,
    )


def test_custom_bisulfite_4(tmp_path_factory):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,0,0",
        "small_mincut3.fastq",
        "small.fastq",
        tmp_path_factory=tmp_path_factory,
    )


@pytest.mark.skipif(
    no_internet("https://ncbi.nlm.nih.gov") or no_import("srastream"),
    "No internet connection or srastream not importable",
)
def test_sra(tmp_path_factory):
    run(
        "-b CTGGAGTTCAGACGTGTGCTCT --max-reads 100",
        "SRR2040662_trimmed.fq",
        sra_accn="SRR2040662",
        tmp_path_factory=tmp_path_factory,
    )
