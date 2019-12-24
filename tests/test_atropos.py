# TODO:
#  test with the --output option
#  test reading from standard input
import os

import pytest
from xphyle import open_

from atropos.console import execute_cli
from atropos.commands.trim.console import TrimCommandConsole
from atropos.utils import ReturnCode

from .utils import (
    run,
    assert_files_equal,
    datapath,
    datapathstr,
    dataurl,
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


def test_example(tmp_path):
    run("-N -b ADAPTER", "example.fa", "example.fa", output_dir=tmp_path)


def test_example_stdout(tmp_path):
    run('-N -b ADAPTER', 'example.fa', 'example.fa', output_dir=tmp_path, stdout=True)


def test_small(tmp_path):
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_empty(tmp_path):
    """empty input"""
    run(
        "-a TTAGACATATCTCCGTCG", "empty.fastq", "empty.fastq",
        output_dir=tmp_path
    )


def test_newlines(tmp_path):
    """DOS/Windows newlines"""
    run(
        "-e 0.12 -b TTAGACATATCTCCGTCG", "dos.fastq", "dos.fastq",
        output_dir=tmp_path
    )


def test_lowercase(tmp_path):
    """lowercase adapter"""
    run(
        "-b ttagacatatctccgtcg", "lowercase.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_rest(tmp_path):
    """-r/--rest-file"""
    rest_tmp = tmp_path / "rest.tmp"
    run(
        ["-b", "ADAPTER", "-N", "-r", rest_tmp], "rest.fa", "rest.fa",
        output_dir=tmp_path
    )
    assert_files_equal(datapath("rest.txt"), rest_tmp)


def test_restfront(tmp_path):
    path = tmp_path / "rest.txt"
    run(
        ["-g", "ADAPTER", "-N", "-r", path], "restfront.fa", "rest.fa",
        output_dir=tmp_path
    )
    assert_files_equal(datapath("restfront.txt"), path)


def test_discard_trimmed(tmp_path):
    """--discard"""
    run(
        "-b TTAGACATATCTCCGTCG --discard-trimmed", "discard.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_discard_untrimmed(tmp_path):
    """--discard-untrimmed"""
    run(
        "-b CAAGAT --discard-untrimmed", "discard-untrimmed.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_plus(tmp_path):
    """test if sequence name after the "+" is retained"""
    run(
        "-e 0.12 -b TTAGACATATCTCCGTCG", "plus.fastq", "plus.fastq",
        output_dir=tmp_path
    )


def test_extensiontxtgz(tmp_path):
    """automatic recognition of "_sequence.txt.gz" extension"""
    run(
        "-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz",
        output_dir=tmp_path
    )


def test_format(tmp_path):
    """the -f/--format parameter"""
    run(
        "-f fastq -b TTAGACATATCTCCGTCG", "small.fastq", "small.myownextension",
        output_dir=tmp_path

    )


def test_minimum_length(tmp_path):
    """-m/--minimum-length"""
    run(
        "-c -m 5 -a 330201030313112312", "minlen.fa", "lengths.fa",
        output_dir=tmp_path
    )


def test_too_short(tmp_path):
    """--too-short-output"""
    run(
        "-c -m 5 -a 330201030313112312 --too-short-output tooshort.tmp.fa",
        "minlen.fa",
        "lengths.fa",
        output_dir=tmp_path,
    )
    assert_files_equal(datapath("tooshort.fa"), "tooshort.tmp.fa")
    os.remove("tooshort.tmp.fa")


def test_too_short_no_primer(tmp_path):
    """--too-short-output and --trim-primer"""
    run(
        "-c -m 5 -a 330201030313112312 --trim-primer --too-short-output tooshort.tmp.fa",
        "minlen.noprimer.fa",
        "lengths.fa",
        output_dir=tmp_path,
    )
    assert_files_equal(datapath("tooshort.noprimer.fa"), "tooshort.tmp.fa")
    os.remove("tooshort.tmp.fa")


def test_maximum_length(tmp_path):
    """-M/--maximum-length"""
    run(
        "-c -M 5 -a 330201030313112312", "maxlen.fa", "lengths.fa",
        output_dir=tmp_path
    )


def test_too_long(tmp_path):
    """--too-long-output"""
    run(
        "-c -M 5 --too-long-output toolong.tmp.fa -a 330201030313112312",
        "maxlen.fa",
        "lengths.fa",
        output_dir=tmp_path,
    )
    assert_files_equal(datapath("toolong.fa"), "toolong.tmp.fa")
    os.remove("toolong.tmp.fa")


def test_length_tag(tmp_path):
    """454 data; -n and --length-tag"""
    run(
        "-n 3 -e 0.1 --length-tag length= "
        "-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "
        "-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
        "454.fa",
        "454.fa",
        output_dir=tmp_path,
    )


def test_overlap_a(tmp_path):
    """-O/--overlap with -a (-c omitted on purpose)"""
    run(
        "-O 10 -a 330201030313112312 -e 0.0 -N", "overlapa.fa", "overlapa.fa",
        output_dir=tmp_path
    )


def test_overlap_b(tmp_path):
    """-O/--overlap with -b"""
    run(
        "-O 10 -b TTAGACATATCTCCGTCG -N", "overlapb.fa", "overlapb.fa",
        output_dir=tmp_path
    )


def test_qualtrim(tmp_path):
    """-q with low qualities"""
    run(
        "-q 10 -a XXXXXX", "lowqual.fastq", "lowqual.fastq",
        output_dir=tmp_path
    )


def test_qualbase(tmp_path):
    """-q with low qualities, using ascii(quality+64) encoding"""
    run(
        "-q 10 --quality-base 64 -a XXXXXX", "illumina64.fastq", "illumina64.fastq",
        output_dir=tmp_path
    )


def test_quality_trim_only(tmp_path):
    """only trim qualities, do not remove adapters"""
    run(
        "-q 10 --quality-base 64", "illumina64.fastq", "illumina64.fastq",
        output_dir=tmp_path
    )


def test_twoadapters(tmp_path):
    """two adapters"""
    run(
        "-a AATTTCAGGAATT -a GTTCTCTAGTTCT", "twoadapters.fasta", "twoadapters.fasta",
        output_dir=tmp_path
    )


def test_polya(tmp_path):
    """poly-A tails"""
    run(
        "-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "polya.fasta",
        "polya.fasta",
        output_dir=tmp_path,
    )


def test_polya_brace_notation(tmp_path):
    """poly-A tails"""
    run(
        "-m 24 -O 10 -a A{35}", "polya.fasta", "polya.fasta",
        output_dir=tmp_path
    )


def test_mask_adapter(tmp_path):
    """mask adapter with N (reads maintain the same length)"""
    run(
        "-b CAAG -n 3 --mask-adapter", "anywhere_repeat.fastq", "anywhere_repeat.fastq",
        output_dir=tmp_path
    )


def test_gz_multiblock(tmp_path):
    """compressed gz file with multiple blocks (created by concatenating two .gz files)
    """
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz",
        output_dir=tmp_path
    )


def test_suffix(tmp_path):
    """-y/--suffix parameter, combined with _F3"""
    run(
        "-c -e 0.12 -a 1=330201030313112312 -y _my_suffix_{name} --strip-f3",
        "suffix.fastq",
        "solid.csfasta",
        qualfile="solid.qual",
        output_dir=tmp_path,
    )


def test_read_wildcard(tmp_path):
    """test wildcards in reads"""
    run(
        "--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa",
        output_dir=tmp_path
    )


def test_adapter_wildcard(tmp_path):
    """wildcards in adapter"""
    for adapter_type, expected in (
        ("-a", "wildcard_adapter.fa"),
        ("-b", "wildcard_adapter_anywhere.fa"),
    ):
        wildcardtmp = tmp_path / "wildcardtmp.txt"
        run(
            "--wildcard-file {0} {1} ACGTNNNACGT".format(wildcardtmp, adapter_type),
            expected,
            "wildcard_adapter.fa",
            output_dir=tmp_path
        )
        with open_(wildcardtmp) as wct:
            lines = wct.readlines()
        lines = [line.strip() for line in lines]
        assert lines == ["AAA 1", "GGG 2", "CCC 3b", "TTT 4b"]


def test_wildcard_n(tmp_path):
    """test 'N' wildcard matching with no allowed errors"""
    run(
        "-e 0 -a GGGGGGG --match-read-wildcards", "wildcardN.fa", "wildcardN.fa",
        output_dir=tmp_path
    )


def test_illumina_adapter_wildcard(tmp_path):
    run(
        "-a VCCGAMCYUCKHRKDCUBBCNUWNSGHCGU", "illumina.fastq", "illumina.fastq.gz",
        output_dir=tmp_path
    )


def test_adapter_front(tmp_path):
    """test adapter in front"""
    run(
        "--front ADAPTER -N", "examplefront.fa", "example.fa",
        output_dir=tmp_path
    )


def test_literal_n(tmp_path):
    """test matching literal 'N's"""
    run(
        "-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta",
        output_dir=tmp_path
    )


def test_literal_n2(tmp_path):
    run(
        "-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta",
        output_dir=tmp_path
    )


def test_literal_n_brace_notation(tmp_path):
    """test matching literal 'N's"""
    run(
        "-N -e 0.2 -a N{14}", "trimN3.fasta", "trimN3.fasta",
        output_dir=tmp_path
    )


def test_literal_n2_brace_notation(tmp_path):
    run(
        "-N -O 1 -g N{14}", "trimN5.fasta", "trimN5.fasta",
        output_dir=tmp_path
    )


def test_anchored_front(tmp_path):
    run(
        "-g ^FRONTADAPT -N", "anchored.fasta", "anchored.fasta",
        output_dir=tmp_path
    )


def test_anchored_front_ellipsis_notation(tmp_path):
    run(
        "-a FRONTADAPT... -N", "anchored.fasta", "anchored.fasta",
        output_dir=tmp_path
    )


def test_anchored_back(tmp_path):
    run(
        "-a BACKADAPTER$ -N", "anchored-back.fasta", "anchored-back.fasta",
        output_dir=tmp_path
    )


def test_anchored_back_no_indels(tmp_path):
    run(
        "-a BACKADAPTER$ -N --no-indels", "anchored-back.fasta", "anchored-back.fasta",
        output_dir=tmp_path
    )


def test_no_indels(tmp_path):
    run(
        "-a TTAGACATAT -g GAGATTGCCA --no-indels", "no_indels.fasta", "no_indels.fasta",
        output_dir=tmp_path
    )


def test_issue_46(tmp_path):
    """issue 46 - IndexError with --wildcard-file"""
    wildcardtmp = tmp_path / "wildcardtmp.txt"
    run(
        "--anywhere=AACGTN --wildcard-file={0}".format(wildcardtmp),
        "issue46.fasta",
        "issue46.fasta",
        output_dir=tmp_path,
    )


def test_strip_suffix(tmp_path):
    run(
        "--strip-suffix _sequence -a XXXXXXX", "stripped.fasta", "simple.fasta",
        output_dir=tmp_path
    )


def test_info_file(tmp_path):
    # The true adapter sequence in the illumina.fastq.gz data set is
    # GCCTAACTTCTTAGACTGCCTTAAGGACGT (fourth base is different)
    #
    infotmp = tmp_path / "infotmp.txt"
    run(
        ["--info-file", infotmp, "-a", "adapt=GCCGAACTTCTTAGACTGCCTTAAGGACGT"],
        "illumina.fastq",
        "illumina.fastq.gz",
        output_dir=tmp_path
    )
    assert_files_equal(cutpath("illumina.info.txt"), infotmp)


def test_info_file_times(tmp_path):
    infotmp = tmp_path / "infotmp.txt"
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
        output_dir=tmp_path
    )
    assert_files_equal(cutpath("illumina5.info.txt"), infotmp)


def test_info_file_fasta(tmp_path):
    infotmp = tmp_path / "infotmp.txt"
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
        output_dir=tmp_path
    )


def test_named_adapter(tmp_path):
    run(
        "-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz",
        output_dir=tmp_path
    )


def test_adapter_with_u(tmp_path):
    run(
        "-a GCCGAACUUCUUAGACUGCCUUAAGGACGU", "illumina.fastq", "illumina.fastq.gz",
        output_dir=tmp_path
    )


def test_no_trim(tmp_path):
    """ --no-trim """
    run(
        "--no-trim --discard-untrimmed -a CCCTAGTTAAAC", "no-trim.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_bzip2(tmp_path):
    """test bzip2 support"""
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.bz2",
        output_dir=tmp_path
    )


@pytest.mark.skipif(lzma is None, reason="no lzma library")
def test_xz(tmp_path):
    """test xz support"""
    run(
        "-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.xz",
        output_dir=tmp_path
    )


def test_qualfile_only():
    retcode = execute_cli(["-sq", datapath("E3M.qual")])
    assert retcode == ReturnCode.ERROR


def test_no_args():
    with redirect_stderr():
        assert execute_cli() != 0


def test_anchored_no_indels(tmp_path):
    """anchored 5' adapter, mismatches only (no indels)"""
    run(
        "-g ^TTAGACATAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
        output_dir=tmp_path,
    )


def test_anchored_no_indels_wildcard_read(tmp_path):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the read
    count as matches.
    """
    run(
        "-g ^TTAGACATAT --match-read-wildcards --no-indels -e 0.1",
        "anchored_no_indels_wildcard.fasta",
        "anchored_no_indels.fasta",
        output_dir=tmp_path,
    )


def test_anchored_no_indels_wildcard_adapt(tmp_path):
    """anchored 5' adapter, mismatches only (no indels), but wildcards in the adapter
    count as matches.
    """
    run(
        "-g ^TTAGACANAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
        output_dir=tmp_path,
    )


def test_unconditional_cut_front(tmp_path):
    run(
        "-u 5", "unconditional-front.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_unconditional_cut_back(tmp_path):
    run(
        "-u -5", "unconditional-back.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_unconditional_cut_both(tmp_path):
    run(
        "-u -5 -u 5", "unconditional-both.fastq", "small.fastq",
        output_dir=tmp_path
    )


def test_untrimmed_output(tmp_path):
    tmp = tmp_path / "untrimmed.tmp.fastq"
    run(
        ["-a", "TTAGACATATCTCCGTCG", "--untrimmed-output", tmp],
        "small.trimmed.fastq",
        "small.fastq",
        output_dir=tmp_path
    )
    assert_files_equal(cutpath("small.untrimmed.fastq"), tmp)


def test_adapter_file(tmp_path):
    run(
        "-a " + dataurl("adapter.fasta"), "illumina.fastq", "illumina.fastq.gz",
        output_dir=tmp_path
    )


def test_adapter_file_5p_anchored(tmp_path):
    run(
        "-N -g " + dataurl("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
        output_dir=tmp_path,
    )


def test_adapter_file_3p_anchored(tmp_path):
    run(
        "-N -a " + dataurl("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
        output_dir=tmp_path,
    )


def test_adapter_file_5p_anchored_no_indels(tmp_path):
    run(
        "-N --no-indels -g " + dataurl("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
        output_dir=tmp_path,
    )


def test_adapter_file_3p_anchored_no_indels(tmp_path):
    run(
        "-N --no-indels -a " + dataurl("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
        output_dir=tmp_path,
    )


def test_demultiplex(tmp_path):
    multiout = str(tmp_path / "tmp-demulti.{name}.fasta")
    params = [
        "-a",
        "first=AATTTCAGGAATT",
        "-a",
        "second=GTTCTCTAGTTCT",
        "-o",
        multiout,
        "-se",
        datapathstr("twoadapters.fasta"),
    ]
    result = TrimCommandConsole.execute(params)
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result[0] == 0
    assert_files_equal(
        cutpath("twoadapters.first.fasta"), multiout.format(name="first")
    )
    assert_files_equal(
        cutpath("twoadapters.second.fasta"), multiout.format(name="second")
    )
    assert_files_equal(
        cutpath("twoadapters.unknown.fasta"), multiout.format(name="unknown")
    )


def test_max_n(tmp_path):
    run("--max-n 0", "maxn0.fasta", "maxn.fasta", output_dir=tmp_path)
    run("--max-n 1", "maxn1.fasta", "maxn.fasta", output_dir=tmp_path)
    run("--max-n 2", "maxn2.fasta", "maxn.fasta", output_dir=tmp_path)
    run("--max-n 0.2", "maxn0.2.fasta", "maxn.fasta", output_dir=tmp_path)
    run("--max-n 0.4", "maxn0.4.fasta", "maxn.fasta", output_dir=tmp_path)


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


def test_twocolor(tmp_path):
    run(
        "--twocolor-trim 22", "nextseq.fastq", "nextseq.fastq",
        output_dir=tmp_path
    )


def test_linked(tmp_path):
    run(
        "-a AAAAAAAAAA...TTTTTTTTTT", "linked.fasta", "linked.fasta",
        output_dir=tmp_path
    )


def test_fasta(tmp_path):
    run(
        "-a TTAGACATATCTCCGTCG", "small.fasta", "small.fastq",
        output_dir=tmp_path
    )


def test_custom_bisulfite_1(tmp_path):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,1", "small.fastq",
        "small.fastq", output_dir=tmp_path
    )


def test_custom_bisulfite_2(tmp_path):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 15,15,1,1",
        "small_mincut1.fastq",
        "small.fastq",
        output_dir=tmp_path,
    )


def test_custom_bisulfite_3(tmp_path):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,0",
        "small_mincut2.fastq",
        "small.fastq",
        output_dir=tmp_path,
    )


def test_custom_bisulfite_4(tmp_path):
    run(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,0,0",
        "small_mincut3.fastq",
        "small.fastq",
        output_dir=tmp_path,
    )


@pytest.mark.skipif(
    no_internet("https://ncbi.nlm.nih.gov") or no_import("srastream"),
    reason="No internet connection or srastream not importable",
)
def test_sra(tmp_path):
    run(
        "-b CTGGAGTTCAGACGTGTGCTCT --max-reads 100",
        "SRR2040662_trimmed.fq",
        sra_accn="SRR2040662",
        output_dir=tmp_path,
    )
