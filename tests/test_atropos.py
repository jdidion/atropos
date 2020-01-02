# TODO:
#  test with the --output option
#  test reading from standard input
import pytest
from xphyle import open_

from atropos.console import execute_cli
from atropos.utils import LOGGING_CONFIG, ReturnCode, no_import

from .utils import no_internet

try:
    import lzma
except ImportError:
    lzma = None


def test_example(run_trimmer):
    run_trimmer("-N -b ADAPTER", "example.fa", "example.fa")


def test_example_stdout(run_trimmer):
    run_trimmer("-N -b ADAPTER", "example.fa", "example.fa", stdout=True)


def test_small(run_trimmer):
    run_trimmer("-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq")


def test_empty(run_trimmer):
    run_trimmer("-a TTAGACATATCTCCGTCG", "empty.fastq", "empty.fastq")


def test_newlines(run_trimmer):
    run_trimmer("-e 0.12 -b TTAGACATATCTCCGTCG", "dos.fastq", "dos.fastq")


def test_lowercase(run_trimmer):
    run_trimmer("-b ttagacatatctccgtcg", "lowercase.fastq", "small.fastq")


def test_rest(run_trimmer, tmp_path):
    run_trimmer(
        ["-b", "ADAPTER", "-N", "-r", "rest.txt"],
        "rest.fa",
        "rest.fa",
        expected_other=["rest.txt"],
    )


def test_restfront(run_trimmer, tmp_path):
    run_trimmer(
        ["-g", "ADAPTER", "-N", "-r", "restfront.txt"],
        "restfront.fa",
        "rest.fa",
        expected_other=["restfront.txt"],
    )


def test_discard_trimmed(run_trimmer):
    run_trimmer(
        "-b TTAGACATATCTCCGTCG --discard-trimmed", "discard.fastq", "small.fastq"
    )


def test_discard_untrimmed(run_trimmer):
    run_trimmer(
        "-b CAAGAT --discard-untrimmed", "discard-untrimmed.fastq", "small.fastq"
    )


def test_plus(run_trimmer):
    run_trimmer("-e 0.12 -b TTAGACATATCTCCGTCG", "plus.fastq", "plus.fastq")


def test_extensiontxtgz(run_trimmer):
    run_trimmer("-b TTAGACATATCTCCGTCG", "s_1_sequence.txt", "s_1_sequence.txt.gz")


def test_format(run_trimmer):
    run_trimmer("-f fastq -b TTAGACATATCTCCGTCG", "small.fastq", "small.myownextension")


def test_minimum_length(run_trimmer):
    run_trimmer("-c -m 5 -a 330201030313112312", "minlen.fa", "lengths.fa")


def test_too_short(run_trimmer):
    run_trimmer(
        "-c -m 5 -a 330201030313112312 --too-short-output tooshort.fa",
        "minlen.fa",
        "lengths.fa",
        expected_other=["tooshort.fa"],
    )


def test_too_short_no_primer(run_trimmer):
    run_trimmer(
        "-c -m 5 -a 330201030313112312 --trim-primer "
        "--too-short-output tooshort.noprimer.fa",
        "minlen.noprimer.fa",
        "lengths.fa",
        expected_other=["tooshort.noprimer.fa"],
    )


def test_maximum_length(run_trimmer):
    run_trimmer("-c -M 5 -a 330201030313112312", "maxlen.fa", "lengths.fa")


def test_too_long(run_trimmer):
    run_trimmer(
        "-c -M 5 --too-long-output toolong.fa -a 330201030313112312",
        "maxlen.fa",
        "lengths.fa",
        expected_other=["toolong.fa"],
    )


def test_length_tag(run_trimmer):
    run_trimmer(
        "-n 3 -e 0.1 --length-tag length= "
        "-b TGAGACACGCAACAGGGGAAAGGCAAGGCACACAGGGGATAGG "
        "-b TCCATCTCATCCCTGCGTGTCCCATCTGTTCCCTCCCTGTCTCA",
        "454.fa",
        "454.fa",
    )


def test_overlap_a(run_trimmer):
    run_trimmer("-O 10 -a 330201030313112312 -e 0.0 -N", "overlapa.fa", "overlapa.fa")


def test_overlap_b(run_trimmer):
    run_trimmer("-O 10 -b TTAGACATATCTCCGTCG -N", "overlapb.fa", "overlapb.fa")


def test_qualtrim(run_trimmer):
    run_trimmer("-q 10 -a XXXXXX", "lowqual.fastq", "lowqual.fastq")


def test_qualbase(run_trimmer):
    run_trimmer(
        "-q 10 --quality-base 64 -a XXXXXX", "illumina64.fastq", "illumina64.fastq"
    )


def test_quality_trim_only(run_trimmer):
    run_trimmer("-q 10 --quality-base 64", "illumina64.fastq", "illumina64.fastq")


def test_twoadapters(run_trimmer):
    run_trimmer(
        "-a AATTTCAGGAATT -a GTTCTCTAGTTCT", "twoadapters.fasta", "twoadapters.fasta"
    )


def test_polya(run_trimmer):
    run_trimmer(
        "-m 24 -O 10 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "polya.fasta",
        "polya.fasta",
    )


def test_polya_brace_notation(run_trimmer):
    run_trimmer("-m 24 -O 10 -a A{35}", "polya.fasta", "polya.fasta")


def test_mask_adapter(run_trimmer):
    run_trimmer(
        "-b CAAG -n 3 --mask-adapter", "anywhere_repeat.fastq", "anywhere_repeat.fastq"
    )


def test_gz_multiblock(run_trimmer):
    run_trimmer("-b TTAGACATATCTCCGTCG", "small.fastq", "multiblock.fastq.gz")


def test_suffix(run_trimmer):
    run_trimmer(
        "-c -e 0.12 -a 1=330201030313112312 -y _my_suffix_{name} --strip-f3",
        "suffix.fastq",
        "solid.csfasta",
        qualfile="solid.qual",
    )


def test_read_wildcard(run_trimmer):
    run_trimmer("--match-read-wildcards -b ACGTACGT", "wildcard.fa", "wildcard.fa")


def test_adapter_wildcard(run_trimmer, tmp_path):
    for adapter_type, expected in (
        ("-a", "wildcard_adapter.fa"),
        ("-b", "wildcard_adapter_anywhere.fa"),
    ):
        wildcardtmp = tmp_path / "wildcardtmp.txt"

        run_trimmer(
            ("--wildcard-file", wildcardtmp, adapter_type, "ACGTNNNACGT"),
            expected,
            "wildcard_adapter.fa",
        )

        with open_(wildcardtmp) as wct:
            lines = [line.strip() for line in wct.readlines()]

        assert lines == ["AAA 1", "GGG 2", "CCC 3b", "TTT 4b"]


def test_wildcard_n(run_trimmer):
    run_trimmer(
        "-e 0 -a GGGGGGG --match-read-wildcards", "wildcardN.fa", "wildcardN.fa"
    )


def test_illumina_adapter_wildcard(run_trimmer):
    run_trimmer(
        "-a VCCGAMCYUCKHRKDCUBBCNUWNSGHCGU", "illumina.fastq", "illumina.fastq.gz"
    )


def test_adapter_front(run_trimmer):
    run_trimmer("--front ADAPTER -N", "examplefront.fa", "example.fa")


def test_literal_n(run_trimmer):
    run_trimmer("-N -e 0.2 -a NNNNNNNNNNNNNN", "trimN3.fasta", "trimN3.fasta")


def test_literal_n2(run_trimmer):
    run_trimmer("-N -O 1 -g NNNNNNNNNNNNNN", "trimN5.fasta", "trimN5.fasta")


def test_literal_n_brace_notation(run_trimmer):
    run_trimmer("-N -e 0.2 -a N{14}", "trimN3.fasta", "trimN3.fasta")


def test_literal_n2_brace_notation(run_trimmer):
    run_trimmer("-N -O 1 -g N{14}", "trimN5.fasta", "trimN5.fasta")


def test_anchored_front(run_trimmer):
    run_trimmer("-g ^FRONTADAPT -N", "anchored.fasta", "anchored.fasta")


def test_anchored_front_ellipsis_notation(run_trimmer):
    run_trimmer("-a FRONTADAPT... -N", "anchored.fasta", "anchored.fasta")


def test_anchored_back(run_trimmer):
    run_trimmer("-a BACKADAPTER$ -N", "anchored-back.fasta", "anchored-back.fasta")


def test_anchored_back_no_indels(run_trimmer):
    run_trimmer(
        "-a BACKADAPTER$ -N --no-indels", "anchored-back.fasta", "anchored-back.fasta"
    )


def test_no_indels(run_trimmer):
    run_trimmer(
        "-a TTAGACATAT -g GAGATTGCCA --no-indels", "no_indels.fasta", "no_indels.fasta"
    )


def test_issue_46(run_trimmer, tmp_path):
    wildcardtmp = tmp_path / "wildcardtmp.txt"
    run_trimmer(
        f"--anywhere=AACGTN --wildcard-file={wildcardtmp}",
        "issue46.fasta",
        "issue46.fasta",
    )


def test_strip_suffix(run_trimmer):
    run_trimmer("--strip-suffix _sequence -a XXXXXXX", "stripped.fasta", "simple.fasta")


def test_info_file(run_trimmer):
    run_trimmer(
        [
            "--info-file",
            "illumina.info.txt",
            "-a",
            "adapt=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        ],
        "illumina.fastq",
        "illumina.fastq.gz",
        expected_other=["illumina.info.txt"],
    )


def test_info_file_times(run_trimmer):
    run_trimmer(
        [
            "--info-file",
            "illumina5.info.txt",
            "--times",
            "2",
            "-a",
            "adapt=GCCGAACTTCTTA",
            "-a",
            "adapt2=GACTGCCTTAAGGACGT",
        ],
        "illumina5.fastq",
        "illumina5.fastq",
        expected_other=["illumina5.info.txt"],
    )


def test_info_file_fasta(run_trimmer, tmp_path):
    infotmp = tmp_path / "infotmp.txt"
    # Just make sure that it runs
    run_trimmer(
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
    )


def test_named_adapter(run_trimmer):
    run_trimmer(
        "-a MY_ADAPTER=GCCGAACTTCTTAGACTGCCTTAAGGACGT",
        "illumina.fastq",
        "illumina.fastq.gz",
    )


def test_adapter_with_u(run_trimmer):
    run_trimmer(
        "-a GCCGAACUUCUUAGACUGCCUUAAGGACGU", "illumina.fastq", "illumina.fastq.gz"
    )


def test_no_trim(run_trimmer):
    run_trimmer(
        "--no-trim --discard-untrimmed -a CCCTAGTTAAAC", "no-trim.fastq", "small.fastq"
    )


def test_bzip2(run_trimmer):
    run_trimmer("-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.bz2")


@pytest.mark.skipif(lzma is None, reason="no lzma library")
def test_xz(run_trimmer):
    run_trimmer("-b TTAGACATATCTCCGTCG", "small.fastq", "small.fastq.xz")


def test_qualfile_only(input_data):
    retcode = execute_cli(["-sq", input_data("E3M.qual")])
    assert retcode == ReturnCode.ERROR


def test_no_args():
    assert execute_cli() != ReturnCode.SUCCESS


def test_anchored_no_indels(run_trimmer):
    run_trimmer(
        "-g ^TTAGACATAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
    )


def test_anchored_no_indels_wildcard_read(run_trimmer):
    run_trimmer(
        "-g ^TTAGACATAT --match-read-wildcards --no-indels -e 0.1",
        "anchored_no_indels_wildcard.fasta",
        "anchored_no_indels.fasta",
    )


def test_anchored_no_indels_wildcard_adapt(run_trimmer):
    run_trimmer(
        "-g ^TTAGACANAT --no-indels -e 0.1",
        "anchored_no_indels.fasta",
        "anchored_no_indels.fasta",
    )


def test_unconditional_cut_front(run_trimmer):
    run_trimmer("-u 5", "unconditional-front.fastq", "small.fastq")


def test_unconditional_cut_back(run_trimmer):
    run_trimmer("-u -5", "unconditional-back.fastq", "small.fastq")


def test_unconditional_cut_both(run_trimmer):
    run_trimmer("-u -5 -u 5", "unconditional-both.fastq", "small.fastq")


def test_untrimmed_output(run_trimmer):
    run_trimmer(
        ["-a", "TTAGACATATCTCCGTCG", "--untrimmed-output", "small.untrimmed.fastq"],
        "small.trimmed.fastq",
        "small.fastq",
        expected_other=["small.untrimmed.fastq"],
    )


def test_adapter_file(run_trimmer, input_data):
    run_trimmer(
        "-a " + input_data.url("adapter.fasta"), "illumina.fastq", "illumina.fastq.gz"
    )


def test_adapter_file_5p_anchored(run_trimmer, input_data):
    run_trimmer(
        "-N -g " + input_data.url("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
    )


def test_adapter_file_3p_anchored(run_trimmer, input_data):
    run_trimmer(
        "-N -a " + input_data.url("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_adapter_file_5p_anchored_no_indels(run_trimmer, input_data):
    run_trimmer(
        "-N --no-indels -g " + input_data.url("prefix-adapter.fasta"),
        "anchored.fasta",
        "anchored.fasta",
    )


def test_adapter_file_3p_anchored_no_indels(run_trimmer, input_data):
    run_trimmer(
        "-N --no-indels -a " + input_data.url("suffix-adapter.fasta"),
        "anchored-back.fasta",
        "anchored-back.fasta",
    )


def test_demultiplex(run_trimmer):
    run_trimmer(
        ["-a", "first=AATTTCAGGAATT", "-a", "second=GTTCTCTAGTTCT"],
        "twoadapters.{name}.fasta",
        "twoadapters.fasta",
        expected_multi=("first", "second", "unknown"),
    )


def test_max_n(run_trimmer):
    run_trimmer("--max-n 0", "maxn0.fasta", "maxn.fasta")
    run_trimmer("--max-n 1", "maxn1.fasta", "maxn.fasta")
    run_trimmer("--max-n 2", "maxn2.fasta", "maxn.fasta")
    run_trimmer("--max-n 0.2", "maxn0.2.fasta", "maxn.fasta")
    run_trimmer("--max-n 0.4", "maxn0.4.fasta", "maxn.fasta")


def test_quiet_is_quiet(capsys, input_data, tmp_path):
    LOGGING_CONFIG.reset()
    tmp_out = tmp_path / "test.fq"
    execute_cli(
        [
            "-o",
            str(tmp_out),
            "--quiet",
            "-a",
            "XXXX",
            "-se",
            str(input_data("illumina.fastq.gz")),
        ]
    )
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_twocolor(run_trimmer):
    run_trimmer("--twocolor-trim 22", "nextseq.fastq", "nextseq.fastq")


def test_linked(run_trimmer):
    run_trimmer("-a AAAAAAAAAA...TTTTTTTTTT", "linked.fasta", "linked.fasta")


def test_fasta(run_trimmer):
    run_trimmer("-a TTAGACATATCTCCGTCG", "small.fasta", "small.fastq")


def test_custom_bisulfite_1(run_trimmer):
    run_trimmer(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,1", "small.fastq", "small.fastq"
    )


def test_custom_bisulfite_2(run_trimmer):
    run_trimmer(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 15,15,1,1",
        "small_mincut1.fastq",
        "small.fastq",
    )


def test_custom_bisulfite_3(run_trimmer):
    run_trimmer(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,1,0",
        "small_mincut2.fastq",
        "small.fastq",
    )


def test_custom_bisulfite_4(run_trimmer):
    run_trimmer(
        "-b TTAGACATATCTCCGTCG -q 0,0 --bisulfite 2,2,0,0",
        "small_mincut3.fastq",
        "small.fastq",
    )


@pytest.mark.skipif(
    (
        no_internet("https://ncbi.nlm.nih.gov")
        or no_import("ngstream")
        or no_import("ngs")
    ),
    reason="No internet connection or ngstream not importable",
)
def test_sra(run_trimmer):
    run_trimmer(
        "-b CTGGAGTTCAGACGTGTGCTCT --max-reads 100",
        "SRR2040662_trimmed.fq",
        sra_accn="SRR2040662",
    )
