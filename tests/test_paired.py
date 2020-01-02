import shutil

import pytest
from xphyle import open_

from atropos.console import execute_cli
from atropos.commands.trim import AlignerType
from atropos.commands.trim.console import TrimCommandConsole
from atropos.utils import ReturnCode


BACK_ALIGNERS = tuple(
    a.name.lower() for a in (AlignerType.ADAPTER, AlignerType.INSERT)
)


def test_paired_separate(run_trimmer):
    run_trimmer(
        "-a TTAGACATAT",
        "paired-separate.1.fastq",
        "paired.1.fastq",
    )
    run_trimmer(
        "-a CAGTGGAGTA",
        "paired-separate.2.fastq",
        "paired.2.fastq",
    )


def test_paired_end_legacy(run_trimmer):
    # The -m 14 filters out one read, which should then also be filtered out in
    # the second output file.
    run_trimmer(
        "-a TTAGACATAT -m 14",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired.m14.1.fastq",
        expected2="paired.m14.2.fastq",
    )


def test_untrimmed_paired_output(run_trimmer):
    run_trimmer(
        [
            "-a",
            "TTAGACATAT",
            "--untrimmed-output",
            "paired-untrimmed.1.fastq",
            "--untrimmed-paired-output",
            "paired-untrimmed.2.fastq",
        ],
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired-trimmed.1.fastq",
        expected2="paired-trimmed.2.fastq",
        expected_other=("paired-untrimmed.1.fastq", "paired-untrimmed.2.fastq")
    )


def test_explicit_format_with_paired(run_trimmer, input_data, tmp_path):
    # Use --format=fastq with input files whose extension is .txt
    txt1 = tmp_path / "paired.1.txt"
    txt2 = tmp_path / "paired.2.txt"
    shutil.copyfile(input_data("paired.1.fastq"), txt1)
    shutil.copyfile(input_data("paired.2.fastq"), txt2)
    run_trimmer(
        "--input-format=fastq -a TTAGACATAT -m 14",
        inpath1=txt1,
        inpath2=txt2,
        expected1="paired.m14.1.fastq",
        expected2="paired.m14.2.fastq",
    )


def test_no_trimming_legacy(input_data):
    # make sure that this doesn't divide by zero
    execute_cli(
        [
            "-a",
            "XXXXX",
            "-o",
            "/dev/null",
            "-p",
            "/dev/null",
            "-pe1",
            input_data("paired.1.fastq"),
            "-pe2",
            input_data("paired.2.fastq"),
        ]
    )


def test_no_trimming(input_data):
    # make sure that this doesn't divide by zero
    execute_cli(
        [
            "-a",
            "XXXXX",
            "-A",
            "XXXXX",
            "-o",
            "/dev/null",
            "-p",
            "/dev/null",
            "-pe1",
            input_data("paired.1.fastq"),
            "-pe2",
            input_data("paired.2.fastq"),
        ]
    )


def test_missing_file(input_data):
    assert ReturnCode.ERROR == execute_cli(
        ["-a", "XX", "--paired-output", "out.fastq", input_data("paired.1.fastq")]
    )


def test_first_too_short(run_trimmer, input_data, tmp_path):
    trunc1 = tmp_path / "truncated.1.fastq"
    out = tmp_path / "out.fastq"
    # Create a truncated file in which the last read is missing
    with open_(input_data("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    with open_(trunc1, "w") as f:
        f.writelines(lines)
    assert ReturnCode.ERROR == execute_cli([
        "-a", "XX", "--paired-output", out, trunc1, input_data("paired.2.fastq")
    ])


def test_second_too_short(run_trimmer, input_data, tmp_path):
    trunc2 = tmp_path / "truncated.2.fastq"
    out = tmp_path / "out.fastq"
    # Create a truncated file in which the last read is missing
    with open_(input_data("paired.2.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    with open_(trunc2, "w") as f:
        f.writelines(lines)
    assert ReturnCode.ERROR == execute_cli([
        "-a", "XX", "--paired-output", out, input_data("paired.1.fastq"), trunc2
    ])


def test_unmatched_read_names(run_trimmer, input_data, tmp_path):
    swapped = tmp_path / "swapped.1.fastq"
    out1 = tmp_path / "out1.fastq"
    out2 = tmp_path / "out2.fastq"

    # Create a file in which reads 2 and are swapped
    with open_(input_data("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[0:4] + lines[8:12] + lines[4:8] + lines[12:]
    with open_(swapped, "w") as f:
        f.writelines(lines)
    result = TrimCommandConsole.execute([
        "-a", "XX", "-o", str(out1), "-p", str(out2), "-pe1", str(swapped), "-pe2",
        str(input_data("paired.2.fastq"))
    ])
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result[0] != 0


def test_legacy_minlength(run_trimmer):
    """Ensure -m is not applied to second read in a pair in legacy mode"""
    run_trimmer(
        "-a XXX -m 27",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired-m27.1.fastq",
        expected2="paired-m27.2.fastq",
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_paired_end(run_trimmer, aligner):
    """single-pass paired-end with -m"""
    run_trimmer(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired_{aligner}.1.fastq",
        expected2="paired_{aligner}.2.fastq",
        aligner=aligner
    )


def test_paired_anchored_back_no_indels(run_trimmer):
    run_trimmer(
        "-a BACKADAPTER$ -A BACKADAPTER$ -N --no-indels",
        inpath1="anchored-back.fasta",
        inpath2="anchored-back.fasta",
        expected1="anchored-back.fasta",
        expected2="anchored-back.fasta",
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_paired_end_qualtrim(run_trimmer, aligner):
    """single-pass paired-end with -q and -m"""
    run_trimmer(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="pairedq.1.fastq",
        expected2="pairedq.2.fastq",
        aligner=aligner,
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_paired_end_qualtrim_swapped(run_trimmer, aligner):
    """single-pass paired-end with -q and -m, but files swapped"""
    run_trimmer(
        "-q 20 -a CAGTGGAGTA -A TTAGACATAT -m 14 --adapter-max-rmp 0.001",
        inpath1="paired.2.fastq",
        inpath2="paired.1.fastq",
        expected1="pairedq.2.fastq",
        expected2="pairedq.1.fastq",
        aligner=aligner,
    )


def test_paired_end_cut(run_trimmer):
    run_trimmer(
        "-u 3 -u -1 -U 4 -U -2",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="pairedu.1.fastq",
        expected2="pairedu.2.fastq",
    )


def test_paired_end_a_only(run_trimmer):
    run_trimmer(
        "-A CAGTGGAGTA",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired-onlyA.1.fastq",
        expected2="paired-onlyA.2.fastq",
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_paired_end_mask_adapter(run_trimmer, aligner):
    """mask adapter with N (reads maintain the same length)"""
    run_trimmer(
        "-a CAAG -A TCGA -n 3 --mask-adapter",
        inpath1="back_repeat.1.fastq",
        inpath2="back_repeat.2.fastq",
        expected1="back_repeat.1.fastq",
        expected2="back_repeat.2.fastq",
        aligner=aligner,
    )


def test_discard_untrimmed(run_trimmer):
    # issue #146
    # the first adapter is a sequence cut out from the first read
    run_trimmer(
        "-a CTCCAGCTTAGACATATC -A XXXXXXXX --discard-untrimmed",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="empty.fastq",
        expected2="empty.fastq",
    )


def test_discard_trimmed(run_trimmer):
    run_trimmer(
        "-A C -O 1 --discard-trimmed",  # applies everywhere
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="empty.fastq",
        expected2="empty.fastq",
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_interleaved(run_trimmer, aligner):
    """single-pass interleaved paired-end with -q and -m"""
    run_trimmer(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath1="interleaved.fastq",
        expected1="interleaved.fastq",
        interleaved_input=True,
        interleaved_output=True,
        aligner=aligner,
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_interleaved_stdout(run_trimmer, aligner):
    """single-pass interleaved paired-end with -q and -m"""
    run_trimmer(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath1="interleaved.fastq",
        expected1="interleaved.fastq",
        interleaved_input=True,
        interleaved_output=True,
        aligner=aligner,
        stdout=True,
    )


def test_interleaved_no_paired_output(tmp_path):
    p1 = tmp_path / "temp-paired.1.fastq"
    p2 = tmp_path / "temp-paired.2.fastq"
    params = "-a XX --interleaved".split()
    params += ["-o", p1, "-p1", p2, "paired.1.fastq", "paired.2.fastq"]
    assert execute_cli(params) == ReturnCode.ERROR


# TODO
# def test_interleaved_input_paired_output():
#     '''single-pass interleaved paired-end with -q and -m, paired output'''
#     run_interleaved2('-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
#         inpath1='interleaved.fastq', expected1='pairedq1.fastq',
#         expected2='pairedq2.fastq'
#     )
@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_pair_filter(run_trimmer, aligner):
    run_trimmer(
        "--pair-filter=both -a TTAGACATAT -A GGAGTA -m 14",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired-filterboth_{aligner}.1.fastq",
        expected2="paired-filterboth_{aligner}.2.fastq",
        aligner=aligner,
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_too_short_paired_output(run_trimmer, aligner, expected_data):
    run_trimmer(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output "
        "paired-too-short.1.fastq --too-short-paired-output paired-too-short.2.fastq",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired_{aligner}.1.fastq",
        expected2="paired_{aligner}.2.fastq",
        aligner=aligner,
        expected_other=("paired-too-short.1.fastq", "paired-too-short.2.fastq")
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_too_long_output(run_trimmer, aligner):
    run_trimmer(
        f"-a TTAGACATAT -A CAGTGGAGTA -M 14 --too-long-output "
        f"paired_{aligner}.1.fastq --too-long-paired-output paired_{aligner}.2.fastq",
        inpath1="paired.1.fastq",
        inpath2="paired.2.fastq",
        expected1="paired-too-short.1.fastq",
        expected2="paired-too-short.2.fastq",
        aligner=aligner,
        expected_other=(
            f"paired_{aligner}.1.fastq", f"paired_{aligner}.2.fastq"
        )
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_too_short_output_paired_option_missing(run_trimmer, aligner, tmp_path):
    with pytest.raises(SystemExit):
        p1 = tmp_path / "temp-too-short.1.fastq"
        run_trimmer(
            f"-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output {p1}",
            inpath1="paired.1.fastq",
            inpath2="paired.2.fastq",
            expected1="paired.1.fastq",
            expected2="paired.2.fastq",
            aligner=aligner,
        )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_custom_bisulfite_1(run_trimmer, aligner):
    run_trimmer(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 -q 0 --bisulfite 2,2,1,1",
        inpath1="paired_bis_{aligner}.1.fastq",
        inpath2="paired_bis_{aligner}.2.fastq",
        expected1="paired_bis1_{aligner}.1.fastq",
        expected2="paired_bis1_{aligner}.2.fastq",
        aligner=aligner,
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_custom_bisulfite_2(run_trimmer, aligner):
    run_trimmer(
        "-a TTAGACATAT -A CAGTGGAGTA -m 10 -q 0 --bisulfite 20,20,1,1;0,0,0,0",
        inpath1="paired_bis_{aligner}.1.fastq",
        inpath2="paired_bis_{aligner}.2.fastq",
        expected1="paired_bis2_{aligner}.1.fastq",
        expected2="paired_bis2_{aligner}.2.fastq",
        aligner=aligner,
    )


def test_no_insert_match(run_trimmer):
    # with -O
    # Note: this fails if you set -e 0.3 because the higher error rate enables a
    # match of 7 bp, even though it has 2 errors. This illustrates why using
    # --adapter-max-rmp is better.
    # run_trimmer('-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
    # -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 7 -m 25 -q 0
    # --trim-n',
    #     inpath1='insert.1.fastq', inpath2='insert.2.fastq',
    #     expected1='insert.1.fastq', expected2='insert.2.fastq',
    #     aligner='insert'
    # )
    # with --adapter-max-rmp
    run_trimmer(
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG "
        "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -e 0.3 "
        "--adapter-max-rmp 0.001 -m 25 -q 0 --trim-n",
        inpath1="insert.1.fastq",
        inpath2="insert.2.fastq",
        expected1="insert.1.fastq",
        expected2="insert.2.fastq",
        aligner="insert",
    )


def test_overwrite(run_trimmer):
    run_trimmer(
        "-w 10,30,10",
        inpath1="lowq.fastq",
        inpath2="highq.fastq",
        expected1="lowq.fastq",
        expected2="highq.fastq",
    )


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_no_writer_process(run_trimmer, aligner):
    # TODO: check contents
    output_dir, _, _ = run_trimmer(
        "--threads 3 --no-writer-process --batch-size 1 "
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG "
        "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        inpath1="big.1.fq",
        inpath2="big.2.fq",
        expected1="out.1.fastq",
        expected2="out.2.fastq",
        aligner=aligner,
        assert_output_equal=False,
    )
    # TODO: If the final worker doesn't get the chance to process any
    #  batches, the last output file is never created.
    for i in (1, 2):
        for j in (0, 1):
            assert (output_dir / f"tmp{i}-out.{i}.{j}.fastq").exists()


@pytest.mark.parametrize("aligner", BACK_ALIGNERS)
def test_summary(run_trimmer, aligner):
    output_dir, infiles, summary = run_trimmer(
        "--threads 2 "
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG "
        "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        inpath1="big.1.fq",
        inpath2="big.2.fq",
        expected1="out.1.fastq",
        expected2="out.2.fastq",
        aligner=aligner,
        assert_output_equal=False,
    )

    assert summary is not None
    assert isinstance(summary, dict)
    assert summary["command"] == "trim"
    assert summary["options"]["orig_args"] == (
        "--threads",
        "2",
        "-a",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        "-A",
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        "--aligner",
        aligner,
        "-pe1",
        str(infiles[0]),
        "-pe2",
        str(infiles[1]),
        "-o",
        str(output_dir / "tmp1-out.1.fastq"),
        "-p",
        str(output_dir / "tmp2-out.2.fastq"),
    )
    assert summary["sample_id"] == "big"
    assert summary["mode"] == "parallel"
    assert summary["threads"] == 2
    # infiles are 100 125 bp PE reads
    expected_record_counts = {0: 100}
    assert summary["record_counts"] == expected_record_counts
    expected_bp_counts = {0: [12500, 12500]}
    assert summary["bp_counts"] == expected_bp_counts
    assert "timing" in summary
    assert "wallclock" in summary["timing"]
    assert summary["timing"]["wallclock"] > 0
    assert "cpu" in summary["timing"]
    assert summary["timing"]["cpu"] > 0


def test_sam(run_trimmer):
    run_trimmer(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --output-format sam",
        "paired_insert.sam",
        "paired.sam",
        interleaved_input=True,
        interleaved_output=True,
        aligner="insert",
    )


# def test_long_reads():
#    run_trimmer(''
#        inpath1='long.1.fq', inpath2='long.2.fq',
#        expected1=, expected2=,
#        aligner=BACK_ALIGNERS)


def test_issue68(run_trimmer):
    run_trimmer(
        "--error-rate 0.20 --insert-match-error-rate 0.30 --minimum-length 20 "
        "--aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
        "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        inpath1="issue68.1.fq",
        inpath2="issue68.2.fq",
        expected1="issue68.1.fq",
        expected2="issue68.2.fq",
        aligner="insert",
    )
