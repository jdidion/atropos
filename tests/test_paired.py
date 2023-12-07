# coding: utf-8
import gzip
import os
import shutil
import pytest
from pytest import raises
from atropos.commands import execute_cli, get_command
from tests.utils import run, files_equal, datapath, cutpath, redirect_stderr

BACK_ALIGNERS = ("adapter", "insert")


def run_paired(
    params,
    in1,
    in2,
    expected1,
    expected2,
    tmp_path,
    aligners=("adapter",),
    callback=None,
    assert_files_equal=True,
    error_on_rc=True,
):
    if type(params) is str:
        params = params.split()
    for aligner in aligners:
        p1 = str(tmp_path / ("tmp1-" + expected1.format(aligner=aligner)))
        p2 = str(tmp_path / ("tmp2-" + expected2.format(aligner=aligner)))
        p = params.copy()
        p += ["--aligner", aligner, "-o", p1, "-p", p2]
        infiles = [datapath(i.format(aligner=aligner)) for i in (in1, in2)]
        for infile_args in zip(("-pe1", "-pe2"), infiles):
            p.extend(infile_args)
        command = get_command("trim")
        result = command.execute(p)
        assert isinstance(result, tuple)
        assert len(result) == 2
        if error_on_rc:
            err = (
                result[1]["exception"]
                if result[1] and "exception" in result[1]
                else None
            )
            # from traceback import format_exception
            if result[0] != 0:
                if err is None:
                    raise AssertionError("Return code {} != 0".format(result[0]))

                else:
                    raise AssertionError(
                        "Return code {} != 0".format(result[0])
                    ) from err["details"][1]

        if assert_files_equal:
            assert files_equal(cutpath(expected1.format(aligner=aligner)), p1)
            assert files_equal(cutpath(expected2.format(aligner=aligner)), p2)
        if callback:
            callback(aligner, infiles, (p1, p2), result)


def run_interleaved(params, inpath, expected, tmp_path, aligners=("adapter",)):
    if type(params) is str:
        params = params.split()
    for aligner in aligners:
        tmp = str(tmp_path / expected.format(aligner=aligner))
        p = params.copy()
        p += [
            "--aligner",
            aligner,
            "-l",
            datapath(inpath.format(aligner=aligner)),
            "-L",
            tmp,
        ]
        command = get_command("trim")
        result = command.execute(p)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert files_equal(cutpath(expected.format(aligner=aligner)), tmp)


# def run_interleaved2(params, inpath, expected1, expected2, aligners=('adapter',)):
#     assert False  # unused function
#     if type(params) is str:
#         params = params.split()
#     for aligner in aligners:
#         with temporary_path('tmp1-' + expected1.format(aligner=aligner)) as p1:
#             with temporary_path('tmp2-' + expected2.format(aligner=aligner)) as p2:
#                 p = params.copy()
#                 p += ['--aligner', aligner, '--interleaved', '-o', p1, '-p', p2]
#                 p += [datapath(inpath.format(aligner=aligner))]
#                 assert execute_cli(p) == 0
#                 assert files_equal(cutpath(expected.format(aligner=aligner)), p1)
#                 assert files_equal(cutpath(expected.format(aligner=aligner)), p2)


def test_paired_separate(tmp_path):
    """test separate trimming of paired-end reads"""
    run(tmp_path, "-a TTAGACATAT", "paired-separate.1.fastq", "paired.1.fastq")
    run(tmp_path, "-a CAGTGGAGTA", "paired-separate.2.fastq", "paired.2.fastq")


def test_paired_end_legacy(tmp_path):
    """--paired-output, not using -A/-B/-G"""
    # The -m 14 filters out one read, which should then also be filtered out in
    # the second output file.
    run_paired(
        "-a TTAGACATAT -m 14",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired.m14.1.fastq",
        expected2="paired.m14.2.fastq",
        tmp_path=tmp_path,
    )


def test_untrimmed_paired_output(tmp_path):
    untrimmed1 = str(tmp_path / "tmp-untrimmed.1.fastq")
    untrimmed2 = str(tmp_path / "tmp-untrimmed.2.fastq")

    def callback(aligner, infiles, outfiles, result):
        assert files_equal(cutpath("paired-untrimmed.1.fastq"), untrimmed1)
        assert files_equal(cutpath("paired-untrimmed.2.fastq"), untrimmed2)

    run_paired(
        [
            "-a",
            "TTAGACATAT",
            "--untrimmed-output",
            untrimmed1,
            "--untrimmed-paired-output",
            untrimmed2,
        ],
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-trimmed.1.fastq",
        expected2="paired-trimmed.2.fastq",
        tmp_path=tmp_path,
        callback=callback,
    )


def test_explicit_format_with_paired(tmp_path):
    # Use --format=fastq with input files whose extension is .txt
    txt1 = str(tmp_path / "paired.1.txt")
    txt2 = str(tmp_path / "paired.2.txt")
    shutil.copyfile(datapath("paired.1.fastq"), txt1)
    shutil.copyfile(datapath("paired.2.fastq"), txt2)
    run_paired(
        "--format=fastq -a TTAGACATAT -m 14",
        in1=txt1,
        in2=txt2,
        expected1="paired.m14.1.fastq",
        expected2="paired.m14.2.fastq",
        tmp_path=tmp_path,
    )


def test_no_trimming_legacy():
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
            datapath("paired.1.fastq"),
            "-pe2",
            datapath("paired.2.fastq"),
        ]
    )


def test_no_trimming():
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
            datapath("paired.1.fastq"),
            "-pe2",
            datapath("paired.2.fastq"),
        ]
    )


def test_missing_file():
    with raises(SystemExit), redirect_stderr():
        execute_cli(
            ["-a", "XX", "--paired-output", "out.fastq", datapath("paired.1.fastq")]
        )


def test_first_too_short(tmp_path):
    trunc1 = str(tmp_path / "truncated.1.fastq")
    # Create a truncated file in which the last read is missing
    with open(datapath("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    with open(trunc1, "w") as f:
        f.writelines(lines)
    with raises(SystemExit), redirect_stderr():
        execute_cli(
            "-a XX --paired-output out.fastq".split()
            + [trunc1, datapath("paired.2.fastq")]
        )


def test_second_too_short(tmp_path):
    trunc2 = str(tmp_path / "truncated.2.fastq")
    # Create a truncated file in which the last read is missing
    with open(datapath("paired.2.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    with open(trunc2, "w") as f:
        f.writelines(lines)
    with raises(SystemExit), redirect_stderr():
        execute_cli(
            "-a XX --paired-output out.fastq".split()
            + [datapath("paired.1.fastq"), trunc2]
        )


def test_unmatched_read_names(tmp_path):
    swapped = str(tmp_path / "swapped.1.fastq")
    try:
        # Create a file in which reads 2 and are swapped
        with open(datapath("paired.1.fastq")) as f:
            lines = f.readlines()
            lines = lines[0:4] + lines[8:12] + lines[4:8] + lines[12:]
        with open(swapped, "w") as f:
            f.writelines(lines)
        with redirect_stderr():
            command = get_command("trim")
            result = command.execute(
                "-a XX -o out1.fastq -p out2.fastq".split()
                + ["-pe1", swapped, "-pe2", datapath("paired.2.fastq")]
            )
            assert isinstance(result, tuple)
            assert len(result) == 2
            assert result[0] != 0
    finally:
        os.remove("out1.fastq")
        os.remove("out2.fastq")


def test_legacy_minlength(tmp_path):
    """Ensure -m is not applied to second read in a pair in legacy mode"""
    run_paired(
        "-a XXX -m 27",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-m27.1.fastq",
        expected2="paired-m27.2.fastq",
        tmp_path=tmp_path,
    )


def test_paired_end(tmp_path):
    """single-pass paired-end with -m"""
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired_{aligner}.1.fastq",
        expected2="paired_{aligner}.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_paired_anchored_back_no_indels(tmp_path):
    run_paired(
        "-a BACKADAPTER$ -A BACKADAPTER$ -N --no-indels",
        in1="anchored-back.fasta",
        in2="anchored-back.fasta",
        expected1="anchored-back.fasta",
        expected2="anchored-back.fasta",
        tmp_path=tmp_path,
    )


def test_paired_end_qualtrim(tmp_path):
    """single-pass paired-end with -q and -m"""
    run_paired(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="pairedq.1.fastq",
        expected2="pairedq.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_paired_end_qualtrim_swapped(tmp_path):
    """single-pass paired-end with -q and -m, but files swapped"""
    run_paired(
        "-q 20 -a CAGTGGAGTA -A TTAGACATAT -m 14 --adapter-max-rmp 0.001",
        in1="paired.2.fastq",
        in2="paired.1.fastq",
        expected1="pairedq.2.fastq",
        expected2="pairedq.1.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_paired_end_cut(tmp_path):
    run_paired(
        "-u 3 -u -1 -U 4 -U -2",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="pairedu.1.fastq",
        expected2="pairedu.2.fastq",
        tmp_path=tmp_path,
    )


def test_paired_end_A_only(tmp_path):
    run_paired(
        "-A CAGTGGAGTA",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-onlyA.1.fastq",
        expected2="paired-onlyA.2.fastq",
        tmp_path=tmp_path,
    )


def test_paired_end_mask_adapter(tmp_path):
    """mask adapter with N (reads maintain the same length)"""
    run_paired(
        "-a CAAG -A TCGA -n 3 --mask-adapter",
        in1="back_repeat.1.fastq",
        in2="back_repeat.2.fastq",
        expected1="back_repeat.1.fastq",
        expected2="back_repeat.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_discard_untrimmed(tmp_path):
    # issue #146
    # the first adapter is a sequence cut out from the first read
    run_paired(
        "-a CTCCAGCTTAGACATATC -A XXXXXXXX --discard-untrimmed",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="empty.fastq",
        expected2="empty.fastq",
        tmp_path=tmp_path,
    )


def test_discard_trimmed(tmp_path):
    run_paired(
        "-A C -O 1 --discard-trimmed",  # applies everywhere
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="empty.fastq",
        expected2="empty.fastq",
        tmp_path=tmp_path,
    )


def test_interleaved(tmp_path):
    """single-pass interleaved paired-end with -q and -m"""
    run_interleaved(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath="interleaved.fastq",
        expected="interleaved.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_interleaved_no_paired_output(tmp_path):
    p1 = str(tmp_path / "temp-paired.1.fastq")
    p2 = str(tmp_path / "temp-paired.2.fastq")
    params = "-a XX --interleaved".split()
    with raises(SystemExit), redirect_stderr():
        params += ["-o", p1, "-p1", p2, "paired.1.fastq", "paired.2.fastq"]
        execute_cli(params)


# TODO
# def test_interleaved_input_paired_output():
#     """single-pass interleaved paired-end with -q and -m, paired output"""
#     run_interleaved2('-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
#         inpath='interleaved.fastq', expected1='pairedq1.fastq', expected2='pairedq2.fastq'
#     )
def test_pair_filter(tmp_path):
    run_paired(
        "--pair-filter=both -a TTAGACATAT -A GGAGTA -m 14",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-filterboth_{aligner}.1.fastq",
        expected2="paired-filterboth_{aligner}.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_too_short_paired_output(tmp_path):
    p1 = str(tmp_path / "temp-too-short.1.fastq")
    p2 = str(tmp_path / "temp-too-short.2.fastq")

    def callback(aligner, infiles, outfiles, result):
        assert files_equal(cutpath("paired-too-short.1.fastq"), p1)
        assert files_equal(cutpath("paired-too-short.2.fastq"), p2)

    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output "
        "{0} --too-short-paired-output {1}".format(p1, p2),
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired_{aligner}.1.fastq",
        expected2="paired_{aligner}.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
        callback=callback,
    )


def test_too_long_output(tmp_path):
    p1 = str(tmp_path / "temp-too-long.1.fastq")
    p2 = str(tmp_path / "temp-too-long.2.fastq")

    def callback(aligner, infiles, outfiles, result):
        assert files_equal(
            cutpath("paired_{aligner}.1.fastq".format(aligner=aligner)), p1
        )
        assert files_equal(
            cutpath("paired_{aligner}.2.fastq".format(aligner=aligner)), p2
        )

    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -M 14 --too-long-output "
        "{0} --too-long-paired-output {1}".format(p1, p2),
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-too-short.1.fastq",
        expected2="paired-too-short.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
        callback=callback,
    )


def test_too_short_output_paired_option_missing(tmp_path):
    p1 = str(tmp_path / "temp-too-short.1.fastq")
    with raises(SystemExit):
        run_paired(
            "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output " "{0}".format(p1),
            in1="paired.1.fastq",
            in2="paired.2.fastq",
            expected1="paired.1.fastq",
            expected2="paired.2.fastq",
            tmp_path=tmp_path,
            aligners=BACK_ALIGNERS,
        )


def test_custom_bisulfite_1(tmp_path):
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 -q 0 --bisulfite 2,2,1,1",
        in1="paired_bis_{aligner}.1.fastq",
        in2="paired_bis_{aligner}.2.fastq",
        expected1="paired_bis1_{aligner}.1.fastq",
        expected2="paired_bis1_{aligner}.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_custom_bisulfite_2(tmp_path):
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 10 -q 0 --bisulfite 20,20,1,1;0,0,0,0",
        in1="paired_bis_{aligner}.1.fastq",
        in2="paired_bis_{aligner}.2.fastq",
        expected1="paired_bis2_{aligner}.1.fastq",
        expected2="paired_bis2_{aligner}.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
    )


def test_no_insert_match(tmp_path):
    # with -O
    # Note: this fails if you set -e 0.3 because the higher error rate enables a
    # match of 7 bp, even though it has 2 errors. This illustrates why using
    # --adapter-max-rmp is better.
    # run_paired('-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 7 -m 25 -q 0 --trim-n',
    #     in1='insert.1.fastq', in2='insert.2.fastq',
    #     expected1='insert.1.fastq', expected2='insert.2.fastq',
    #     aligners=('insert',)
    # )
    # with --adapter-max-rmp
    run_paired(
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -e 0.3 --adapter-max-rmp 0.001 -m 25 -q 0 --trim-n",
        in1="insert.1.fastq",
        in2="insert.2.fastq",
        expected1="insert.1.fastq",
        expected2="insert.2.fastq",
        tmp_path=tmp_path,
        aligners=("insert",),
    )


def test_overwrite(tmp_path):
    run_paired(
        "-w 10,30,10",
        in1="lowq.fastq",
        in2="highq.fastq",
        expected1="lowq.fastq",
        expected2="highq.fastq",
        tmp_path=tmp_path,
    )


def test_no_writer_process(tmp_path):
    def check_multifile(aligner, infiles, outfiles, result):
        assert os.path.basename(outfiles[0]) == "tmp1-out.1.fastq"
        assert os.path.basename(outfiles[1]) == "tmp2-out.2.fastq"
        tmpdir = os.path.dirname(outfiles[0])
        assert tmpdir == os.path.dirname(outfiles[1])
        # If a worker doesn't get the chance to process any batches, its output file is
        # never created. For now we just make sure that at least one pair of files was created.
        assert any(
            [
                (
                    os.path.exists(os.path.join(tmpdir, "tmp1-out.1.0.fastq"))
                    and os.path.exists(os.path.join(tmpdir, "tmp2-out.2.0.fastq"))
                ),
                (
                    os.path.exists(os.path.join(tmpdir, "tmp1-out.1.1.fastq"))
                    and os.path.exists(os.path.join(tmpdir, "tmp2-out.2.1.fastq"))
                ),
                (
                    os.path.exists(os.path.join(tmpdir, "tmp1-out.1.2.fastq"))
                    and os.path.exists(os.path.join(tmpdir, "tmp2-out.2.2.fastq"))
                ),
            ]
        )
        # TODO: check contents

    run_paired(
        "--threads 3 --no-writer-process --batch-size 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        in1="big.1.fq",
        in2="big.2.fq",
        expected1="out.1.fastq",
        expected2="out.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
        assert_files_equal=False,
        callback=check_multifile,
    )


def test_summary(tmp_path):
    def check_summary(aligner, infiles, outfiles, result):
        summary = result[1]
        assert summary is not None
        assert isinstance(summary, dict)
        assert summary["command"] == "trim"
        assert summary["options"]["orig_args"] == [
            "--threads",
            "2",
            "-a",
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
            "-A",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
            "--aligner",
            aligner,
            "-o",
            outfiles[0],
            "-p",
            outfiles[1],
            "-pe1",
            infiles[0],
            "-pe2",
            infiles[1],
        ]
        assert summary["sample_id"] == "big"
        assert summary["mode"] == "parallel"
        assert summary["threads"] == 2
        # infiles are 100 125 bp PE reads
        expected_record_counts = {}
        expected_record_counts[0] = 100
        assert summary["record_counts"] == expected_record_counts
        expected_bp_counts = {}
        expected_bp_counts[0] = [12500, 12500]
        assert summary["bp_counts"] == expected_bp_counts
        assert "timing" in summary
        assert "wallclock" in summary["timing"]
        assert summary["timing"]["wallclock"] > 0
        assert "cpu" in summary["timing"]
        assert summary["timing"]["cpu"] > 0

    run_paired(
        "--threads 2 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        in1="big.1.fq",
        in2="big.2.fq",
        expected1="out.1.fastq",
        expected2="out.2.fastq",
        tmp_path=tmp_path,
        aligners=BACK_ALIGNERS,
        assert_files_equal=False,
        callback=check_summary,
    )


class IssueTests:
    def test_issue68(self, tmp_path):
        run_paired(
            "--error-rate 0.20 --insert-match-error-rate 0.30 --minimum-length 20 --aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
            in1="issue68.1.fq",
            in2="issue68.2.fq",
            expected1="issue68.1.fq",
            expected2="issue68.2.fq",
            tmp_path=tmp_path,
            aligners=["insert"],
        )

    @pytest.mark.timeout(10)
    def test_issue122(self, tmp_path):
        # test that the empty fastq.gz files are valid gzip files
        def callback(aligner, infiles, outfiles, result):
            for out in outfiles:
                with gzip.open(out) as z:
                    assert z.read() == ""

        run_paired(
            "--threads 2 --preserve-order --no-default-adapters -a TTAGACATAT -A CAGTGGAGTA",
            in1="empty.fastq",
            in2="empty.fastq",
            expected1="empty.fastq.gz",
            expected2="empty.fastq.gz",
            tmp_path=tmp_path,
        )
