import os
import shutil
from typing import Callable, Iterable, Union

from pytest import raises
from xphyle import open_

from atropos.console import execute_cli
from atropos.commands.trim.console import TrimCommandConsole
from atropos.utils import ReturnCode

from .utils import (
    run,
    assert_files_equal,
    datapath,
    cutpath,
    redirect_stderr,
    intercept_stdout,
)

BACK_ALIGNERS = ("adapter", "insert")


def run_paired(
    params: Union[str, list],
    in1: str,
    in2: str,
    expected1: str,
    expected2: str,
    aligners: Iterable[str] = ("adapter",),
    callback: Callable = None,
    files_equal: bool = True,
    error_on_rc: bool = True,
    output_dir=None,
):
    if type(params) is str:
        params = params.split()

    for aligner in aligners:
        p1 = output_dir / f"tmp1-{expected1.format(aligner=aligner)}"
        p2 = output_dir / f"tmp2-{expected2.format(aligner=aligner)}"
        p = params.copy()
        p += ["--aligner", aligner, "-o", str(p1), "-p", str(p2)]

        infiles = []

        for i, infile in enumerate((in1, in2), 1):
            inpath = datapath(infile.format(aligner=aligner))
            infiles.append(inpath)
            p.append(f"-pe{i}")
            p.append(str(inpath))

        result = TrimCommandConsole.execute(p)

        assert isinstance(result, tuple)
        assert len(result) == 2

        if error_on_rc:
            err = (
                result[1]["exception"]
                if result[1] and "exception" in result[1]
                else None
            )
            if result[0] != 0:
                if err is None:
                    raise AssertionError("Return code {} != 0".format(result[0]))
                else:
                    raise AssertionError(
                        "Return code {} != 0".format(result[0])
                    ) from err["details"][1]

        if files_equal:
            assert_files_equal(cutpath(expected1.format(aligner=aligner)), p1)
            assert_files_equal(cutpath(expected2.format(aligner=aligner)), p2)

        if callback:
            callback(aligner, infiles, (p1, p2), result)


def run_interleaved(
    params: Union[str, list],
    inpath: str,
    expected: str,
    aligners: Iterable[str] = ("adapter",),
    stdout: bool = False,
    error_on_rc: bool = True,
    output_dir=None,
):
    if type(params) is str:
        params = params.split()

    for aligner in aligners:
        tmp = output_dir / expected.format(aligner=aligner)
        p = params.copy()
        p += ["--aligner", aligner, "-l", str(datapath(inpath.format(aligner=aligner)))]

        if stdout:
            # Output is going to stdout, so we need to redirect it to the
            # temp file
            with intercept_stdout() as stdout:
                # print(params)
                result = TrimCommandConsole.execute(p)
                with open(tmp, "wt") as out:
                    out.write(stdout.getvalue())
        else:
            p.extend(["-L", str(tmp)])
            result = TrimCommandConsole.execute(p)

        assert isinstance(result, tuple)
        assert len(result) == 2

        if error_on_rc:
            err = (
                result[1]["exception"]
                if result[1] and "exception" in result[1]
                else None
            )
            if result[0] != 0:
                raise AssertionError("Return code {} != 0".format(result[0])) from err[
                    "details"
                ][1]

        assert_files_equal(cutpath(expected.format(aligner=aligner)), tmp)


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
    run(
        "-a TTAGACATAT",
        "paired-separate.1.fastq",
        "paired.1.fastq",
        output_dir=tmp_path,
    )
    run(
        "-a CAGTGGAGTA",
        "paired-separate.2.fastq",
        "paired.2.fastq",
        output_dir=tmp_path,
    )


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
        output_dir=tmp_path,
    )


def test_untrimmed_paired_output(tmp_path):
    untrimmed1 = tmp_path / "tmp-untrimmed.1.fastq"
    untrimmed2 = tmp_path / "tmp-untrimmed.2.fastq"

    def callback(aligner, infiles, outfiles, result):
        assert_files_equal(cutpath("paired-untrimmed.1.fastq"), untrimmed1)
        assert_files_equal(cutpath("paired-untrimmed.2.fastq"), untrimmed2)

    run_paired(
        [
            "-a",
            "TTAGACATAT",
            "--untrimmed-output",
            str(untrimmed1),
            "--untrimmed-paired-output",
            str(untrimmed2),
        ],
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-trimmed.1.fastq",
        expected2="paired-trimmed.2.fastq",
        callback=callback,
        output_dir=tmp_path,
    )


def test_explicit_format_with_paired(tmp_path):
    # Use --format=fastq with input files whose extension is .txt
    txt1 = tmp_path / "paired.1.txt"
    txt2 = tmp_path / "paired.2.txt"
    shutil.copyfile(datapath("paired.1.fastq"), txt1)
    shutil.copyfile(datapath("paired.2.fastq"), txt2)
    run_paired(
        "--input-format=fastq -a TTAGACATAT -m 14",
        in1=str(txt1),
        in2=str(txt2),
        expected1="paired.m14.1.fastq",
        expected2="paired.m14.2.fastq",
        output_dir=tmp_path,
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
    with redirect_stderr():
        assert ReturnCode.ERROR == execute_cli(
            ["-a", "XX", "--paired-output", "out.fastq", datapath("paired.1.fastq")]
        )


def test_first_too_short(tmp_path):
    trunc1 = tmp_path / "truncated.1.fastq"
    out = tmp_path / "out.fastq"
    # Create a truncated file in which the last read is missing
    with open_(datapath("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    with open_(trunc1, "w") as f:
        f.writelines(lines)
    with redirect_stderr():
        assert ReturnCode.ERROR == execute_cli([
            "-a", "XX", "--paired-output", str(out),
            str(trunc1), str(datapath("paired.2.fastq"))
        ])


def test_second_too_short(tmp_path):
    trunc2 = tmp_path / "truncated.2.fastq"
    out = tmp_path / "out.fastq"
    # Create a truncated file in which the last read is missing
    with open_(datapath("paired.2.fastq")) as f:
        lines = f.readlines()
        lines = lines[:-4]
    with open_(trunc2, "w") as f:
        f.writelines(lines)
    with redirect_stderr():
        assert ReturnCode.ERROR == execute_cli([
            "-a", "XX", "--paired-output", str(out),
            str(datapath("paired.1.fastq")), str(trunc2)
        ])


def test_unmatched_read_names(tmp_path):
    swapped = tmp_path / "swapped.1.fastq"
    out1 = tmp_path / "out1.fastq"
    out2 = tmp_path / "out2.fastq"

    # Create a file in which reads 2 and are swapped
    with open_(datapath("paired.1.fastq")) as f:
        lines = f.readlines()
        lines = lines[0:4] + lines[8:12] + lines[4:8] + lines[12:]
    with open_(swapped, "w") as f:
        f.writelines(lines)
    with redirect_stderr():
        result = TrimCommandConsole.execute([
            "-a", "XX", "-o", str(out1), "-p", str(out2),
            "-pe1", str(swapped), "-pe2", str(datapath("paired.2.fastq"))
        ])
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert result[0] != 0


def test_legacy_minlength(tmp_path):
    """Ensure -m is not applied to second read in a pair in legacy mode"""
    run_paired(
        "-a XXX -m 27",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-m27.1.fastq",
        expected2="paired-m27.2.fastq",
        output_dir=tmp_path,
    )


def test_paired_end(tmp_path):
    """single-pass paired-end with -m"""
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired_{aligner}.1.fastq",
        expected2="paired_{aligner}.2.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
    )


def test_paired_anchored_back_no_indels(tmp_path):
    run_paired(
        "-a BACKADAPTER$ -A BACKADAPTER$ -N --no-indels",
        in1="anchored-back.fasta",
        in2="anchored-back.fasta",
        expected1="anchored-back.fasta",
        expected2="anchored-back.fasta",
        output_dir=tmp_path,
    )


def test_paired_end_qualtrim(tmp_path):
    """single-pass paired-end with -q and -m"""
    run_paired(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="pairedq.1.fastq",
        expected2="pairedq.2.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
    )


def test_paired_end_qualtrim_swapped(tmp_path):
    """single-pass paired-end with -q and -m, but files swapped"""
    run_paired(
        "-q 20 -a CAGTGGAGTA -A TTAGACATAT -m 14 --adapter-max-rmp 0.001",
        in1="paired.2.fastq",
        in2="paired.1.fastq",
        expected1="pairedq.2.fastq",
        expected2="pairedq.1.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
    )


def test_paired_end_cut(tmp_path):
    run_paired(
        "-u 3 -u -1 -U 4 -U -2",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="pairedu.1.fastq",
        expected2="pairedu.2.fastq",
        output_dir=tmp_path,
    )


def test_paired_end_a_only(tmp_path):
    run_paired(
        "-A CAGTGGAGTA",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-onlyA.1.fastq",
        expected2="paired-onlyA.2.fastq",
        output_dir=tmp_path,
    )


def test_paired_end_mask_adapter(tmp_path):
    """mask adapter with N (reads maintain the same length)"""
    run_paired(
        "-a CAAG -A TCGA -n 3 --mask-adapter",
        in1="back_repeat.1.fastq",
        in2="back_repeat.2.fastq",
        expected1="back_repeat.1.fastq",
        expected2="back_repeat.2.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
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
        output_dir=tmp_path,
    )


def test_discard_trimmed(tmp_path):
    run_paired(
        "-A C -O 1 --discard-trimmed",  # applies everywhere
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="empty.fastq",
        expected2="empty.fastq",
        output_dir=tmp_path,
    )


def test_interleaved(tmp_path):
    """single-pass interleaved paired-end with -q and -m"""
    run_interleaved(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath="interleaved.fastq",
        expected="interleaved.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
    )


def test_interleaved_stdout(tmp_path):
    """single-pass interleaved paired-end with -q and -m"""
    run_interleaved(
        "-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90",
        inpath="interleaved.fastq",
        expected="interleaved.fastq",
        aligners=BACK_ALIGNERS,
        stdout=True,
        output_dir=tmp_path
    )


def test_interleaved_no_paired_output(tmp_path):
    p1 = tmp_path / "temp-paired.1.fastq"
    p2 = tmp_path / "temp-paired.2.fastq"
    params = "-a XX --interleaved".split()
    with redirect_stderr():
        params += ["-o", p1, "-p1", p2, "paired.1.fastq", "paired.2.fastq"]
        assert execute_cli(params) == ReturnCode.ERROR


# TODO
# def test_interleaved_input_paired_output():
#     '''single-pass interleaved paired-end with -q and -m, paired output'''
#     run_interleaved2('-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
#         inpath='interleaved.fastq', expected1='pairedq1.fastq',
#         expected2='pairedq2.fastq'
#     )
def test_pair_filter(tmp_path):
    run_paired(
        "--pair-filter=both -a TTAGACATAT -A GGAGTA -m 14",
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-filterboth_{aligner}.1.fastq",
        expected2="paired-filterboth_{aligner}.2.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
    )


def test_too_short_paired_output(tmp_path):
    p1 = tmp_path / "temp-too-short.1.fastq"
    p2 = tmp_path / "temp-too-short.2.fastq"

    def callback(aligner, infiles, outfiles, result):
        assert_files_equal(cutpath("paired-too-short.1.fastq"), p1)
        assert_files_equal(cutpath("paired-too-short.2.fastq"), p2)

    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output "
        "{0} --too-short-paired-output {1}".format(p1, p2),
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired_{aligner}.1.fastq",
        expected2="paired_{aligner}.2.fastq",
        aligners=BACK_ALIGNERS,
        callback=callback,
        output_dir=tmp_path,
    )


def test_too_long_output(tmp_path):
    p1 = tmp_path / "temp-too-long.1.fastq"
    p2 = tmp_path /"temp-too-long.2.fastq"

    def callback(aligner, infiles, outfiles, result):
        assert_files_equal(
            cutpath("paired_{aligner}.1.fastq".format(aligner=aligner)), p1
        )
        assert_files_equal(
            cutpath("paired_{aligner}.2.fastq".format(aligner=aligner)), p2
        )

    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -M 14 --too-long-output "
        "{0} --too-long-paired-output {1}".format(p1, p2),
        in1="paired.1.fastq",
        in2="paired.2.fastq",
        expected1="paired-too-short.1.fastq",
        expected2="paired-too-short.2.fastq",
        aligners=BACK_ALIGNERS,
        callback=callback,
        output_dir=tmp_path,
    )


def test_too_short_output_paired_option_missing(tmp_path):
    with raises(SystemExit):
        p1 = tmp_path / "temp-too-short.1.fastq"
        run_paired(
            "-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output " "{0}".format(p1),
            in1="paired.1.fastq",
            in2="paired.2.fastq",
            expected1="paired.1.fastq",
            expected2="paired.2.fastq",
            aligners=BACK_ALIGNERS,
            output_dir=tmp_path,
        )


def test_custom_bisulfite_1(tmp_path):
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 -q 0 --bisulfite 2,2,1,1",
        in1="paired_bis_{aligner}.1.fastq",
        in2="paired_bis_{aligner}.2.fastq",
        expected1="paired_bis1_{aligner}.1.fastq",
        expected2="paired_bis1_{aligner}.2.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
    )


def test_custom_bisulfite_2(tmp_path):
    run_paired(
        "-a TTAGACATAT -A CAGTGGAGTA -m 10 -q 0 --bisulfite 20,20,1,1;0,0,0,0",
        in1="paired_bis_{aligner}.1.fastq",
        in2="paired_bis_{aligner}.2.fastq",
        expected1="paired_bis2_{aligner}.1.fastq",
        expected2="paired_bis2_{aligner}.2.fastq",
        aligners=BACK_ALIGNERS,
        output_dir=tmp_path,
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
        aligners=("insert",),
        output_dir=tmp_path,
    )


def test_overwrite(tmp_path):
    run_paired(
        "-w 10,30,10",
        in1="lowq.fastq",
        in2="highq.fastq",
        expected1="lowq.fastq",
        expected2="highq.fastq",
        output_dir=tmp_path,
    )


def test_no_writer_process(tmp_path):
    def check_multifile(aligner, infiles, outfiles, result):
        assert os.path.basename(outfiles[0]) == "tmp1-out.1.fastq"
        assert os.path.basename(outfiles[1]) == "tmp2-out.2.fastq"
        tmpdir = os.path.dirname(outfiles[0])
        assert tmpdir == os.path.dirname(outfiles[1])
        # TODO: If the final worker doesn't get the chance to process any
        # batches, the last output file is never created.
        assert os.path.exists(os.path.join(tmpdir, "tmp1-out.1.0.fastq"))
        assert os.path.exists(os.path.join(tmpdir, "tmp1-out.1.1.fastq"))
        # assert os.path.exists(os.path.join(tmpdir, 'tmp1-out.1.2.fastq'))
        assert os.path.exists(os.path.join(tmpdir, "tmp2-out.2.0.fastq"))
        assert os.path.exists(os.path.join(tmpdir, "tmp2-out.2.1.fastq"))

    # assert os.path.exists(os.path.join(tmpdir, 'tmp2-out.2.2.fastq'))
    # TODO: check contents
    run_paired(
        "--threads 3 --no-writer-process --batch-size 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        in1="big.1.fq",
        in2="big.2.fq",
        expected1="out.1.fastq",
        expected2="out.2.fastq",
        aligners=BACK_ALIGNERS,
        files_equal=False,
        callback=check_multifile,
        output_dir=tmp_path,
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
            str(outfiles[0]),
            "-p",
            str(outfiles[1]),
            "-pe1",
            str(infiles[0]),
            "-pe2",
            str(infiles[1]),
        ]
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

    run_paired(
        "--threads 2 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        in1="big.1.fq",
        in2="big.2.fq",
        expected1="out.1.fastq",
        expected2="out.2.fastq",
        aligners=BACK_ALIGNERS,
        files_equal=False,
        callback=check_summary,
        output_dir=tmp_path,
    )


def test_sam(tmp_path):
    run_interleaved(
        "-a TTAGACATAT -A CAGTGGAGTA -m 14 --output-format sam",
        "paired.sam",
        "paired_insert.sam",
        aligners=["insert"],
        output_dir=tmp_path,
    )


# def test_long_reads():
#    run_paired(''
#        in1='long.1.fq', in2='long.2.fq',
#        expected1=, expected2=,
#        aligners=BACK_ALIGNERS)


def test_issue68(tmp_path):
    run_paired(
        "--error-rate 0.20 --insert-match-error-rate 0.30 --minimum-length 20 --aligner insert -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        in1="issue68.1.fq",
        in2="issue68.2.fq",
        expected1="issue68.1.fq",
        expected2="issue68.2.fq",
        aligners=["insert"],
        output_dir=tmp_path,
    )
