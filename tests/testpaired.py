# coding: utf-8
from __future__ import print_function, division, absolute_import

import shutil
import os
from nose.tools import raises
from atropos.scripts import atropos
from .utils import run, files_equal, datapath, cutpath, redirect_stderr, temporary_path

BACK_ALIGNERS = ('adapter', 'insert')

def run_paired(params, in1, in2, expected1, expected2, aligners=('adapter',), callback=None):
    if type(params) is str:
        params = params.split()
    for aligner in aligners:
        with temporary_path('tmp1-' + expected1.format(aligner=aligner)) as p1:
            with temporary_path('tmp2-' + expected2.format(aligner=aligner)) as p2:
                print(aligner)
                p = params.copy()
                p += ['--aligner', aligner, '-o', p1, '-p', p2]
                p += [datapath(in1.format(aligner=aligner)), datapath(in2.format(aligner=aligner))]
                assert atropos.main(p) is None
                assert files_equal(cutpath(expected1.format(aligner=aligner)), p1)
                assert files_equal(cutpath(expected2.format(aligner=aligner)), p2)
                if callback:
                    callback(aligner)

def run_interleaved(params, inpath, expected, aligners=('adapter',)):
    if type(params) is str:
        params = params.split()
    for aligner in aligners:
        with temporary_path(expected.format(aligner=aligner)) as tmp:
            p = params.copy()
            p += ['--aligner', aligner, '--interleaved', '-o', tmp, datapath(inpath.format(aligner=aligner))]
            assert atropos.main(p) is None
            assert files_equal(cutpath(expected.format(aligner=aligner)), tmp)

def run_interleaved2(params, inpath, expected1, expected2, aligners=('adapter',)):
    assert False  # unused function
    if type(params) is str:
        params = params.split()
    for aligner in aligners:
        with temporary_path('tmp1-' + expected1.format(aligner=aligner)) as p1:
            with temporary_path('tmp2-' + expected2.format(aligner=aligner)) as p2:
                p = params.copy()
                p += ['--aligner', aligner, '--interleaved', '-o', p1, '-p', p2]
                p += [datapath(inpath.format(aligner=aligner))]
                assert atropos.main(p) is None
                assert files_equal(cutpath(expected.format(aligner=aligner)), p1)
                assert files_equal(cutpath(expected.format(aligner=aligner)), p2)

def test_paired_separate():
    '''test separate trimming of paired-end reads'''
    run('-a TTAGACATAT', 'paired-separate.1.fastq', 'paired.1.fastq')
    run('-a CAGTGGAGTA', 'paired-separate.2.fastq', 'paired.2.fastq')

def test_paired_end_legacy():
    '''--paired-output, not using -A/-B/-G'''
    # the -m 14 filters out one read, which should then also be filtered out in the second output file
    # -q 10 should not change anything: qualities in file 1 are high enough,
    # qualities in file 2 should not be inspected.
    run_paired('-a TTAGACATAT -m 14 -q 10',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='paired.m14.1.fastq', expected2='paired.m14.2.fastq'
    )

def test_untrimmed_paired_output():
    with temporary_path("tmp-untrimmed.1.fastq") as untrimmed1:
        with temporary_path("tmp-untrimmed.2.fastq") as untrimmed2:
            def callback(aligner):
                assert files_equal(cutpath('paired-untrimmed.1.fastq'), untrimmed1)
                assert files_equal(cutpath('paired-untrimmed.2.fastq'), untrimmed2)
            run_paired(['-a', 'TTAGACATAT',
                '--untrimmed-output', untrimmed1,
                '--untrimmed-paired-output', untrimmed2],
                in1='paired.1.fastq', in2='paired.2.fastq',
                expected1='paired-trimmed.1.fastq', expected2='paired-trimmed.2.fastq',
                callback=callback
            )

def test_explicit_format_with_paired():
    # Use --format=fastq with input files whose extension is .txt
    with temporary_path("paired.1.txt") as txt1:
        with temporary_path("paired.2.txt") as txt2:
            shutil.copyfile(datapath("paired.1.fastq"), txt1)
            shutil.copyfile(datapath("paired.2.fastq"), txt2)
            run_paired('--format=fastq -a TTAGACATAT -m 14',
                in1=txt1, in2=txt2,
                expected1='paired.m14.1.fastq',
                expected2='paired.m14.2.fastq'
            )

def test_no_trimming_legacy():
    # make sure that this doesn't divide by zero
    atropos.main(['-a', 'XXXXX', '-o', '/dev/null', '-p', '/dev/null', datapath('paired.1.fastq'), datapath('paired.2.fastq')])

def test_no_trimming():
    # make sure that this doesn't divide by zero
    atropos.main(['-a', 'XXXXX', '-A', 'XXXXX', '-o', '/dev/null', '-p', '/dev/null', datapath('paired.1.fastq'), datapath('paired.2.fastq')])

@raises(SystemExit)
def test_missing_file():
    with redirect_stderr():
        atropos.main(['-a', 'XX', '--paired-output', 'out.fastq', datapath('paired.1.fastq')])

@raises(SystemExit)
def test_first_too_short():
    with temporary_path("truncated.1.fastq") as trunc1:
        # Create a truncated file in which the last read is missing
        with open(datapath('paired.1.fastq')) as f:
            lines = f.readlines()
            lines = lines[:-4]
        with open(trunc1, 'w') as f:
            f.writelines(lines)
        with redirect_stderr():
            atropos.main('-a XX --paired-output out.fastq'.split() + [trunc1, datapath('paired.2.fastq')])

@raises(SystemExit)
def test_second_too_short():
    with temporary_path("truncated.2.fastq") as trunc2:
        # Create a truncated file in which the last read is missing
        with open(datapath('paired.2.fastq')) as f:
            lines = f.readlines()
            lines = lines[:-4]
        with open(trunc2, 'w') as f:
            f.writelines(lines)
        with redirect_stderr():
            atropos.main('-a XX --paired-output out.fastq'.split() + [datapath('paired.1.fastq'), trunc2])

@raises(SystemExit)
def test_unmatched_read_names():
    with temporary_path("swapped.1.fastq") as swapped:
        try:
            # Create a file in which reads 2 and are swapped
            with open(datapath('paired.1.fastq')) as f:
                lines = f.readlines()
                lines = lines[0:4] + lines[8:12] + lines[4:8] + lines[12:]
            with open(swapped, 'w') as f:
                f.writelines(lines)
            with redirect_stderr():
                print('here')
                atropos.main('-a XX -o out1.fastq --paired-output out2.fastq'.split() + [swapped, datapath('paired.2.fastq')])
        finally:
            os.remove('out1.fastq')
            os.remove('out2.fastq')
            
def test_legacy_minlength():
    '''Ensure -m is not applied to second read in a pair in legacy mode'''
    run_paired('-a XXX -m 27',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='paired-m27.1.fastq', expected2='paired-m27.2.fastq'
    )

def test_paired_end():
    '''single-pass paired-end with -m'''
    run_paired('-a TTAGACATAT -A CAGTGGAGTA -m 14',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='paired_{aligner}.1.fastq', expected2='paired_{aligner}.2.fastq',
        aligners=BACK_ALIGNERS
    )

def test_paired_anchored_back_no_indels():
    run_paired("-a BACKADAPTER$ -A BACKADAPTER$ -N --no-indels",
        in1='anchored-back.fasta', in2='anchored-back.fasta',
        expected1='anchored-back.fasta', expected2="anchored-back.fasta"
    )

def test_paired_end_qualtrim():
    '''single-pass paired-end with -q and -m'''
    run_paired('-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='pairedq.1.fastq', expected2='pairedq.2.fastq',
        aligners=BACK_ALIGNERS
    )

def test_paired_end_qualtrim_swapped():
    '''single-pass paired-end with -q and -m, but files swapped'''
    run_paired('-q 20 -a CAGTGGAGTA -A TTAGACATAT -m 14',
        in1='paired.2.fastq', in2='paired.1.fastq',
        expected1='pairedq.2.fastq', expected2='pairedq.1.fastq',
        aligners=BACK_ALIGNERS
    )

def test_paired_end_cut():
    run_paired('-u 3 -u -1 -U 4 -U -2',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='pairedu.1.fastq', expected2='pairedu.2.fastq'
    )

def test_paired_end_A_only():
    run_paired('-A CAGTGGAGTA',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='paired-onlyA.1.fastq', expected2='paired-onlyA.2.fastq'
    )

def test_discard_untrimmed():
    # issue #146
    # the first adapter is a sequence cut out from the first read
    run_paired('-a CTCCAGCTTAGACATATC -A XXXXXXXX --discard-untrimmed',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='empty.fastq', expected2='empty.fastq'
    )

def test_discard_trimmed():
    run_paired('-A C -O 1 --discard-trimmed',  # applies everywhere
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='empty.fastq', expected2='empty.fastq'
    )

def test_interleaved():
    '''single-pass interleaved paired-end with -q and -m'''
    run_interleaved('-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
        inpath='interleaved.fastq', expected='interleaved.fastq',
        aligners=BACK_ALIGNERS
    )

@raises(SystemExit)
def test_interleaved_no_paired_output():
    with temporary_path("temp-paired.1.fastq") as p1:
        with temporary_path("temp-paired.2.fastq") as p2:
            params = '-a XX --interleaved'.split()
            with redirect_stderr():
                params += [ '-o', p1, '-p1', p2, 'paired.1.fastq', 'paired.2.fastq']
                atropos.main(params)

"""
# TODO
def test_interleaved_input_paired_output():
    '''single-pass interleaved paired-end with -q and -m, paired output'''
    run_interleaved2('-q 20 -a TTAGACATAT -A CAGTGGAGTA -m 14 -M 90',
        inpath='interleaved.fastq', expected1='pairedq1.fastq', expected2='pairedq2.fastq'
    )
"""

def test_pair_filter():
    run_paired('--pair-filter=both -a TTAGACATAT -A GGAGTA -m 14',
        in1='paired.1.fastq', in2='paired.2.fastq',
        expected1='paired-filterboth_{aligner}.1.fastq', expected2='paired-filterboth_{aligner}.2.fastq',
        aligners=BACK_ALIGNERS
    )

def test_too_short_paired_output():
    with temporary_path("temp-too-short.1.fastq") as p1:
        with temporary_path("temp-too-short.2.fastq") as p2:
            def callback(aligner):
                assert files_equal(cutpath('paired-too-short.1.fastq'), p1)
                assert files_equal(cutpath('paired-too-short.2.fastq'), p2)
            run_paired('-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output '
                '{0} --too-short-paired-output {1}'.format(p1, p2),
                in1='paired.1.fastq', in2='paired.2.fastq',
                expected1='paired_{aligner}.1.fastq', expected2='paired_{aligner}.2.fastq',
                aligners=BACK_ALIGNERS, callback=callback
            )

def test_too_long_output():
    with temporary_path("temp-too-long.1.fastq") as p1:
        with temporary_path("temp-too-long.2.fastq") as p2:
            def callback(aligner):
                assert files_equal(cutpath('paired_{aligner}.1.fastq'.format(aligner=aligner)), p1)
                assert files_equal(cutpath('paired_{aligner}.2.fastq'.format(aligner=aligner)), p2)
            run_paired('-a TTAGACATAT -A CAGTGGAGTA -M 14 --too-long-output '
                '{0} --too-long-paired-output {1}'.format(p1, p2),
                in1='paired.1.fastq', in2='paired.2.fastq',
                expected1='paired-too-short.1.fastq', expected2='paired-too-short.2.fastq',
                aligners=BACK_ALIGNERS, callback=callback
            )

@raises(SystemExit)
def test_too_short_output_paired_option_missing():
    with temporary_path("temp-too-short.1.fastq") as p1:
        run_paired('-a TTAGACATAT -A CAGTGGAGTA -m 14 --too-short-output '
            '{0}'.format(p1),
            in1='paired.1.fastq', in2='paired.2.fastq',
            expected1='paired.1.fastq', expected2='paired.2.fastq',
            aligners=BACK_ALIGNERS
        )

def test_custom_bisulfite_1():
    run_paired('-a TTAGACATAT -A CAGTGGAGTA -m 14 -q 0 --bisulfite 2,2,1,1',
        in1='paired_bis_{aligner}.1.fastq', in2='paired_bis_{aligner}.2.fastq',
        expected1='paired_bis1_{aligner}.1.fastq', expected2='paired_bis1_{aligner}.2.fastq',
        aligners=BACK_ALIGNERS
    )

def test_custom_bisulfite_2():
    run_paired('-a TTAGACATAT -A CAGTGGAGTA -m 10 -q 0 --bisulfite 20,20,1,1;0,0,0,0',
        in1='paired_bis_{aligner}.1.fastq', in2='paired_bis_{aligner}.2.fastq',
        expected1='paired_bis2_{aligner}.1.fastq', expected2='paired_bis2_{aligner}.2.fastq',
        aligners=BACK_ALIGNERS
    )
