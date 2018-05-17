# coding: utf-8
from contextlib import contextmanager
from importlib import import_module
import os
import sys
import traceback
import urllib.request
from atropos.commands import get_command


@contextmanager
def redirect_stderr():
    "Send stderr to stdout. Nose doesn't capture stderr, yet."
    old_stderr = sys.stderr
    sys.stderr = sys.stdout
    yield

    sys.stderr = old_stderr


@contextmanager
def temporary_path(name):
    directory = os.path.join(os.path.dirname(__file__), 'testtmp')
    if not os.path.isdir(directory):
        os.mkdir(directory)
    path = os.path.join(directory, name)
    yield path

    if os.path.exists(path):
        os.remove(path)


def datapath(path):
    return os.path.join(os.path.dirname(__file__), 'data', path)


def cutpath(path):
    return os.path.join(os.path.dirname(__file__), 'cut', path)


def files_equal(path1, path2):
    # return os.system("diff -u {0} {1}".format(path1, path2)) == 0
    with open(path1, 'r') as i1, open(path2, 'r') as i2:
        print("<[{}]>".format(i1.read()))
        print("<[{}]>".format(i2.read()))
    from subprocess import check_output, CalledProcessError

    try:
        check_output("diff -u {0} {1}".format(path1, path2), shell=True)
        return True

    except CalledProcessError as e:
        print("Diff: <{}>".format(e.output.decode("utf-8")))
        return False


def run(
    params,
    expected,
    inpath=None,
    inpath2=None,
    qualfile=None,
    interleaved_input=False,
    interleaved_output=False,
    sra_accn=None,
):
    if type(params) is str:
        params = params.split()
    with temporary_path(expected) as tmp_fastaq:
        if sra_accn:
            params += ['-sra', sra_accn]
        elif interleaved_input:
            params += ['-l', inpath]
        elif inpath2:
            params += ['-pe1', datapath(inpath)]
            params += ['-pe2', datapath(inpath2)]
        else:
            params += ['-se', datapath(inpath)]
            if qualfile:
                params += ['-sq', datapath(qualfile)]
        if interleaved_output:
            params += ['-L', tmp_fastaq]
        else:
            params += ['-o', tmp_fastaq]  # TODO not parallelizable
        # print(params)
        command = get_command('trim')
        retcode, summary = command.execute(params)
        assert summary is not None
        assert isinstance(summary, dict)
        if 'exception' in summary and summary['exception'] is not None:
            assert retcode != 0
            err = summary['exception']
            traceback.print_exception(* err['details'])
            raise Exception("Unexpected error: {}".format(err['message']))

        else:
            assert retcode == 0
        # TODO redirect standard output
        assert os.path.exists(tmp_fastaq)
        assert files_equal(cutpath(expected), tmp_fastaq)




# TODO diff log files
def approx_equal(a, b, tol):
    return abs(a - b) <= tol


def no_internet(url="https://github.com"):
    """Test whether there's no internet connection available.
    """
    try:
        urllib.request.urlopen(url).info()
        return False

    except:
        return True


def no_import(lib):
    """Test whether a library is importable.
    """
    try:
        mod = import_module(lib)
        return mod is None

    except:
        return True
