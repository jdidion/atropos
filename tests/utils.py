# coding: utf-8
from contextlib import contextmanager
from importlib import import_module
import os
from subprocess import check_output, CalledProcessError
import sys
import tempfile
import traceback
import urllib.request
from atropos.commands import get_command
from atropos.io import xopen


@contextmanager
def redirect_stderr():
    "Send stderr to stdout. Nose doesn't capture stderr, yet."
    old_stderr = sys.stderr
    sys.stderr = sys.stdout
    yield

    sys.stderr = old_stderr


def datapath(path):
    return os.path.join(os.path.dirname(__file__), "data", path)


def cutpath(path):
    return os.path.join(os.path.dirname(__file__), "cut", path)


def files_equal(path1, path2):
    temp1 = tempfile.mkstemp()[1]
    temp2 = tempfile.mkstemp()[1]
    try:
        with xopen(path1, "r") as i1, xopen(path2, "r") as i2:
            # write contents to temp files in case the files are compressed
            content1 = i1.read()
            content2 = i2.read()
        print("<[{}]>".format(content1))
        print("<[{}]>".format(content2))
        with open(temp1, "w") as out:
            out.write(content1)
        with open(temp2, "w") as out:
            out.write(content2)
        check_output("diff -u {0} {1}".format(temp1, temp2), shell=True)
        return True
    except CalledProcessError as e:
        print("Diff: <{}>".format(e.output.decode("utf-8")))
        return False
    finally:
        os.remove(temp1)
        os.remove(temp2)


def run(
    tmp_path,
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
    tmp_fastaq = str(tmp_path / expected)
    if sra_accn:
        params += ["-sra", sra_accn]
    elif interleaved_input:
        params += ["-l", inpath]
    elif inpath2:
        params += ["-pe1", datapath(inpath)]
        params += ["-pe2", datapath(inpath2)]
    else:
        params += ["-se", datapath(inpath)]
        if qualfile:
            params += ["-sq", datapath(qualfile)]
    if interleaved_output:
        params += ["-L", tmp_fastaq]
    else:
        params += ["-o", tmp_fastaq]  # TODO not parallelizable
    # print(params)
    command = get_command("trim")
    retcode, summary = command.execute(params)
    assert summary is not None
    assert isinstance(summary, dict)
    if "exception" in summary and summary["exception"] is not None:
        assert retcode != 0
        err = summary["exception"]
        traceback.print_exception(*err["details"])
        raise Exception("Unexpected error: {}".format(err["message"]))

    else:
        assert retcode == 0
    # TODO redirect standard output
    assert os.path.exists(tmp_fastaq)
    assert files_equal(cutpath(expected), tmp_fastaq)


# TODO diff log files
def approx_equal(a, b, tol):
    return abs(a - b) <= tol


def no_internet(url="https://github.com"):
    """Test whether there's no internet connection available."""
    try:
        urllib.request.urlopen(url).info()
        return False

    except:
        return True


def no_import(lib):
    """Test whether a library is importable."""
    try:
        mod = import_module(lib)
        return mod is None

    except:
        return True
