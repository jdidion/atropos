# coding: utf-8
from __future__ import print_function, division, absolute_import

import sys, os
from contextlib import contextmanager
from atropos.scripts import atropos

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
    os.remove(path)


def datapath(path):
    return os.path.join(os.path.dirname(__file__), 'data', path)


def cutpath(path):
    return os.path.join(os.path.dirname(__file__), 'cut', path)


def files_equal(path1, path2):
    #return os.system("diff -u {0} {1}".format(path1, path2)) == 0
    with open(path1,'r') as i1, open(path2, 'r') as i2:
        print("<[{}]>".format(i1.read()))
        print("<[{}]>".format(i2.read()))
    from subprocess import check_output, CalledProcessError
    try:
        check_output("diff -u {0} {1}".format(path1, path2), shell=True)
        return True
    except CalledProcessError as e:
        print("Diff: <{}>".format(e.output.decode("utf-8")))
        return False


def run(params, expected, inpath, inpath2=None):
    if type(params) is str:
        params = params.split()
    with temporary_path(expected) as tmp_fastaq:
        params += ['-o', tmp_fastaq ] # TODO not parallelizable
        params += [ datapath(inpath) ]
        if inpath2:
            params += [ datapath(inpath2) ]
        print(params)
        assert atropos.main(params) is None
        # TODO redirect standard output
        assert files_equal(cutpath(expected), tmp_fastaq)
    # TODO diff log files

def approx_equal(a, b, tol):
    return abs(a-b) <= tol
