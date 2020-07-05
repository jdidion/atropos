import os
from subprocess import CalledProcessError
import tempfile
import urllib.request

import subby


def assert_files_equal(path1, path2, **format_vars):
    def cmp(p1):
        try:
            subby.sub("diff -u {0} {1}".format(p1, path2))
        except CalledProcessError as e:
            raise AssertionError(
                f"Files not equal: {path1} != {path2}\nDiff: <{e.output}>"
            )

    if format_vars:
        # read the file, substitute in vars, and write it out to a tempfile
        with open(path1, "rt") as inp:
            formatted = inp.read().format(**format_vars)
            tmppath = tempfile.mkstemp()[1]
            with open(tmppath, "wt") as out:
                out.write(formatted)
            try:
                cmp(tmppath)
            finally:
                os.remove(tmppath)
    else:
        cmp(path1)


def no_internet(url="https://github.com", timeout: int = 10):
    """
    Tests whether there's no internet connection available.
    """
    try:
        urllib.request.urlopen(url, timeout=timeout).info()
        return False
    except:
        return True
