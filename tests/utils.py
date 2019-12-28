from importlib import import_module
import urllib.request

import subby


def assert_files_equal(path1, path2):
    try:
        subby.sub("diff -u {0} {1}".format(path1, path2))
    except subby.CalledProcessError as e:
        raise AssertionError(
            f"Files not equal: {path1} != {path2}\n"
            f"Diff: <{e.output.decode('utf-8')}>"
        )


def no_internet(url="https://github.com"):
    """
    Tests whether there's no internet connection available.
    """
    try:
        urllib.request.urlopen(url).info()
        return False
    except:
        return True


def no_import(lib):
    """
    Tests whether a library is importable.
    """
    try:
        mod = import_module(lib)
        return mod is None
    except ImportError:
        return True
