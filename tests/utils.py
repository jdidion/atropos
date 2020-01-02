from subprocess import CalledProcessError
import urllib.request

import subby


def assert_files_equal(path1, path2):
    try:
        subby.sub("diff -u {0} {1}".format(path1, path2))
    except CalledProcessError as e:
        raise AssertionError(
            f"Files not equal: {path1} != {path2}\nDiff: <{e.output}>"
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
