import os
from subprocess import CalledProcessError
import tempfile
import urllib.request

import subby

from xphyle import open_


def assert_files_equal(path1, path2, **format_vars):
    with open_(path1, 'r') as i1, open_(path2, 'r') as i2:
        # write contents to temp files in case the files are compressed
        content1 = i1.read()
        content2 = i2.read()

    temp1 = tempfile.mkstemp()[1]
    temp2 = tempfile.mkstemp()[1]

    try:
        with open(temp1, "w") as out:
            if format_vars:
                out.write(content1.format(**format_vars))
            else:
                out.write(content1)
        with open(temp2, "w") as out:
            out.write(content2)
        subby.sub("diff -u {0} {1}".format(temp1, temp2))
        return True
    except CalledProcessError as e:
        raise AssertionError(
            f"Files not equal: {path1} != {path2}\nDiff: <{e.output}>"
        )
    finally:
        os.remove(temp1)
        os.remove(temp2)


def no_internet(url="https://github.com", timeout: int = 10):
    """
    Tests whether there's no internet connection available.
    """
    try:
        urllib.request.urlopen(url, timeout=timeout).info()
        return False
    except:
        return True
