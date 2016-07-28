"""
Open compressed files transparently.
"""
import gzip
import sys
import io
import os
import re
from subprocess import Popen, PIPE

try:
    import bz2
except ImportError:
    bz2 = None

try:
    import lzma
except ImportError:
    lzma = None

def open_output(filename, mode='w'):
    """
    Replacement for the "open" function that can also open files that have
    been compressed with gzip, bzip2 or xz. If the filename is '-', standard
    output (mode 'w') or input (mode 'r') is returned. If the filename ends
    with .gz, the file is opened with a pipe to the gzip program. If that
    does not work, then gzip.open() is used (the gzip module is slower than
    the pipe to the gzip program). If the filename ends with .bz2, it's
    opened as a bz2.BZ2File. Otherwise, the regular open() is used.

    mode can be: 'a', 'ab', 'wt', or 'wb'
    Instead of 'rt' and 'wt', 'r' and 'w' can be used as abbreviations.

    In Python 2, the 't' and 'b' characters are ignored.

    Append mode ('a') is unavailable with BZ2 compression and will raise an error.
    """
    if mode == 'w':
        mode = 'wt'
    elif mode == 'a':
        mode = 'at'
    if mode not in ('wt', 'wb', 'at', 'ab'):
        raise ValueError("mode '{0}' not supported".format(mode))
    if not isinstance(filename, str):
        raise ValueError("the filename must be a string")

    # standard input and standard output handling
    if filename == '-':
        return dict(
            wt=sys.stdout,
            wb=sys.stdout.buffer)[mode]
    
    return open(filename, mode)

def xopen(filename, mode='r'):
    """
    Replacement for the "open" function that can also open files that have
    been compressed with gzip, bzip2 or xz. If the filename is '-', standard
    output (mode 'w') or input (mode 'r') is returned. If the filename ends
    with .gz, the file is opened with a pipe to the gzip program. If that
    does not work, then gzip.open() is used (the gzip module is slower than
    the pipe to the gzip program). If the filename ends with .bz2, it's
    opened as a bz2.BZ2File. Otherwise, the regular open() is used.

    mode can be: 'rt', 'rb', 'a', 'wt', or 'wb'
    Instead of 'rt' and 'wt', 'r' and 'w' can be used as abbreviations.

    In Python 2, the 't' and 'b' characters are ignored.

    Append mode ('a') is unavailable with BZ2 compression and will raise an error.
    """
    if mode == 'r':
        mode = 'rt'
    elif mode == 'w':
        mode = 'wt'
    elif mode == 'a':
        mode = 'at'
    if mode not in ('rt', 'rb', 'wt', 'wb', 'at', 'ab'):
        raise ValueError("mode '{0}' not supported".format(mode))
    if not isinstance(filename, str):
        raise ValueError("the filename must be a string")

    # standard input and standard output handling
    if filename == '-':
        return dict(
            rt=sys.stdin,
            wt=sys.stdout,
            rb=sys.stdin.buffer,
            wb=sys.stdout.buffer)[mode]

    if filename.endswith('.bz2'):
        if bz2 is None:
            raise ImportError("Cannot open bz2 files: The bz2 module is not available")
        if 't' in mode:
            return io.TextIOWrapper(bz2.BZ2File(filename, mode[0]))
        else:
            return bz2.BZ2File(filename, mode)
    elif filename.endswith('.xz'):
        if lzma is None:
            raise ImportError("Cannot open xz files: The lzma module is not available "
                "(use Python 3.3 or newer)")
        return lzma.open(filename, mode)
    elif filename.endswith('.gz'):
        if 't' in mode:
            # gzip.open in Python 3.2 does not support modes 'rt' and 'wt''
            return io.TextIOWrapper(gzip.open(filename, mode[0]))
        else:
            if 'r' in mode:
                return io.BufferedReader(gzip.open(filename, mode))
            else:
                return io.BufferedWriter(gzip.open(filename, mode))
    else:
        return open(filename, mode)

compressors = {
    ".gz"  : gzip,
    ".bz2" : bz2,
    ".xz"  : lzma
}

def get_compressor(filename):
    for ext, compressor in compressors.items():
        if filename.endswith(ext):
            return compressor
    return None
