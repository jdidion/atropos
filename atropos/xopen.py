"""
Open compressed files transparently.
"""
import sys
import io
import os

from .compression import get_file_opener

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

def xopen(filename, mode='r', use_system=True):
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
    
    file_opener = get_file_opener(filename)
    if file_opener:
        return file_opener(filename, mode, use_system)
    else:
         return open(filename, mode)
