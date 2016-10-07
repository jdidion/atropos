"""
Open compressed files transparently.
"""
import errno
import io
import os
import sys

from .compression import get_file_opener

STDOUT = '-'
STDERR = '_'

def abspath(path):
    return os.path.abspath(os.path.expanduser(path))
    
def resolve_path(path, parent=None):
    """Resolves the absolute path of the specified file.
    
    Args:
        path (str): Path to resolve.
        parent (str): The directory containing ``path`` if ``path`` is relative.
    
    Returns:
        The absolute path.
    
    Raises:
        IOError: if the path does not exist.
    """
    apath = abspath(path)
    if not os.path.exists(apath) and parent is not None:
        apath = abspath(os.path.join(parent, path))
    if not os.path.exists(apath):
        raise IOError(errno.ENOENT, "%s does not exist" % apath, apath)
    return apath

def check_path(path, ptype=None, access=None):
    """Checks that a path exists, is of the specified type, and allows the specified access.
    
    Args:
        ptype: 'f' for file or 'd' for directory.
        access (int): One of the access values from :module:`os`
    
    Raises:
        IOError if the path does not exist, is not of the specified type, or doesn't allow the
        specified access.
    """
    if ptype == 'f' and not path.startswith("/dev/") and not os.path.isfile(path):
        raise IOError(errno.EISDIR, "{} is not a file".format(path), path)
    elif ptype == 'd' and not os.path.isdir(path):
        raise IOError(errno.ENOTDIR, "{} is not a directory".format(path), path)
    elif not os.path.exists(path):
        raise IOError(errno.ENOENT, "{} does not exist".format(path), path)
    if access is not None and not os.access(path, access):
        raise IOError(errno.EACCES, "{} is not accessable".format(path), path)
    return path

def check_writeable(p, ptype=None):
    if p in (STDOUT, STDERR):
        return p
    p = abspath(p)
    try:
        path = resolve_path(p)
        check_path(path, ptype, os.W_OK)
    except IOError:
        dirpath = os.path.dirname(p)
        if os.path.exists(dirpath):
            check_path(dirpath, "d", os.W_OK)
        else:
            os.makedirs(dirpath)
        path = os.path.join(dirpath, os.path.basename(p))
    return path

def open_output(filename, mode='w', context_wrapper=False):
    """
    Replacement for the "open" function that is only for writing text files.
    If the filename is '-', standard output (mode 'w').
    
    mode can be: 'a', 'ab', 'wt', or 'wb'
    Instead of 'wt', 'w' can be used as an abbreviation.
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
    if filename in (STDOUT, STDERR):
        fh = sys.stdout if filename == STDOUT else sys.stderr
        if mode == 'wb':
            fh = fh.buffer
        if context_wrapper:
            class StdWrapper(object):
                def __init__(self, fh):
                    self.fh = fh
                def __enter__(self):
                    return self.fh
                def __exit__(exception_type, exception_value, traceback):
                    pass
            fh = StdWrapper(fh)
    else:
        filename = check_writeable(filename, 'f')
        fh = open(filename, mode)
    
    return fh

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
    if filename in (STDOUT, STDERR):
        if 'r' in mode:
            fh = sys.stdin
        else:
            fh = sys.stdout if filename == STDOUT else sys.stderr
        if 'b' in mode:
            fh = fh.buffer
        return fh
    
    file_opener = get_file_opener(filename)
    if file_opener:
        return file_opener(filename, mode, use_system)
    else:
         return open(filename, mode)
