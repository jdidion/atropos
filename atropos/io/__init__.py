"""Open compressed files transparently.
"""
import errno
import io
import os
import sys

from atropos.io.compression import get_file_opener

STDOUT = '-'
STDERR = '_'

def abspath(path):
    """Returns the user home-resolved absolute path.
    """
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
    """Checks that a path exists, is of the specified type, and allows the
    specified access.
    
    Args:
        ptype: 'f' for file or 'd' for directory.
        access (int): One of the access values from :module:`os`
    
    Raises:
        IOError if the path does not exist, is not of the specified type, or
        doesn't allow the specified access.
    """
    if (
            ptype == 'f' and not path.startswith("/dev/") and
            not os.path.isfile(path)):
        raise IOError(errno.EISDIR, "{} is not a file".format(path), path)
    elif ptype == 'd' and not os.path.isdir(path):
        raise IOError(errno.ENOTDIR, "{} is not a directory".format(path), path)
    elif not os.path.exists(path):
        raise IOError(errno.ENOENT, "{} does not exist".format(path), path)
    if access is not None and not os.access(path, access):
        raise IOError(errno.EACCES, "{} is not accessable".format(path), path)
    return path

def check_writeable(rawpath, ptype=None):
    """Resolves the absolute path. Raises an IOError if the path is not
    writable.
    
    Args:
        rawpath: The path to resolve/check.
        ptype: The path type (f=file, d=directory).
    """
    if rawpath in (STDOUT, STDERR):
        return rawpath
    rawpath = abspath(rawpath)
    try:
        path = resolve_path(rawpath)
        check_path(path, ptype, os.W_OK)
    except IOError:
        dirpath = os.path.dirname(rawpath)
        if os.path.exists(dirpath):
            check_path(dirpath, 'd', os.W_OK)
        else:
            os.makedirs(dirpath)
        path = os.path.join(dirpath, os.path.basename(rawpath))
    return path

def open_output(filename, mode='w', context_wrapper=False):
    """Replacement for the "open" function that is only for writing text files.
    If the filename is '-', standard output (mode 'w').
    
    Args:
        filename: The file to open.
        mode: The file open mode; can be: 'a', 'ab', 'wt', or 'wb'.
            Instead of 'wt', 'w' can be used as an abbreviation.
        context_wrapper: Whether to wrap the file in a context manager.
    
    Returns:
        The opened file.
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
        fileobj = sys.stdout if filename == STDOUT else sys.stderr
        if mode == 'wb':
            fileobj = fileobj.buffer
        if context_wrapper:
            class StdWrapper(object):
                """Context manager for stdout/stderr that is no-op on exit.
                """
                def __init__(self, fileobj):
                    self.fileobj = fileobj
                def __enter__(self):
                    return self.fileobj
                def __exit__(self, exception_type, exception_value, traceback):
                    pass
            fileobj = StdWrapper(fileobj)
    else:
        filename = check_writeable(filename, 'f')
        fileobj = open(filename, mode)
    
    return fileobj

def xopen(filename, mode='r', use_system=True):
    """Replacement for the "open" function that can also open files that have
    been compressed with gzip, bzip2 or xz. If the filename is '-', standard
    output (mode 'w') or input (mode 'r') is returned. If the filename ends
    with .gz, the file is opened with a pipe to the gzip program. If that
    does not work, then gzip.open() is used (the gzip module is slower than
    the pipe to the gzip program). If the filename ends with .bz2, it's
    opened as a bz2.BZ2File. Otherwise, the regular open() is used.
    
    Args:
        filename: The file to open.
        mode: The file open mode. Can be: 'rt', 'rb', 'a', 'wt', or 'wb'.
            Instead of 'rt' and 'wt', 'r' and 'w' can be used as abbreviations.
            Append mode ('a') is unavailable with BZ2 compression and will raise
            an error.
        use_system: Whether to use the system compression/decompression program.
    
    Returns:
        The opened file.
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
            fileobj = sys.stdin
        else:
            fileobj = sys.stdout if filename == STDOUT else sys.stderr
        if 'b' in mode:
            fileobj = fileobj.buffer
        return fileobj
    
    file_opener = get_file_opener(filename)
    if file_opener:
        return file_opener(filename, mode, use_system=use_system)
    else:
        return open(filename, mode)
