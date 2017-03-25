"""File compression/decompression functions.
"""
import bz2
import gzip
import io
import lzma
import os
from subprocess import Popen, PIPE

COMPRESSORS = {
    ".gz"  : gzip,
    ".bz2" : bz2,
    ".xz"  : lzma
}
"""Mapping of file extension to python compression library."""

class GzipWriter:
    """Wrapper for a process that uses the system gzip program to compress
    bytes.
    
    Args:
        path: The path of the output file.
        mode: The file open mode.
    """
    def __init__(self, path, mode='w'):
        self.name = path
        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, 'w')
        self.closed = False
        try:
            # Setting close_fds to True is necessary due to
            # http://bugs.python.org/issue12786
            self.process = Popen(
                [get_program_path('gzip')], stdin=PIPE, stdout=self.outfile,
                stderr=self.devnull, close_fds=True)
        except IOError:
            self.outfile.close()
            self.devnull.close()
            raise
    
    def readable(self):
        return False
    
    def writable(self):
        return True
    
    def seekable(self):
        return False
    
    def write(self, arg):
        self.process.stdin.write(arg)
    
    def flush(self):
        self.process.stdin.flush()
    
    def close(self):
        self.closed = True
        self.process.stdin.close()
        retcode = self.process.wait()
        self.outfile.close()
        self.devnull.close()
        if retcode != 0:
            raise IOError(
                "Output gzip process terminated with exit code {0}".format(
                    retcode))
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.close()

class GzipReader:
    """Wrapper for a process that uses the system gzip program to decompress
    bytes.
    
    Args:
        path: The path of the input file.
    """
    def __init__(self, path):
        self.name = path
        self.process = Popen([get_program_path('gzip'), '-cd', path], stdout=PIPE)
        self.closed = False
    
    def readable(self):
        return True
    
    def writable(self):
        return False
    
    def seekable(self):
        return False
    
    def flush(self):
        pass
    
    def close(self):
        if self.closed:
            return
        self.closed = True
        retcode = self.process.poll()
        if retcode is None:
            # still running
            self.process.terminate()
        self._raise_if_error()
    
    def __iter__(self):
        for line in self.process.stdout:
            yield line
        self.process.wait()
        self._raise_if_error()
    
    def _raise_if_error(self):
        """Raise EOFError if process is not running anymore and the exit code
        is nonzero.
        """
        retcode = self.process.poll()
        if retcode is not None and retcode != 0:
            raise EOFError(
                "gzip process returned non-zero exit code {0}. Is the "
                "input file truncated or corrupt?".format(retcode))
    
    def read(self, *args):
        data = self.process.stdout.read(*args)
        if len(args) == 0 or args[0] <= 0:
            # wait for process to terminate until we check the exit code
            self.process.wait()
        self._raise_if_error()
        return data
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_info):
        self.close()

def can_use_system_compression():
    """Whether the system gzip program is available.
    """
    return get_program_path("gzip") is not None

def get_compressor(filename):
    """Returns the python compression library for a file based on its extension.
    """
    ext = os.path.splitext(filename)[1]
    if ext in COMPRESSORS:
        return COMPRESSORS[ext]
    return None

def open_gzip_file(filename, mode, use_system=True):
    """Open a gzip file, preferring the system gzip program if `use_system`
    is True, falling back to the gzip python library.
    
    Args:
        mode: The file open mode.
        use_system: Whether to try to use the system gzip program.
    """
    if use_system:
        try:
            if 'r' in mode:
                gzfile = GzipReader(filename)
            else:
                gzfile = GzipWriter(filename)
            if 't' in mode:
                gzfile = io.TextIOWrapper(gzfile)
            return gzfile
        except:
            pass
    
    gzfile = gzip.open(filename, mode)
    if 'b' in mode:
        if 'r' in mode:
            gzfile = io.BufferedReader(gzfile)
        else:
            gzfile = io.BufferedWriter(gzfile)
    return gzfile

def open_bzip_file(filename, mode, **kwargs):
    """Open a bzip file.
    """
    if 't' in mode:
        return io.TextIOWrapper(bz2.BZ2File(filename, mode[0]))
    else:
        return bz2.BZ2File(filename, mode)

def open_lzma_file(filename, mode, **kwargs):
    """Open a LZMA (xz) file.
    """
    return lzma.open(filename, mode)

FILE_OPENERS = {
    ".gz"  : open_gzip_file,
    ".bz2" : open_bzip_file,
    ".xz"  : open_lzma_file,
}
"""Mapping of file extensions to file opener functions."""

def get_file_opener(filename):
    """Returns the file opener for a filename based on its extension.
    """
    ext = os.path.splitext(filename)[1]
    if ext in FILE_OPENERS:
        return FILE_OPENERS[ext]
    return None

PROGRAM_CACHE = {}
"""Cache for resolved program paths."""

def get_program_path(program):
    """Get the path of a program on the system. Returns from the cache if it
    was previously resolved, otherwise finds it in the system path and caches
    it.
    
    Args:
        program: The command name to resolve.
    
    Returns:
        The fully resolved path of the command.
    """
    if program in PROGRAM_CACHE:
        return PROGRAM_CACHE[program]
    
    def is_exe(fpath):
        """Returns True if `fpath` is executable.
        """
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
    exe_file = None
    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            exe_file = program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                break
    
    PROGRAM_CACHE[program] = exe_file
    return exe_file

def open_compressed_file(filename, mode):
    """Open a compressed file, determining the compression type based on the
    file name.
    
    Args:
        filename: The file to open.
        mode: The file open mode.
    
    Returns:
        The opened file.
    """
    ext = os.path.splitext(filename)
    opener = get_file_opener(ext)
    if not opener:
        raise ValueError("{} is not a recognized compression format")
    return opener(filename, mode)

def splitext_compressed(name):
    """Split the filename and extensions of a file that potentially has two
    extensions - one for the file type (e.g. 'fq') and one for the compression
    type (e.g. 'gz').
    
    Args:
        name: The filename.
    
    Returns:
        A tuple (name, ext1, ext2), where ext1 is the filetype extension and
        ext2 is the compression type extension, or None.
    """
    ext1 = ext2 = None
    for ext in COMPRESSORS:
        if name.endswith(ext):
            ext2 = ext
            name = name[:-len(ext)]
            break
    name, ext1 = os.path.splitext(name)
    return (name, ext1, ext2)
