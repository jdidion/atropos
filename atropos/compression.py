import io
import os
import sys
from subprocess import Popen, PIPE

# Inline compression

import gzip
import bz2
import lzma

compressors = {
    ".gz"  : gzip,
    ".bz2" : bz2,
    ".xz"  : lzma
}

def get_compressor(filename):
    ext = os.path.splitext(filename)[1]
    if ext in compressors:
        return compressors[ext]
    return None

# File compression
# Using a system process can be faster, but only works
# on a system with gzip installed.

class GzipWriter:
    def __init__(self, path, mode='w'):
        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, 'w')
        self.closed = False
        try:
            # Setting close_fds to True is necessary due to
            # http://bugs.python.org/issue12786
            self.process = Popen([get_program_path('gzip')], stdin=PIPE, stdout=self.outfile,
                stderr=self.devnull, close_fds=True)
        except IOError as e:
            self.outfile.close()
            self.devnull.close()
            raise
    
    def readable(self): return False
    def writable(self): return True
    def seekable(self): return False
    
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
            raise IOError("Output gzip process terminated with exit code {0}".format(retcode))

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

class GzipReader:
    def __init__(self, path):
        self.process = Popen([get_program_path('gzip'), '-cd', path], stdout=PIPE)
        self.closed = False
    
    def readable(self): return True
    def writable(self): return False
    def seekable(self): return False
    def flush(self): pass
    
    def close(self):
        self.close = True
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
        """
        Raise EOFError if process is not running anymore and the
        exit code is nonzero.
        """
        retcode = self.process.poll()
        if retcode is not None and retcode != 0:
            raise EOFError("gzip process returned non-zero exit code {0}. Is the "
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

def open_gzip_file(filename, mode, use_system=True):
    if use_system:
        try:
            if 'r' in mode:
                z = GzipReader(filename)
            else:
                z = GzipWriter(filename)
            if 't' in mode:
                z = io.TextIOWrapper(z)
            return z
        except:
            pass
    
    z = gzip.open(filename, mode)
    if 'b' in mode:
        if 'r' in mode:
            z = io.BufferedReader(z)
        else:
            z = io.BufferedWriter(z)
    return z

# TODO: experiment with system process-based bzip and lzma readers/writers

def open_bzip_file(filename, mode, use_system=False):
    if 't' in mode:
        return io.TextIOWrapper(bz2.BZ2File(filename, mode[0]))
    else:
        return bz2.BZ2File(filename, mode)

def open_lzma_file(filename, mode, use_system=False):
    return lzma.open(filename, mode)

file_openers = {
    ".gz"  : open_gzip_file,
    ".bz2" : open_bzip_file,
    ".xz"  : open_lzma_file
}

def can_use_system_compression():
    return get_program_path("gzip") is not None

def get_file_opener(filename):
    ext = os.path.splitext(filename)[1]
    if ext in file_openers:
        return file_openers[ext]
    return None

program_cache = {}

def get_program_path(program):
    if program in program_cache:
        return program_cache[program]
    
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
    exe_file = None
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            exe_file = program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                break
    
    program_cache[program] = exe_file
    return exe_file

def open_compressed_file(filename, mode):
    ext = os.path.splitext(filename)
    opener = get_file_opener(ext)
    if not opener:
        raise Exception("{} is not a recognized compression format")
    return opener(filename, mode)
