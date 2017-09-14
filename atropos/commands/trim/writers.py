"""Classes for formatting and writing trimmed reads to output.
"""
import sys
from atropos.io import STDOUT, xopen, open_output
from atropos.io.compression import splitext_compressed
from atropos.io.seqio import create_seq_formatter
from .filters import NoFilter

class Writers(object):
    """Manages writing to one or more outputs.
    
    Args:
        force_create: Whether empty output files should be created.
    """
    def __init__(self, force_create=[]):
        self.writers = {}
        self.force_create = force_create
        self.suffix = None
    
    def get_writer(self, file_desc, compressed=False):
        """Create the writer for a file descriptor if it does not already
        exist.
        
        Args:
            file_desc: File descriptor. If `compressed==True`, this is a tuple
                (path, mode), otherwise it's only a path.
            compressed: Whether data has already been compressed.
        
        Returns:
            The writer.
        """
        if compressed:
            path, mode = file_desc
        else:
            path = file_desc
        
        if path not in self.writers:
            if self.suffix:
                real_path = add_suffix_to_path(path, self.suffix)
            else:
                real_path = path
            # TODO: test whether O_NONBLOCK allows non-blocking write to NFS
            if compressed:
                self.writers[path] = open_output(real_path, mode)
            else:
                self.writers[path] = xopen(real_path, "w")
        
        return self.writers[path]
    
    def write_result(self, result, compressed=False):
        """Write results to output.
        
        Args:
            result: Dict with keys being file descriptors and values being data
                (either bytes or strings). Strings are expected to already have
                appropriate line-endings.
            compressed: Whether data has already been compressed.
        """
        for file_desc, data in result.items():
            self.write(file_desc, data, compressed)
    
    def write(self, file_desc, data, compressed=False):
        """Write data to output. If the specified path has not been seen before,
        the output is opened.
        
        Args:
            file_desc: File descriptor. If `compressed==True`, this is a tuple
                (path, mode), otherwise it's only a path.
            data: The data to write.
            compressed: Whether data has already been compressed.
        """
        self.get_writer(file_desc, compressed).write(data)
    
    def close(self):
        """Close all outputs.
        """
        for path in self.force_create:
            if path not in self.writers and path != STDOUT:
                with open_output(path, "w"):
                    pass
        for writer in self.writers.values():
            if writer not in (sys.stdout, sys.stderr):
                writer.close()

class Formatters(object):
    """Manages multiple formatters.
    
    Args:
        output: The output file name template.
        seq_formatter_args: Additional arguments to pass to the formatter
            constructor.
    """
    def __init__(self, output, seq_formatter_args):
        self.output = output
        self.multiplexed = output is not None and '{name}' in output
        self.seq_formatter_args = seq_formatter_args
        self.seq_formatters = {}
        self.mux_formatters = {}
        self.info_formatters = []
        self.discarded = 0
    
    def add_seq_formatter(self, filter_type, file1, file2=None):
        """Add a formatter.
        
        Args:
            filter_type: The type of filter that triggers writing with the
                formatter.
            file1, file2: The output file(s).
        """
        self.seq_formatters[filter_type] = create_seq_formatter(
            file1, file2, **self.seq_formatter_args)
    
    def add_info_formatter(self, formatter):
        """Add a formatter for one of the delimited detail files
        (rest, info, wildcard).
        """
        self.info_formatters.append(formatter)
    
    def get_mux_formatter(self, name):
        """Returns the formatter associated with the given name (barcode) when
        running in multiplexed mode.
        """
        assert self.multiplexed
        if name not in self.mux_formatters:
            path = self.output.format(name=name)
            self.mux_formatters[name] = create_seq_formatter(
                path, **self.seq_formatter_args)
        return self.mux_formatters[name]
    
    def get_seq_formatters(self):
        """Returns a set containing all formatters that have handled at least
        one record.
        """
        return (
            set(f for f in self.seq_formatters.values() if f.written > 0) |
            set(f for f in self.mux_formatters.values() if f.written > 0))
    
    def format(self, result, dest, read1, read2=None):
        """Format read(s) and add to a result dict. Also writes info records
        to any registered info formatters.
        
        Args:
            result: The result dict.
            dest: The destination (filter type).
            read1, read2: The read(s).
        """
        if self.multiplexed and (dest == NoFilter) and read1.match:
            name = read1.match.adapter.name
            formatter = self.get_mux_formatter(name)
            formatter.format(result, read1, read2)
        elif dest in self.seq_formatters:
            self.seq_formatters[dest].format(result, read1, read2)
        else:
            self.discarded += 1
        
        for fmtr in self.info_formatters:
            fmtr.format(result, read1)
            if read2:
                fmtr.format(result, read2)
    
    def summarize(self):
        """Returns a summary dict.
        """
        seq_formatters = self.get_seq_formatters()
        return dict(
            records_written=sum(f.written for f in seq_formatters),
            bp_written=[
                sum(f.read1_bp for f in seq_formatters),
                sum(f.read2_bp for f in seq_formatters)
            ])

class DelimFormatter(object):
    """Base class for formatters that write to a delimited file.
    
    Args:
        path: The output file path.
        delim: The field delimiter.
    """
    def __init__(self, path, delim=' '):
        self.path = path
        self.delim = delim
    
    def format(self, result, read):
        """Format a read and add it to `result`.
        """
        raise NotImplementedError()
    
    def _format(self, result, fields):
        result[self.path].append("".join((
            self.delim.join(str(f) for f in fields),
            "\n")))

class RestFormatter(DelimFormatter):
    """Rest file formatter.
    """
    def format(self, result, read):
        if read.match:
            rest = read.match.rest()
            if len(rest) > 0:
                self._format(result, (rest, read.name))

class InfoFormatter(DelimFormatter):
    """Info file formatter.
    """
    def __init__(self, path):
        super(InfoFormatter, self).__init__(path, delim='\t')
    
    def format(self, result, read):
        if read.match:
            for match_info in read.match_info:
                self._format(result, match_info[0:11])
        else:
            seq = read.sequence
            qualities = read.qualities if read.qualities is not None else ''
            self._format(result, (read.name, -1, seq, qualities))

class WildcardFormatter(DelimFormatter):
    """Wildcard file formatter.
    """
    def format(self, result, read):
        if read.match:
            self._format(result, (read.match.wildcards(), read.name))

def add_suffix_to_path(path, suffix):
    """
    Add the suffix (str or int) after the file name but
    before the extension.
    """
    name, ext1, ext2 = splitext_compressed(path)
    return "{}{}{}{}".format(name, suffix, ext1, ext2 or "")
