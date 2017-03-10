from atropos.trim.filters import NoFilter
from atropos.io import STDOUT, xopen, open_output

class Formatters(object):
    def __init__(self, output, seq_formatter_args):
        self.output = output
        self.multiplexed = output is not None and '{name}' in output
        self.seq_formatter_args = seq_formatter_args
        self.seq_formatters = {}
        self.mux_formatters = {}
        self.info_formatters = []
        self.discarded = 0
    
    def add_seq_formatter(self, filter_type, file1, file2=None):
        self.seq_formatters[filter_type] = create_seq_formatter(
            file1, file2, **self.seq_formatter_args)
    
    def add_info_formatter(self, formatter):
        self.info_formatters.append(formatter)
    
    def get_mux_formatter(self, name):
        assert self.multiplexed
        if name not in self.mux_formatters:
            path = self.output.format(name=name)
            self.mux_formatters[name] = create_seq_formatter(
                path, **self.seq_formatter_args)
        return self.mux_formatters[name]
    
    def get_seq_formatters(self):
        return (
            set(f for f in self.seq_formatters.values() if f.written > 0) |
            set(f for f in self.mux_formatters.values() if f.written > 0)
        )
    
    def summary(self):
        seq_formatters = self.get_seq_formatters()
        written = sum(f.written for f in seq_formatters)
        read_bp = [
            sum(f.read1_bp for f in seq_formatters),
            sum(f.read2_bp for f in seq_formatters)
        ]
        return (written, read_bp)
    
    def format(self, result, dest, read1, read2=None):
        if self.multiplexed and (dest == NoFilter) and read1.match:
            name = read1.match.adapter.name
            formatter = self.get_mux_formatter(name)
            formatter.format(result, read1, read2)
        elif dest in self.seq_formatters:
            self.seq_formatters[dest].format(result, read1, read2)
        else:
            self.discarded += 1
        
        for f in self.info_formatters:
            f.format(result, read1)
            if read2:
                f.format(result, read2)

class DelimFormatter(object):
    def __init__(self, path, delim=' '):
        self.path = path
        self.delim = delim
    
    def _format(self, result, fields):
        result[self.path].append("".join((
            self.delim.join(str(f) for f in fields),
            "\n"
        )))

class RestFormatter(DelimFormatter):
    def format(self, result, read):
        if read.match:
            rest = read.match.rest()
            if len(rest) > 0:
                self._format(result, (rest, read.name))

class InfoFormatter(DelimFormatter):
    def __init__(self, path):
        super(InfoFormatter, self).__init__(path, delim='\t')
    
    def format(self, result, read):
        if read.match:
            for m in read.match_info:
                self._format(result, m[0:11])
        else:
            seq = read.sequence
            qualities = read.qualities if read.qualities is not None else ''
            self._format(result, (read.name, -1, seq, qualities))

class WildcardFormatter(DelimFormatter):
    def format(self, result, read):
        if read.match:
            self._format(result, (read.match.wildcards(), read.name))

## Writing sequences to files ##

class Writers(object):
    def __init__(self, force_create):
        self.writers = {}
        self.force_create = force_create
        self.suffix = None
    
    def write_result(self, result, compressed=False):
        """
        Write results to file.
        
        result :: dict with keys being file descriptors
        and values being data (either bytes or strings).
        Strings are expected to already have appropriate
        line-endings.
        compressed :: whether data has already been compressed.
        """
        for file_desc, data in result.items():
            self.write(file_desc, data, compressed)
    
    def write(self, file_desc, data, compressed=False):
        """
        Write data to file.
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
        self.writers[path].write(data)
    
    def close(self):
        for path in self.force_create:
            if path not in self.writers and path != STDOUT:
                with open_output(path, "w"):
                    pass
        for writer in self.writers.values():
            if writer != sys.stdout:
                writer.close()
