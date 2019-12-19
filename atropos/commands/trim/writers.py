from abc import ABCMeta, abstractmethod
from pathlib import Path
import sys
from typing import (
    Any, Dict, IO, List, Optional, Sequence as SequenceType, Set, Tuple, Type, Union
)

from xphyle import STDOUT, xopen, open_
from xphyle.types import ModeArg

from atropos.commands.trim.filters import Filter, NoFilter
from atropos.io.sequence import Sequence
from atropos.io.formatters import Formatter, create_seq_formatter
from atropos.utils.paths import splitext_compressed
from atropos.utils.collections import Summarizable


FileDesc = Union[Path, Tuple[Path, ModeArg]]


class Writers:
    """
    Manages writing to one or more outputs.
    """

    def __init__(
        self,
        force_create: SequenceType[Path] = (),
        compression_format: Optional[str] = None
    ):
        """
        Args:
            force_create: Whether empty output files should be created.
        """
        self.writers = {}
        self.force_create = force_create
        self.suffix = None
        self.compression_format = compression_format

    def get_writer(self, file_desc: FileDesc, compressed: bool = False) -> IO:
        """
        Create the writer for a file descriptor if it does not already exist.

        Args:
            file_desc: File descriptor. If `compressed==True`, this is a tuple
                (path, mode), otherwise it's only a path.
            compressed: Whether data has already been compressed.

        Returns:
            The writer.
        """
        if compressed:
            path, mode = file_desc
            compression = False
        else:
            path = file_desc
            mode = "w"
            compression = self.compression_format

        if path not in self.writers:
            if self.suffix:
                real_path = add_suffix_to_path(path, self.suffix)
            else:
                real_path = path

            # TODO: test whether O_NONBLOCK allows non-blocking write to NFS
            self.writers[path] = xopen(real_path, mode, compression=compression)

        return self.writers[path]

    def write_result(self, result, compressed: bool = False) -> None:
        """
        Writes results to output.

        Args:
            result: Dict with keys being file descriptors and values being data
                (either bytes or strings). Strings are expected to already have
                appropriate line-endings.
            compressed: Whether data has already been compressed.
        """
        for file_desc, data in result.items():
            self.write(file_desc, data, compressed)

    def write(self, file_desc: FileDesc, data, compressed: bool = False):
        """
        Writes data to output. If the specified path has not been seen before,
        the output is opened.

        Args:
            file_desc: File descriptor. If `compressed==True`, this is a tuple
                (path, mode), otherwise it's only a path.
            data: The data to write.
            compressed: Whether data has already been compressed.
        """
        self.get_writer(file_desc, compressed).write(data)

    def close(self):
        """
        Closes all outputs.
        """
        for path in self.force_create:
            if path not in self.writers and path != STDOUT:
                with open_(path, "w"):
                    pass

        for writer in self.writers.values():
            if writer not in (sys.stdout, sys.stderr):
                writer.close()


class DetailFormatter(metaclass=ABCMeta):
    """
    Base class for formatters that write to a delimited file.

    Args:
        path: The output file path.
        delim: The field delimiter.
    """

    def __init__(self, path, delim=" "):
        self.path = path
        self.delim = delim

    @abstractmethod
    def format(self, result: dict, read: Sequence):
        """
        Format a read and add it to `result`.
        """

    def _format(self, result: dict, fields: SequenceType[Any]):
        if self.path not in result:
            result[self.path] = []

        result[self.path].append(
            "".join((self.delim.join(str(f) for f in fields), "\n"))
        )


class RestFormatter(DetailFormatter):
    """
    Rest file formatter.
    """

    def format(self, result: dict, read: Sequence):
        if read.match:
            rest = read.match.rest()

            if len(rest) > 0:
                self._format(result, (rest, read.name))


class InfoFormatter(DetailFormatter):
    """
    Info file formatter.
    """

    def __init__(self, path: Path):
        super().__init__(path, delim="\t")

    def format(self, result: dict, read: Sequence):
        if read.match:
            for match_info in read.match_info:
                self._format(result, match_info[0:11])
        else:
            seq = read.sequence
            qualities = read.qualities if read.qualities is not None else ""
            self._format(result, (read.name, -1, seq, qualities))


class WildcardFormatter(DetailFormatter):
    """
    Wildcard file formatter.
    """

    def format(self, result: dict, read: Sequence):
        if read.match:
            self._format(result, (read.match.wildcards(), read.name))


class Formatters(Summarizable):
    """
    Manages multiple formatters.
    """

    def __init__(self, output: str, seq_formatter_args: dict):
        """
        Args:
            output: The output file name template.
            seq_formatter_args: Additional arguments to pass to the formatter
                constructor.
        """
        self.output = output
        self.multiplexed = output is not None and "{name}" in output.name
        self.seq_formatter_args = seq_formatter_args
        self.seq_formatters: Dict[Type[Filter], Formatter] = {}
        self.mux_formatters: Dict[str, Formatter] = {}
        self.detail_formatters: List[DetailFormatter] = []
        self.discarded = 0

    def add_seq_formatter(
        self, filter_type: Type[Filter], file1: Path, file2: Optional[Path] = None
    ):
        """
        Adds a formatter.

        Args:
            filter_type: The type of filter that triggers writing with the
                formatter.
            file1:
            file2: The output file(s).
        """
        self.seq_formatters[filter_type] = create_seq_formatter(
            file1, file2, **self.seq_formatter_args
        )

    def add_detail_formatter(self, formatter: DetailFormatter):
        """
        Add a formatter for one of the delimited detail files (rest, info, wildcard).
        """
        self.detail_formatters.append(formatter)

    def get_mux_formatter(self, name: str) -> Formatter:
        """
        Returns the formatter associated with the given name (barcode) when running
        in multiplexed mode.
        """
        if not self.multiplexed:
            raise ValueError(
                "Cannot call 'get_mux_formatter' for non-multiplexed data"
            )

        if name not in self.mux_formatters:
            path = Path(str(self.output).format(name=name))
            self.mux_formatters[name] = create_seq_formatter(
                path, **self.seq_formatter_args
            )

        return self.mux_formatters[name]

    def get_seq_formatters(self) -> Set[Formatter]:
        """
        Returns a set containing all formatters that have handled at least one record.
        """
        return set(f for f in self.seq_formatters.values() if f.written > 0) | set(
            f for f in self.mux_formatters.values() if f.written > 0
        )

    def format(
        self,
        result: dict,
        dest: Type[Filter],
        read1: Sequence,
        read2: Optional[Sequence] = None
    ):
        """
        Formats read(s) and add to a result dict. Also writes info records to any
        registered info formatters.

        Args:
            result: The result dict.
            dest: The destination (filter type).
            read1: The first read.
            read2: The second read.
        """
        if self.multiplexed and (dest == NoFilter) and read1.match:
            name = read1.match.adapter.name
            formatter = self.get_mux_formatter(name)
            formatter.format(result, read1, read2)
        elif dest in self.seq_formatters:
            self.seq_formatters[dest].format(result, read1, read2)
        else:
            self.discarded += 1

        for fmtr in self.detail_formatters:
            fmtr.format(result, read1)
            if read2:
                fmtr.format(result, read2)

    def summarize(self) -> dict:
        """
        Returns a summary dict.
        """
        seq_formatters = self.get_seq_formatters()
        return dict(
            records_written=sum(f.written for f in seq_formatters),
            bp_written=[
                sum(f.read1_bp for f in seq_formatters),
                sum(f.read2_bp for f in seq_formatters),
            ],
        )


def add_suffix_to_path(path: Path, suffix: str) -> str:
    """
    Add the suffix (str or int) after the file name but
    before the extension.
    """
    name, ext1, ext2 = splitext_compressed(path)
    return f"{name}{suffix}{ext1}{ext2 or ''}"
