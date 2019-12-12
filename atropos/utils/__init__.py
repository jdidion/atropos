from enum import Enum, IntEnum
import errno
import os
from pathlib import Path
import sys
from typing import Callable, Iterable, Iterator, Optional, Tuple, TypeVar

from loguru import logger
import pokrok as pk
from xphyle import STDOUT, STDERR
from xphyle.formats import FORMATS

from atropos.errors import AtroposError


T = TypeVar("T")


class LoggingConfig:
    """
    Sets up logging and retain state of whether setup has already been done.
    Specifically, don't do setup again if this function is being called externally
    such as from unit tests.
    """
    def __init__(self) -> None:
        self._was_setup = False

    @property
    def was_setup(self) -> bool:
        return self._was_setup

    def setup(
        self,
        log_path: Path,
        log_level: str,
        fmt: str = "%(asctime)s %(levelname)s: %(message)s",
        force: bool = False
    ):
        if self._was_setup and not force:
            raise ValueError("Logging was already setup")

        if log_path == STDOUT:
            sink = sys.stdout
        elif log_path == STDERR:
            sink = sys.stderr
        else:
            sink = log_path

        logger.add(sink=sink, level=log_level, format=fmt)

        self._was_setup = True


LOGGING_CONFIG = LoggingConfig()


class ReturnCode(IntEnum):
    SUCCESS = 0
    ERROR = 1
    TERMINATED = 130


class classproperty(classmethod):
    """
    Decorator for read-only class property.
    """
    def __init__(self, f):
        super().__init__(f)
        self.f = f

    def __get__(self, *args, **kwargs):
        _, owner = args
        return self.f(owner)


class Magnitude(Enum):
    """Enumeration of some common number scales.
    """

    G = ("billion", "giga", 1e9)
    M = ("million", "mega", 1e6)
    K = ("thousand", "kilo", 1e3)

    def __init__(self, label, prefix, divisor):
        self.label = label
        self.prefix = prefix
        self.divisor = divisor


def enumerate_range(
    collection: Iterable[T], start, end, step=1
) -> Iterable[Tuple[int, T]]:
    """
    Generates an indexed series:  (0,coll[0]), (1,coll[1]) ...

    Only up to (start-end+1) items will be yielded.

    Args:
        collection: The collection to enumerate.
        start: The starting index.
        end: The ending index.
        step: The amount to increment the index each iteration (defaults to 1).

    Yields:
        (index, item) tuples.
    """
    idx = start
    itr = iter(collection)
    while idx < end:
        yield idx, next(itr)
        idx += step


def truncate_string(string: Optional[str], max_len: int = 100) -> Optional[str]:
    """
    Shortens string to at most `max_len` characters, appending "..." if necessary.
    """
    if string is None:
        return None
    if len(string) > max_len:
        string = string[: max_len - 3] + "..."
    return string


def splitext_compressed(name: str) -> Tuple[str, str, str]:
    """
    Splits the filename and extensions of a file that potentially has two extensions -
    one for the file type (e.g. 'fq') and one for the compression type (e.g. 'gz').

    Args:
        name: The filename.

    Returns:
        A tuple (name, ext1, ext2), where ext1 is the filetype extension and
        ext2 is the compression type extension, or None.
    """
    ext2 = None
    for ext in FORMATS.list_compression_formats(with_sep=True):
        if name.endswith(ext):
            ext2 = ext
            name = name[: -len(ext)]
            break

    name, ext1 = os.path.splitext(name)
    return name, ext1, ext2


def run_interruptible(func: Callable, *args, **kwargs) -> ReturnCode:
    """
    Runs a function, gracefully handling keyboard interrupts.

    Args:
        func: The function to execute.
        args, kwargs: Positional and keyword arguments to pass to `func`.

    Returns:
        A (unix-style) return code (0=normal, anything else is an error).
    """
    retcode = ReturnCode.SUCCESS
    try:
        func(*args, **kwargs)
    except KeyboardInterrupt:
        logger.error("Interrupted")
        retcode = ReturnCode.TERMINATED
    except IOError as err:
        if err.errno == errno.EPIPE:
            retcode = ReturnCode.ERROR
        else:
            raise
    except (AtroposError, EOFError):
        logger.exception("Atropos error")
        retcode = ReturnCode.ERROR
    except:  # pylint: disable=broad-except
        logger.exception("Unknown error")
        retcode = ReturnCode.ERROR
    return retcode


def create_progress_reader(
    reader: Iterable[T],
    progress_type: str = "bar",
    batch_size: int = 1,
    max_items: Optional[int] = None,
    **kwargs
) -> Iterator[T]:
    """
    Wraps an iterable in a progress bar of the specified type.

    Args:
        reader: The iterable to wrap.
        progress_type: msg = a custom progress bar that reports via log
            messages; bar = use a ProgressBar (from the progressbar library)
            or tqdm.
        max_items: Max number of items, if known in advance.
        batch_size: The number of records in each iterable item (iterable is
            typically a BatchReader).
        kwargs: Additional arguments to pass to the progress bar constructor.

    Returns:
        A wrapped iterable. If `progress_type == 'bar'` and neither of the
        supported libraries are available, a warning is logged and the unwrapped
        reader is returned.
    """
    if progress_type == "msg":
        pk.set_plugins(["logging"])

    if max_items:
        widgets = [
            pk.Widget.COUNTER,
            pk.Widget.PERCENT,
            pk.Widget.ELAPSED,
            pk.Widget.BAR,
            pk.Widget.ETA,
        ]
    else:
        widgets = [pk.Widget.COUNTER, pk.Widget.ELAPSED, pk.Widget.SPINNER]

    return pk.progress_iter(
        reader,
        desc="Processed",
        size=max_items,
        unit="records",
        multiplier=batch_size,
        style=pk.Style(widgets),
        **kwargs,
    )
