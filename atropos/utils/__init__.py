from enum import Enum, IntEnum
import errno
from importlib import import_module
from pathlib import Path
import sys
from typing import Callable, Iterable, Optional, Tuple, TypeVar

from loguru import logger
from xphyle import STDOUT, STDERR

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

    def reset(self):
        if self._was_setup:
            logger.remove()
            self._was_setup = False
            return True
        else:
            return False


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


def no_import(lib):
    """
    Tests whether a library is importable.
    """
    try:
        mod = import_module(lib)
        return mod is None
    except ImportError:
        return True
