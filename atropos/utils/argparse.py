from argparse import ArgumentParser, ArgumentTypeError, HelpFormatter, Namespace
from numbers import Number
import operator
from pathlib import Path
import re
import textwrap
from typing import Callable, Sequence, List, Type, TypeVar, Generic, Union, Optional
from urllib.parse import urlparse

from xphyle.paths import (
    STDOUT,
    STDERR,
    check_path,
    as_path,
    resolve_path,
    check_readable_file,
    check_writable_file,
)
from xphyle.types import PathLike, PathType, Permission

from atropos.utils import Magnitude


class AtroposArgumentParser(ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._groups = {}

    def add_group(
        self,
        name: str,
        title: Optional[str] = None,
        description: Optional[str] = None,
        mutex: bool = False,
        required: bool = False,
    ):
        """
        Add a group to the parser. The group will be stored under `name` and can
        later be retrieved via `get_group`.
        """
        if name in self._groups:
            raise ValueError(f"Group already exists: {name}")

        if mutex:
            group = self.add_mutually_exclusive_group(required=required)
        else:
            group = self.add_argument_group(title or name, description)
        self._groups[name] = group
        return group

    def get_group(self, name: str):
        """
        If a group has already been created with `name`, return the group, otherwise
        create a new group with that name.
        """
        if name in self._groups:
            return self._groups[name]
        else:
            return self.add_group(name)


class ParagraphHelpFormatter(HelpFormatter):
    """
    HelpFormatter that wraps text in paragraphs.
    """

    def _fill_text(self, text, width, indent):
        text = re.sub("[ \t]{2,}", " ", text)
        paragraphs = [
            textwrap.fill(p, width, initial_indent=indent, subsequent_indent=indent)
            for p in re.split("\n\n", text)
        ]
        return "\n\n".join(paragraphs)


T = TypeVar("T")


class CompositeType:
    """A composite of multiple data types.
    """

    def __init__(self, *types):
        self.types = types

    def __call__(self, string):
        result = string
        for datatype in self.types:
            result = datatype(result)
        return result


class ComparisonValidator(Generic[T]):
    """
    Validator that compares an argument against an expected value.
    """

    def __init__(self, rhs: T, oper: Callable[[T, T], bool], expected: bool = True):
        self.rhs = rhs
        self.oper = oper
        self.expected = expected

    def __call__(self, lhs: T):
        if self.oper(lhs, self.rhs) != self.expected:
            raise ArgumentTypeError(
                f"{self.oper}({lhs}, {self.rhs}) != {self.expected}"
            )

        return lhs


class CharList:
    """
    Parses a string into a list of characters and ensures they are all in the choices
    tuple.
    """

    def __init__(self, choices):
        self.choices = set(choices)

    def __call__(self, string):
        chars = list(string)

        if not all(char in self.choices for char in chars):
            raise ArgumentTypeError(
                f"One or more characters in {chars} not in {self.choices}"
            )

        return chars


class Delimited(Generic[T]):
    """
    Splits a string argument using a delimiter.
    """

    def __init__(
        self,
        delim: str = ",",
        data_type: Optional[Callable[[str], T]] = None,
        choices: Optional[Sequence[Union[str, T]]] = None,
        min_len: Optional[int] = None,
        max_len: Optional[int] = None,
    ):
        self.delim = delim
        self.data_type = data_type
        self.choices = choices
        self.min_len = min_len
        self.max_len = max_len

    def __call__(self, arg: Union[str, List[T]]) -> List[T]:
        if isinstance(arg, str):
            vals = arg.split(self.delim) if self.delim else [arg]
        else:
            vals = arg
        if vals[0] == "*" and self.choices is not None:
            vals = self.choices
        if self.data_type:
            vals = [self.data_type(v) for v in vals]

        if self.min_len and len(vals) < self.min_len:
            raise ArgumentTypeError(f"There must be at least {self.min_len} values")

        if self.max_len and len(vals) > self.max_len:
            raise ArgumentTypeError(f"There can be at most {self.max_len} values")

        return vals


class EnumChoice(Generic[T]):
    def __init__(self, enum_class: Callable[[str], T] = None):
        self.enum_class = enum_class

    def __call__(self, arg: str) -> T:
        return self.enum_class(arg)


class AccessiblePath:
    """
    Tests that a path is accessible.
    """

    def __init__(self, path_type: PathType, mode: Permission):
        self.path_type = path_type
        self.mode = mode

    def __call__(self, path: Union[str, PathLike]) -> Path:
        path = as_path(path)

        if self.path_type == PathType.FILE and path in (STDOUT, STDERR):
            return path

        try:
            return Path(check_path(path, self.path_type, self.mode))
        except IOError:
            raise ArgumentTypeError(
                f"Path {path} is not an accessible {self.path_type} with mode "
                f"{self.mode}"
            )


def existing_path(path: Union[str, PathLike]) -> Path:
    """
    Tests that a path exists.
    """
    path = as_path(path)

    if path == STDOUT:
        return path

    try:
        return Path(resolve_path(path))
    except IOError:
        raise ArgumentTypeError(f"Path {path} does not exist")


def readable_file(path):
    try:
        return check_readable_file(resolve_path(path))
    except IOError:
        raise ArgumentTypeError(
            f"Path {path} does not exist, is not a file, or is not readable"
        )


def writable_file(path):
    try:
        return check_writable_file(as_path(path, "w"))
    except IOError:
        raise ArgumentTypeError(f"Path {path} is not writable")


def readwritable_file(path):
    """Tests that a file is both readable and writeable."""
    path = as_path(path)

    if path.exists():
        path = readable_file(path)

    return writable_file(path)


def readable_url(url):
    """
    Validator for a URL that must be readable.

    Args:
        url: The URL to validate.
    """
    parsed = urlparse(url)
    scheme = parsed.scheme or "file"
    if scheme == "file":
        filename = readable_file(parsed.path)
        return "file:" + str(filename)
    else:
        return url


str_list = Delimited[str](data_type=str)
"""Comma-delimited list of strings."""

INT_OR_STR_RE = re.compile(r"([\d.]+)([KkMmGg]?)")


def int_or_str(arg):
    """
    Similar to int(), but accepts K, M, and G abbreviations.
    """
    if arg is None or isinstance(arg, int):
        return arg
    elif isinstance(arg, str):
        match = INT_OR_STR_RE.match(arg.upper())
        num, mult = match.groups()
        if mult:
            return int(float(num) * Magnitude[mult].divisor)
        else:
            return int(num)
    else:
        raise ArgumentTypeError(f"Unsupported type {arg}")


def positive(type_=int, inclusive=False):
    """
    Tests that a number is greater than (or equal to, if ``inclusive=True``) zero.
    """
    oper = operator.ge if inclusive else operator.gt
    return CompositeType(type_, ComparisonValidator(0, oper))


N = TypeVar("N", bound=Number)


def between(min_val: N = None, max_val: N = None, type_: Type[N] = int):
    """
    Returns a CompositeType that validates `min_val <= x <= max_val`.
    """
    return CompositeType(
        type_,
        ComparisonValidator(min_val, operator.ge),
        ComparisonValidator(max_val, operator.le),
    )


probability = between(0, 1, float)
"""A float between 0-1 (inclusive)."""
