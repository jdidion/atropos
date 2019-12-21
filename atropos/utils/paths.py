import os
from pathlib import Path
from typing import Optional, Tuple, Union

from xphyle import FORMATS
from xphyle.paths import as_path, check_path
from xphyle.types import ModeAccess, PathLike, PathType


def splitext_compressed(name: Union[str, Path]) -> Tuple[str, str, str]:
    """
    Splits the filename and extensions of a file that potentially has two extensions -
    one for the file type (e.g. 'fq') and one for the compression type (e.g. 'gz').

    Args:
        name: The filename.

    Returns:
        A tuple (name, ext1, ext2), where ext1 is the filetype extension and
        ext2 is the compression type extension, or None.
    """
    name_str = str(name)
    ext2 = None

    for ext in FORMATS.list_extensions(with_sep=True):
        if name_str.endswith(ext):
            ext2 = ext
            name_str = name_str[:-len(ext)]
            break

    name_str, ext1 = os.path.splitext(name_str)

    return name_str, ext1, ext2


def as_readable_path(path: Union[str, PathLike], exists: bool = True) -> Optional[Path]:
    if path is None:
        if exists:
            raise ValueError("'path' is None")
        else:
            return None

    p = as_path(path)

    if p.exists():
        check_path(p, PathType.FILE, ModeAccess.READ)
    elif exists:
        raise ValueError(f"Path {p} does not exist")

    return p
