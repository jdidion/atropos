from enum import IntFlag
from pathlib import Path
from typing import IO, Union

from atropos.errors import UnknownFileTypeError
from atropos.utils.paths import splitext_compressed


class InputRead(IntFlag):
    SINGLE = 0
    READ1 = 1
    READ2 = 2
    PAIRED = 1 | 2


def guess_format_from_name(path: Union[str, Path, IO], raise_on_failure: bool = False):
    """
    Detects file format based on the file name.

    Args:
        path: The filename to guess.
        raise_on_failure: Whether to raise an exception if the filename cannot
            be detected.

    Returns:
        The format name.
    """
    name = None

    if isinstance(path, str):
        name = path
    elif hasattr(path, "name"):  # seems to be an open file-like object
        name = path.name

    if name:
        name, ext1, _ = splitext_compressed(name)
        ext = ext1.lower()
        if ext in (".fasta", ".fa", ".fna", ".csfasta", ".csfa"):
            return "fasta"
        elif ext in (".fastq", ".fq") or (ext == ".txt" and name.endswith("_sequence")):
            return "fastq"
        elif ext in (".sam", ".bam"):
            return ext[1:]

    if raise_on_failure:
        raise UnknownFileTypeError(
            "Could not determine whether file {path!r} is FASTA or FASTQ: file "
            "name extension {ext!r} not recognized"
        )
