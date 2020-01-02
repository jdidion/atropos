from enum import IntFlag, Enum
from pathlib import Path
from typing import IO, Set, Union

from atropos.errors import UnknownFileTypeError
from atropos.utils.paths import splitext_compressed


class InputRead(IntFlag):
    SINGLE = 0
    READ1 = 1
    READ2 = 2
    PAIRED = 1 | 2


class SequenceFileType(Enum):
    FASTA = {".fasta", ".fa", ".fna", ".csfasta", ".csfa"}, True
    FASTQ = {".fastq", ".fq"}, True
    SAM = {".sam"}, True
    BAM = {".bam"}, False
    SRA_FASTQ = {}, False
    FASTA_QUAL = {}, False

    def __init__(self, exts: Set[str], can_output: bool):
        self.exts = exts
        self.can_output = can_output

    @staticmethod
    def guess_from_name(
        path: Union[str, Path, IO],
        output: bool,
        raise_on_failure: bool = False,
    ):
        """
            Detects file format based on the file name.

            Args:
                path: The filename to guess.
                output: Whether to only guess formats that can be output.
                raise_on_failure: Whether to raise an exception if the filename cannot
                    be detected.

            Returns:
                SequenceFileType.
            """
        name = None

        if isinstance(path, str):
            name = path
        elif hasattr(path, "name"):
            name = path.name

        if name:
            name, ext1, _ = splitext_compressed(name)
            ext = ext1.lower()

            for t in SequenceFileType:
                if output and not t.can_output:
                    continue
                if ext in t.exts:
                    return t

            if ext == ".txt" and name.endswith("_sequence"):
                return SequenceFileType.FASTQ

        if raise_on_failure:
            raise UnknownFileTypeError(
                f"Could not determine whether file {path!r} is FASTA, FASTQ, or SAM/BAM"
            )
