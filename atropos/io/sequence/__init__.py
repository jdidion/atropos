from typing import Optional

from atropos.errors import FormatError
from atropos.utils import truncate_string
from atropos.utils.ngs import Alphabet

from ._sequence import Sequence


class ColorspaceSequence(Sequence):
    """
    Sequence object for colorspace reads.
    """

    def __init__(
        self,
        name: str,
        sequence: str,
        qualities: Optional[str] = None,
        primer: Optional[str] = None,
        name2: str = "",
        original_length: Optional[int] = None,
        match=None,
        match_info=None,
        clipped=None,
        insert_overlap: bool = False,
        merged: bool = False,
        corrected: int = 0,
        umi: Optional[str] = None,
        alphabet: Optional[Alphabet] = None,
    ):
        # In colorspace, the first character is the last nucleotide of the
        # primer base and the second character encodes the transition from the
        # primer base to the first real base of the read.
        if primer is None:
            self.primer = sequence[0:1]
            sequence = sequence[1:]
        else:
            self.primer = primer

        if self.primer not in ("A", "C", "G", "T"):
            raise FormatError(
                f"Primer base is {self.primer!r} in read {truncate_string(name)!r}, "
                f"but it should be one of A, C, G, T."
            )

        if qualities is not None and len(sequence) != len(qualities):
            rname = truncate_string(name)
            raise FormatError(
                f"In read named {rname!r}: length of colorspace quality sequence "
                f"({len(qualities)}) and length of read ({len(sequence)}) do not match "
                f"(primer is: {self.primer!r})"
            )

        super().__init__(
            name,
            sequence,
            qualities,
            name2,
            original_length,
            match,
            match_info,
            clipped,
            insert_overlap,
            merged,
            corrected,
            umi,
            alphabet,
        )

    def __repr__(self):
        fmt_str = "<ColorspaceSequence(name={0!r}, primer={1!r}, sequence={2!r}{3})>"

        if self.qualities is None:
            qstr = ""
        else:
            qstr = f", qualities={truncate_string(self.qualities)!r}"

        return fmt_str.format(
            truncate_string(self.name),
            self.primer,
            truncate_string(self.sequence),
            qstr,
        )

    def __getitem__(self, key):
        return self.__class__(
            self.name,
            self.sequence[key],
            self.qualities[key] if self.qualities is not None else None,
            self.primer,
            self.name2,
            self.original_length,
            self.match,
            self.match_info,
            self.clipped,
            self.insert_overlap,
            self.merged,
            self.corrected,
            self.umi,
        )
