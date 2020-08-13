import math
from typing import Dict, Iterable, Optional, Sequence

from atropos.errors import NotInAlphabetError


class Alphabet:
    def __init__(
        self, name: str,
        valid_characters: Iterable[str],
        default_character: Optional[str] = None
    ):
        self.name = name
        self.valid_characters = set(valid_characters)
        self.default_character = default_character
        if self.default_character not in self.valid_characters:
            self.valid_characters.add(self.default_character)

    def __contains__(self, character: str) -> bool:
        return character in self.valid_characters

    def contains_all(self, string: str) -> bool:
        return set(string) <= self.valid_characters

    def validate(self, character: str) -> None:
        """
        Raises NotInAlphabetError if the character is not in the alphabet.
        """
        if character not in self:
            raise NotInAlphabetError(character)

    def validate_string(self, string: str) -> None:
        """
        Raises NotInAlphabetError if any character in 'string' is not in the alphabet.
        """
        diff = set(string) - self.valid_characters
        if diff:
            raise NotInAlphabetError("".join(diff))

    def resolve(self, character: str) -> str:
        """
        Returns 'character' if it's in the alphabet, otherwise the alphabet's default
        character.
        """
        if character in self.valid_characters:
            return character
        elif self.default_character:
            return self.default_character
        else:
            raise ValueError(
                f"Alphabet {self.name} does not contain character" f" {character}"
            )

    def resolve_string(self, string: str) -> str:
        """
        Returns a new string with any non-alphabet characters replaced with the
        default character.
        """
        return "".join(self.resolve(c) for c in string)


def build_iso_nucleotide_table() -> Dict[str, str]:
    """
    Generates a dict mapping ISO nucleotide characters to their complements, in both
    upper and lower case.

    TODO: the nucleotide table should be implemented as an alphabet.
    """
    nuc = {
        "A": "T",
        "C": "G",
        "R": "Y",
        "S": "S",
        "W": "W",
        "K": "M",
        "B": "V",
        "D": "H",
        "N": "N",
    }
    for base, comp in tuple(nuc.items()):
        nuc[comp] = base
        nuc[base.lower()] = comp.lower()
        nuc[comp.lower()] = base.lower()
    return nuc


BASE_COMPLEMENTS = build_iso_nucleotide_table()
DNA_BASES = frozenset("ACGT")
IUPAC_BASES = frozenset(("X",) + tuple(BASE_COMPLEMENTS.keys()))
"""Valid IUPAC bases, plus 'X'"""
GC_BASES = frozenset("CGRYSKMBDHVN")
"""IUPAC bases that include C or G."""
LOG2 = math.log(2)


ALPHABETS = {
    a.name: a
    for a in (
        Alphabet("dna", DNA_BASES, None),
        Alphabet("iupac", IUPAC_BASES, "N"),
        Alphabet("colorspace", "0123", None)
    )
}


def complement(seq: str) -> str:
    """
    Returns the complement of nucleotide sequence `seq`.
    """
    return "".join(BASE_COMPLEMENTS[base] for base in seq)


def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of nucleotide sequence `seq`.
    """
    return "".join(BASE_COMPLEMENTS[base] for base in reversed(seq))


def sequence_complexity(seq: str) -> float:
    """
    Computes a simple measure of sequence complexity.

    Args:
        seq: The sequence to measure.

    Returns:
        Complexity, as a value [0,2], where 0 = a homopolymer and
        2 = completely random.
    """
    seq = seq.upper()
    seqlen = float(len(seq))
    term = 0
    for base in ("A", "C", "G", "T"):
        count = seq.count(base)
        if count > 0:
            frac = count / seqlen
            term += frac * math.log(frac) / LOG2
    return -term


def qual2int(qual: str, base: int = 33) -> int:
    """
    Converts a quality charater to a phred-scale int.

    Args:
        qual: The quality value.
        base: The offset of the first quality value (Old Illumina = 64,
            new Illumina = 33).

    Returns:
        An iterable of integer qualities.
    """
    return ord(qual) - base


def quals2ints(quals: Sequence[str], base: int = 33) -> Sequence[int]:
    """
    Converts an iterable of quality characters to phred-scale ints.

    Args:
        quals: The qualities.
        base: The offset of the first quality value (Old Illumina = 64,
            new Illumina = 33).

    Returns:
        A tuple of integer qualities.
    """
    return tuple(ord(q) - base for q in quals)


def qual2prob(qchar: str) -> float:
    """
    Converts a quality char to a probability.
    """
    return 10 ** (-qual2int(qchar) / 10)
