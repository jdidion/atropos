class AtroposError(Exception):
    """
    Base class for Atropos-specific errors.
    """


class FormatError(AtroposError):
    """
    Raised when an input file (FASTA or FASTQ) is malformatted.
    """


class UnknownFileTypeError(AtroposError):
    """
    Raised when file type could not be detected.
    """


class MulticoreError(AtroposError):
    """
    Base error for parallel processes.
    """


class NotInAlphabetError(AtroposError):
    """
    Error that indicates an invalid character in an input sequence.
    """
    def __init__(self, character: str):
        super().__init__()
        self.character = character
