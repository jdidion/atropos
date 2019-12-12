class AtroposError(Exception):
    """
    Base class for Atropos-specific errors.
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
