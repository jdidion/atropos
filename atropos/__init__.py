# coding: utf-8
"""Top-level Atropos package.
"""
import sys

from atropos._version import get_versions

__version__ = get_versions()['version']
del get_versions


class AtroposError(Exception):
    """Base class for Atropos-specific errors.
    """
    pass


def check_importability() -> bool:  # pragma: no cover
    """Check that cython modules haev been compile.

    Returns:
        True if a cython-compiled module can be imported, otherwise False.
    """
    try:
        import atropos.align._align  # pylint: disable=unused-variable
        return True
    except ImportError as err:
        if 'undefined symbol' in str(err):
            return False
        else:
            raise
