# coding: utf-8
"""Top-level Atropos package.
"""
import os
import sys

try:
    from importlib.metadata import version as _get_version
    __version__ = _get_version("atropos")
except Exception:
    __version__ = "0.0.0"

class AtroposError(Exception):
    """Base class for Atropos-specific errors.
    """
    pass

def check_importability():  # pragma: no cover
    """Check that cython modules haev been compile.
    """
    try:
        import atropos.align._align # pylint: disable=unused-variable
    except ImportError as err:
        if 'undefined symbol' in str(err):
            print("""
ERROR: A required extension module could not be imported because it is
incompatible with your system. A quick fix is to recompile the extension
modules with the following command:

    {0} setup.py build_ext -i

See the documentation for alternative ways of installing the program.

The original error message follows.
""".format(sys.executable))
        raise
