"""
Build Cython extensions for atropos.

This setup.py is retained solely for building C extensions from .pyx files.
All project metadata lives in pyproject.toml.

Cython is run when:
* no pre-generated C sources are found,
* or the pre-generated C sources are out of date,
* or when --cython is given on the command line.
"""
import os

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist


MIN_CYTHON_VERSION = "0.25.2"


def out_of_date(_extensions):
    """
    Check whether any pyx source is newer than the corresponding generated
    C source or whether any C source is missing.
    """
    for extension in _extensions:
        for pyx in extension.sources:
            path, ext = os.path.splitext(pyx)
            if ext not in (".pyx", ".py"):
                continue
            if extension.language == "c++":
                csource = path + ".cpp"
            else:
                csource = path + ".c"
            # When comparing modification times, allow five seconds slack:
            # If the installation is being run from pip, modification
            # times are not preserved and therefore depends on the order in
            # which files were unpacked.
            if not os.path.exists(csource) or (
                os.path.getmtime(pyx) > os.path.getmtime(csource) + 5
            ):
                return True
    return False


def no_cythonize(_extensions, **_ignore):
    """
    Change file extensions from .pyx to .c or .cpp.

    Copied from Cython documentation
    """
    for extension in _extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources


def check_cython_version():
    """Exit if Cython was not found or is too old."""
    from packaging.version import Version

    try:
        from Cython import __version__ as cyversion
    except ImportError:
        raise RuntimeError(
            "Cython is not installed. Install at least Cython version "
            + str(MIN_CYTHON_VERSION)
            + " to continue."
        )
    if Version(cyversion) < Version(MIN_CYTHON_VERSION):
        raise RuntimeError(
            "Your Cython is at version '{}' but at least version '{}' "
            "is required.".format(cyversion, MIN_CYTHON_VERSION)
        )


extensions = [
    Extension("atropos.align._align", sources=["atropos/align/_align.pyx"]),
    Extension(
        "atropos.commands.trim._qualtrim",
        sources=["atropos/commands/trim/_qualtrim.pyx"],
    ),
    Extension("atropos.io._seqio", sources=["atropos/io/_seqio.pyx"]),
]


class BuildExt(_build_ext):
    def run(self):
        # If we encounter a PKG-INFO file, then this is likely a .tar.gz/.zip
        # file retrieved from PyPI that already includes the pre-cythonized
        # extension modules, and then we do not need to run cythonize().
        if os.path.exists("PKG-INFO"):
            no_cythonize(extensions)
        else:
            # Otherwise, this is a "developer copy" of the code, and then the
            # only sensible thing is to require Cython to be installed.
            check_cython_version()
            from Cython.Build import cythonize

            self.extensions = cythonize(self.extensions)
        _build_ext.run(self)


class SDist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize

        check_cython_version()
        cythonize(extensions)
        _sdist.run(self)


setup(
    ext_modules=extensions,
    cmdclass={"build_ext": BuildExt, "sdist": SDist},
)
