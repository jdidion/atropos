"""
Build atropos.

Cython is run when
* no pre-generated C sources are found,
* or the pre-generated C sources are out of date,
* or when --cython is given on the command line.
"""
import codecs
import os.path
import sys

from setuptools import setup, Extension, find_packages
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext

import versioneer

MIN_CYTHON_VERSION = "0.25.2"


if sys.version_info < (3, 3):
    sys.stdout.write("At least Python 3.3 is required.\n")
    sys.exit(1)


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
    """Exit if Cython was not found or is too old"""
    try:
        from Cython import __version__ as cyversion
    except ImportError:
        sys.stdout.write(
            "ERROR: Cython is not installed. Install at least Cython version " +
            str(MIN_CYTHON_VERSION) + " to continue.\n")
        sys.exit(1)
    if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
        sys.stdout.write(
            "ERROR: Your Cython is at version '{}' but at least version '{}' "
            "is required.\n".format(cyversion, MIN_CYTHON_VERSION))
        sys.exit(1)


extensions = [
    Extension("atropos.align._align", sources=["atropos/align/_align.pyx"]),
    Extension(
        "atropos.commands.trim._qualtrim",
        sources=["atropos/commands/trim/_qualtrim.pyx"]
    ),
    Extension("atropos.io._seqio", sources=["atropos/io/_seqio.pyx"]),
]

cmdclass = versioneer.get_cmdclass()
versioneer_build_ext = cmdclass.get("build_ext", _build_ext)
versioneer_sdist = cmdclass.get("sdist", _sdist)


class BuildExt(versioneer_build_ext):
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


cmdclass["build_ext"] = BuildExt


class SDist(versioneer_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        check_cython_version()
        cythonize(extensions)
        versioneer_sdist.run(self)


cmdclass["sdist"] = SDist


setup(
    name="atropos",
    version=versioneer.get_version(),
    cmdclass=cmdclass,
    author="John Didion",
    author_email="john.didion@nih.gov",
    url="https://atropos.readthedocs.org/",
    description="trim adapters from high-throughput sequencing reads",
    long_description=codecs.open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "README.md"
        ),
        "rb",
        "utf-8"
    ).read(),
    long_description_content_type="text/markdown",
    license="MIT",
    ext_modules=extensions,
    packages=find_packages(),
    scripts=["bin/atropos"],
    package_data={
        "atropos": [
            "adapters/*.fa",
            "commands/**/templates/*"
        ]
    },
    install_requires=["cython>={}".format(MIN_CYTHON_VERSION)],
    tests_require=["pytest"],  # , "jinja2", "pysam"],
    extras_require={
        "progressbar": ["progressbar2"],
        "tqdm": ["tqdm"],
        "khmer": ["khmer"],
        "pysam": ["pysam"],
        "jinja": ["jinja2"],
        "sra": ["srastream>=0.1.3"]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "License :: Public Domain",
        "Natural Language :: English",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6"
    ]
)
