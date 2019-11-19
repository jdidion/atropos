"""
Build atropos.

Cython is run when
* no pre-generated C sources are found,
* or the pre-generated C sources are out of date,
* or when --cython is given on the command line.
"""
from pathlib import Path
import sys

from setuptools import setup, Extension, find_packages
from distutils.version import LooseVersion
from distutils.command.sdist import sdist as _sdist
from distutils.command.build_ext import build_ext as _build_ext

import versioneer


MIN_CYTHON_VERSION = "0.25.2"


version_info = sys.version_info
if sys.version_info < (3, 6):
    sys.stdout.write("At least Python 3.6 is required.\n")
    sys.exit(1)


cmdclass = versioneer.get_cmdclass()
VersioneerBuildExt = cmdclass.get("build_ext", _build_ext)
VersioneerSdist = cmdclass.get("sdist", _sdist)


class BuildExtCython(VersioneerBuildExt):
    def run(self):
        # If we encounter a PKG-INFO file, then this is likely a .tar.gz/.zip
        # file retrieved from PyPI that already includes the pre-cythonized
        # extension modules, and then we do not need to run cythonize().
        if (Path.cwd() / "PKG-INFO").exists():
            no_cythonize(extensions)
        else:
            # Otherwise, this is a 'developer copy' of the code, and then the
            # only sensible thing is to require Cython to be installed.
            check_cython_version()
            from Cython.Build import cythonize
            self.extensions = cythonize(self.extensions)

        _build_ext.run(self)


cmdclass["build_ext"] = BuildExtCython


class SdistCython(VersioneerSdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        check_cython_version()
        cythonize(extensions)
        _sdist.run(self)


cmdclass["sdist"] = SdistCython


def out_of_date(_extensions):
    """
    Check whether any pyx source is newer than the corresponding generated
    C source or whether any C source is missing.
    """
    for extension in _extensions:
        for pyx in extension.sources:
            pyx_path = Path(pyx)
            if pyx_path.suffix not in (".pyx", ".py"):
                continue
            if extension.language == "c++":
                csource = pyx_path.with_suffix(".cpp")
            else:
                csource = pyx_path.with_suffix(".c")
            # When comparing modification times, allow five seconds slack:
            # If the installation is being run from pip, modification
            # times are not preserved and therefore depends on the order in
            # which files were unpacked.
            if not csource.exists() or (
                pyx_path.stat().st_mtime > csource.stat().st_mtime + 5
            ):
                return True

    return False


def check_cython_version():
    """Exit if Cython was not found or is too old"""
    try:
        from Cython import __version__ as cyversion
    except ImportError:
        sys.stdout.write(
            f"ERROR: Cython is not installed. Install at least Cython version "
            f"{MIN_CYTHON_VERSION} to continue.\n"
        )
        sys.exit(1)
    if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
        sys.stdout.write(
            f"ERROR: Your Cython is at version '{cyversion}', but at least "
            f"version {MIN_CYTHON_VERSION} is required.\n"
        )
        sys.exit(1)


def no_cythonize(_extensions, **_ignore):
    """
    Change file extensions from .pyx to .c or .cpp.

    Copied from Cython documentation
    """
    for extension in _extensions:
        sources = []
        for sfile in extension.sources:
            sfile_path = Path(sfile)
            if sfile_path.suffix in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = sfile_path.with_suffix(ext)
            sources.append(sfile)
        extension.sources[:] = list(str(src) for src in sources)


with open(
    Path(__file__).parent.absolute() / "README.md", encoding="utf-8"
) as f:
    readme = f.read()

# Define extensions to be Cythonized
extensions = [
    Extension("atropos.align._align", sources=["atropos/align/_align.pyx"]),
    Extension(
        "atropos.commands.trim._qualtrim",
        sources=["atropos/commands/trim/_qualtrim.pyx"],
    ),
    Extension("atropos.io._seqio", sources=["atropos/io/_seqio.pyx"]),
]

# TODO: load these from requirements*.txt file
install_requirements = [
    f"Cython>={MIN_CYTHON_VERSION}",
    "xphyle>=4.0.0-rc.0",
    "pokrok>=0.1.0"
]
test_requirements = [
    "pytest",  # 'jinja2', 'pysam',
]
extra_requirements = {
    "progressbar": ["progressbar2"],
    "tqdm": ["tqdm"],
    "khmer": ["khmer"],
    "pysam": ["pysam"],
    "jinja": ["jinja2"],
    "sra": ["ngstream"],
}

setup(
    name="atropos",
    version=versioneer.get_version(),
    cmdclass=cmdclass,
    author="John Didion",
    author_email="github@didion.net",
    url="https://atropos.readthedocs.org/",
    description="Specific, sensitive, and speedy trimming of NGS reads.",
    long_description=readme,
    long_description_content_type="text/markdown",
    license="MIT",
    ext_modules=extensions,
    packages=find_packages(),
    package_data={"atropos": ["adapters/*.fa", "commands/**/templates/*"]},
    install_requires=install_requirements,
    tests_require=test_requirements,
    extras_require=extra_requirements,
    entry_points={"console_scripts": ["atropos=atropos.__main__:main"]},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "License :: Public Domain",
        "Natural Language :: English",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
