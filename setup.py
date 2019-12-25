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
from distutils.command.sdist import sdist
from distutils.command.build_ext import build_ext

version_info = sys.version_info
if sys.version_info < (3, 6):
    sys.stdout.write("At least Python 3.6 is required.\n")
    sys.exit(1)


MIN_CYTHON_VERSIONS = {
    (3, 6): "0.25.2",
    (3, 7): "0.29",
    (3, 8): "0.29.14"
}
min_cython_version = MIN_CYTHON_VERSIONS[version_info[0:2]]


class BuildExtCython(build_ext):
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

        super().run()


class SdistCython(sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        check_cython_version()
        cythonize(extensions)
        super().run()


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
            f"{min_cython_version} to continue.\n"
        )
        sys.exit(1)
    if LooseVersion(cyversion) < LooseVersion(min_cython_version):
        sys.stdout.write(
            f"ERROR: Your Cython is at version '{cyversion}', but at least "
            f"version {min_cython_version} is required.\n"
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
    Extension("atropos.aligner._aligner", sources=["atropos/aligner/_aligner.pyx"]),
    Extension(
        "atropos.commands.trim._qualtrim",
        sources=["atropos/commands/trim/_qualtrim.pyx"],
    ),
    Extension(
        "atropos.io.sequence._sequence",
        sources=["atropos/io/sequence/_sequence.pyx"]
    ),
    Extension(
        "atropos.io.readers._readers",
        sources=["atropos/io/readers/_readers.pyx"]
    ),
]

install_requirements = [
    f"Cython>={min_cython_version}",
    "loguru>=0.4.0",
    "pokrok>=0.2.0",
    "xphyle>=4.2.1",
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
    "sra": ["ngstream>=0.2.2"],
}

setup(
    name="atropos",
    use_scm_version=True,
    cmdclass={
        "build_ext": BuildExtCython,
        "sdist": SdistCython
    },
    author="John Didion",
    author_email="github@didion.net",
    url="https://atropos.readthedocs.org/",
    description="Specific, sensitive, and speedy trimming of NGS reads.",
    long_description=readme,
    long_description_content_type="text/markdown",
    license="MIT",
    ext_modules=extensions,
    packages=find_packages(),
    package_data={
        "atropos": ["adapters/*.fa", "commands/**/templates/*"]
    },
    setup_requires=["setuptools_scm"],
    install_requires=install_requirements,
    tests_require=test_requirements,
    extras_require=extra_requirements,
    entry_points={
        "console_scripts": [
            "atropos=atropos.__main__:main"
        ],
        "atropos.commands": [
            "detect = atropos.commands.detect.console:DetectCommandConsole",
            "error = atropos.commands.error.console:ErrorCommandConsole",
            "qc = atropos.commands.qc.console:QcCommandConsole",
            "trim = atropos.commands.trim.console:TrimCommandConsole",
        ]
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
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        #"Programming Language :: Python :: 3.8",
    ],
)
