# coding: utf-8
import os
import sys
from atropos.io import xopen, open_output
from atropos.io.compression import get_compressor

base = "tests/data/small.fastq"
files = [base + ext for ext in ["", ".gz", ".bz2", ".xz"]]


def test_context_manager():
    major, minor = sys.version_info[0:2]
    for name in files:
        if major == 2 and minor == 6:
            continue  # Py26 compression libraries do not support context manager protocol.

        with xopen(name, "rt") as f:
            lines = list(f)
            assert len(lines) == 12
            assert lines[5] == "AGCCGCTANGACGGGTTGGCCCTTAGACGTATCT\n", name


def test_append(tmp_path):
    for ext in ["", ".gz"]:  # BZ2 does NOT support append
        text = "AB"
        reference = text + text
        filename = "truncated.fastq" + ext
        mode = "a"
        if ext != "":
            mode = "ab"
            text = text.encode()
            reference = text + text
            text = get_compressor(filename).compress(
                text
            )  # On Py3, need to send BYTES, not unicode
        print("Trying ext=%s" % ext)
        path = tmp_path / filename
        try:
            os.unlink(path)
        except OSError:
            pass
        with open_output(str(path), mode) as f:
            f.write(text)
        print(path)
        with open_output(str(path), mode) as f:
            f.write(text)
        with xopen(str(path), "r") as f:
            try:
                reference = reference.decode("utf-8")
            except AttributeError:
                pass
            for appended in f:
                assert appended == reference


def test_xopen_text():
    for name in files:
        f = None
        try:
            f = xopen(name, "rt")
            lines = list(f)
            assert len(lines) == 12
            assert lines[5] == "AGCCGCTANGACGGGTTGGCCCTTAGACGTATCT\n", name
        finally:
            if f is not None:
                f.close()


def test_xopen_binary():
    for name in files:
        f = None
        try:
            f = xopen(name, "rb")
            lines = list(f)
            assert len(lines) == 12
            assert lines[5] == b"AGCCGCTANGACGGGTTGGCCCTTAGACGTATCT\n", name
        finally:
            if f is not None:
                f.close()
