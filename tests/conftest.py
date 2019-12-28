from pathlib import Path
import traceback
from typing import List, Optional, Sequence, Tuple, Union

import pytest

from atropos.commands.trim.console import TrimCommandConsole

from .utils import assert_files_equal


class DatadirWrapper:
    def __init__(self, datadir, subdir: str):
        self.datadir = datadir
        self.subdir = subdir

    def __call__(self, name: str, **kwargs):
        name = name.format(**kwargs)
        return self.datadir[f"{self.subdir}/{name}"]

    def url(self, name: str, **kwargs):
        return f"file:{self(name, **kwargs)}"


@pytest.fixture(scope="function")
def input_data(datadir):
    return DatadirWrapper(datadir, "input")


@pytest.fixture(scope="function")
def expected_data(datadir):
    return DatadirWrapper(datadir, "expected")


@pytest.fixture(scope="function")
def run_trimmer(tmp_path: Path, capsys, input_data, expected_data):
    def _run_atropos_trim(
        params: Union[str, list],
        expected1=None,
        inpath1=None,
        inpath2=None,
        qualfile=None,
        interleaved_input=False,
        expected2=None,
        interleaved_output=False,
        aligner: str = "adapter",
        sra_accn=None,
        stdout=False,
        expected_multi: Optional[Sequence[str]] = None,
        expected_other: Optional[Sequence[str]] = None,
        assert_output_equal: bool = True,
        error_on_rc: bool = True,
    ) -> Tuple[Path, List[Path], dict]:
        if isinstance(params, str):
            param_list = params.split()
        else:
            param_list = list(params)

        to_compare = []

        output = tmp_path / expected1
        output2 = None
        if expected2:
            output2 = tmp_path / expected2

        if expected_multi:
            if stdout:
                raise ValueError("'stdout' and 'expected_multi' are mutually exclusive")

            for name in expected_multi:
                multi_path = Path(str(output).format(name=name, aligner=aligner))
                to_compare.append(
                    (expected_data(multi_path.name), multi_path)
                )
        elif assert_output_equal:
            to_compare.append((expected_data(expected1, aligner=aligner), output))

            if output2:
                to_compare.append(
                    (expected_data(expected2, aligner=aligner), output2)
                )

        if expected_other:
            for other in expected_other:
                i = param_list.index(other)
                other_path = tmp_path / other
                param_list[i] = other_path
                to_compare.append((expected_data(other, aligner=aligner), other_path))

        infiles = [
            input_data(infile, aligner=aligner) if infile else None
            for infile in (inpath1, inpath2)
        ]

        if sra_accn:
            param_list.extend(("--accession", sra_accn))
        elif interleaved_input:
            param_list.extend(("-l", infiles[0]))
        elif inpath2:
            param_list.extend(("-pe1", infiles[0]))
            param_list.extend(("-pe2", infiles[1]))
        else:
            param_list.extend(("-se", infiles[0]))
            if qualfile:
                param_list.extend(("-sq", infiles[0]))

        param_list.extend(("--aligner", aligner))

        if stdout:
            # Output is going to stdout, so we need to redirect it to the temp file
            result = TrimCommandConsole.execute(
                tuple(str(p) for p in param_list)
            )
            captured = capsys.readouterr()
            with open(output, "wt") as out:
                out.write(captured.out)
        else:
            if interleaved_output:
                param_list.extend(("-L", output))
            else:
                param_list.extend(("-o", output))  # TODO not parallelizable
                if output2:
                    param_list.extend(("-p", output2))

            result = TrimCommandConsole.execute(
                tuple(str(p) for p in param_list)
            )

        retcode, summary = result

        assert summary is not None
        assert isinstance(summary, dict)

        if "exception" in summary and summary["exception"] is not None:
            assert retcode != 0
            err = summary["exception"]
            traceback.print_exception(*err["details"])
            if error_on_rc:
                raise AssertionError(
                    f"Return code {retcode} != 0; error: {err['message']}"
                )
        elif error_on_rc:
            assert retcode == 0

        for expected_path, actual_path in to_compare:
            assert actual_path.exists()
            assert_files_equal(expected_path, actual_path)

        return tmp_path, infiles, summary

    return _run_atropos_trim
