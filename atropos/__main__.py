import sys
from atropos import check_importability
from atropos.commands import execute_cli
from typing import Sequence


def main(args: Sequence[str] = sys.argv[1:]) -> None:
    """Main method.

    Args:
        args: Command-line arguments.
    """
    if not check_importability():
        print(
            f"""
ERROR: A required extension module could not be imported because it is
incompatible with your system. A quick fix is to recompile the extension
modules with the following command:

{sys.executable} setup.py build_ext -i

See the documentation for alternative ways of installing the program.

The original error message follows.
""")
        rc = 1
    else:
        rc = execute_cli(args)
    sys.exit(rc)
