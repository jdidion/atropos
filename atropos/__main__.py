import sys
from typing import Optional, Sequence

from atropos.console import run_atropos


def main(args: Optional[Sequence[str]] = None):
    """
    Main method, designed to be called by the system.

    Always calls `sys.exit()`. If you want to run Atropos programatically, call
    `run()` instead.

    Args:
        args: Command-line arguments. If None, `sys.argv[1:]` is used.
    """
    sys.exit(run_atropos(args))


if __name__ == "__main__":
    main()
