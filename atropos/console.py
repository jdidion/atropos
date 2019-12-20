"""
Atropos version {}

usage: atropos [--config <config file>] <command> [options]

commands
--------
{}

optional arguments:
  -h, --help                show this help message and exit
  --config <config file>    provide options in a config file

Use "atropos <command> --help" to see all options for a specific command.
See http://atropos.readthedocs.org/ for full documentation.

Atropos is a fork of Cutadapt 1.10 (
https://github.com/marcelm/cutadapt/tree/2f3cc0717aa9ff1e0326ea6bcb36b712950d4999)
by John Didion, et al., "Atropos: sensitive, specific, and speedy trimming of
NGS reads, submitted.

Cutadapt (https://github.com/marcelm/cutadapt) was developed by Marcel Martin,
"Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads,"
EMBnet Journal, 2011, 17(1):10-12.
"""
import sys
from typing import Dict, Optional, Sequence, TypeVar, cast

from loguru import logger
import pkg_resources
from xphyle import open_

from atropos import __version__
from atropos.commands.console import CommandConsole
from atropos.utils import ReturnCode
from atropos.utils.paths import as_readable_path


HELP_TEXT = __doc__
T = TypeVar("T")


def run_atropos(args: Optional[Sequence[str]] = None) -> ReturnCode:
    """
    Run atropos.

    Args:
        args: Command-line arguments. If None, `sys.argv[1:]` is used.

    Returns:
        The return code.
    """
    if check_importability():
        if args is None:
            args = sys.argv[1:]
        return execute_cli(args)
    else:
        logger.error(
            f"""
ERROR: A required extension module could not be imported because it is
incompatible with your system. A quick fix is to recompile the extension
modules with the following command:

{sys.executable} setup.py build_ext -i

See the documentation for alternative ways of installing the program.

The original error message follows.
"""
        )
        return ReturnCode.ERROR


def check_importability() -> bool:  # pragma: no cover
    """
    Check that Cython modules have been compiled.

    Returns:
        True if a cython-compiled module can be imported, otherwise False.
    """
    try:
        import atropos.aligner  # pylint: disable=unused-variable
        return True
    except ImportError as err:
        if "undefined symbol" in str(err):
            return False
        else:
            raise


def execute_cli(args: Sequence[str] = ()) -> ReturnCode:
    """
    Executes the Atropos command-line interface.

    The first argument is expected to be the command name. If not (i.e. args is
    empty or the first argument starts with a '-'), the 'trim' command is
    assumed. If the first argument is '-h' or '--help', the command-level help
    is printed.

    Args:
        args: Command-line arguments.

    Returns:
        The return code.
    """
    commands = discover_commands()

    if len(args) == 0 or args[0] in ("-h", "--help"):
        print(
            HELP_TEXT.format(
                __version__,
                "\n".join(cmd.get_help() for cmd in commands.values())
            ),
            file=sys.stderr
        )
        return ReturnCode.ERROR

    if args[0] == "--config":
        config_path = as_readable_path(args[1])
        args = args[2:]

        with open_(config_path, "rt") as config_file:
            config_args = list(
                token for line in config_file for token in line.rstrip().split()
            )
    else:
        config_args = None

    def parse_command(_args):
        if not _args or _args[0][0] == "-":
            return "trim", _args
        else:
            return _args[0], _args[1:]

    if len(args) == 0:
        command_name, args = parse_command(config_args)
    else:
        command_name, args = parse_command(args)
        if config_args:
            args = config_args + args

    try:
        command = commands[command_name]
        retcode, summary = command.execute(args)
        if "exception" in summary:
            logger.opt(exception=summary["exception"]["details"]).error(
                f"Error executing command {command_name}"
            )
        return retcode
    except:
        logger.exception(f"Error executing command: {command_name}")
        return ReturnCode.ERROR


def discover_commands() -> Dict[str, CommandConsole]:
    """
    Discover Commands defined as entrypoints.
    """
    return dict(
        (cmd_class.name, cmd_class)
        for cmd_class in (
            cast(CommandConsole, entry_point.load(True))
            for entry_point in pkg_resources.iter_entry_points(group="atropos.commands")
        )
    )
