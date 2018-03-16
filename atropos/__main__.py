import sys
from atropos import check_importability
from atropos.commands import execute_cli


def main(args=sys.argv[1:]):
    """Main method.

    Args:
        args: Command-line arguments.
    """
    check_importability()
    sys.exit(execute_cli(args))


if __name__ == '__main__':
    main()
