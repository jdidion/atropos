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
from importlib import import_module
import logging
import os
from pkgutil import walk_packages
import re
import textwrap
from atropos import __version__

class Command(object):
    """Contains information about a command package.
    
    A command package consists of at least 3 modules:
    1. The top level module (__init__.py), which contains the logic to run the
    command. It must have a top-level CommandRunner class that extends
    atropos.commands.base.BaseCommandRunner.
    2. The cli module (cli.py), which configures the ArgumentParser and
    validates command line options. It must have a top-level CommandParser
    class that extends atropos.commands.base.BaseCommandParser.
    3. The reports module (reports.py), which generates reports from the
    command summary. It must have a top-level ReportGenerator class that
    extends atropos.commands.reports.BaseReportGenerator.
    
    Args:
        name: The command name. Must match the name of a submodule of
            atropos.commands, or all of the other arguments must be
            specified.
        module, cli_module, report_module: The absolute names of the three
            modules described above. If None, they are determined automatically
            from the command name.
    """
    def __init__(self, name, module=None, cli_module=None, report_module=None):
        self.name = name
        self.package = module or 'atropos.commands.{}'.format(name)
        self.cli_module = cli_module or '{}.cli'.format(self.package)
        self.report_module = report_module or '{}.reports'.format(self.package)
    
    def execute(self, args=()):
        """Parse command line arguments, execute the command, and generate
        summary reports.
        
        Returns:
            The return code.
        """
        options = self.parse_args(args)
        retcode, summary = self.run_command(options)
        if retcode == 0 and options.report_file:
            logging.getLogger().debug(
                "Writing report to %s", options.report_file)
            self.generate_reports(summary, options)
        else:
            logging.getLogger().debug("Not generating report file")
        return retcode, summary
    
    def get_command_parser_class(self):
        """Returns the CommandParser class within the cli module.
        """
        mod = import_module(self.cli_module)
        return mod.CommandParser
    
    @property
    def usage(self):
        """Returns the command's usage string.
        """
        return self.get_command_parser_class().usage
    
    @property
    def description(self):
        """Returns the command's description string.
        """
        return self.get_command_parser_class().description
    
    def get_help(self, fmt="* {name}: {description}", wrap=80, indent=2):
        """Returns a string to include in the command help.
        """
        helpstr = fmt.format(
            name=self.name, description=self.description.strip())
        if wrap:
            helpstr = "\n".join(textwrap.wrap(
                re.sub(r'\s+', ' ', helpstr), wrap,
                subsequent_indent=' ' * indent))
        return helpstr
    
    def parse_args(self, args):
        """Parse the command line options.
        
        Returns:
            A Namespace-like object.
        """
        parser_class = self.get_command_parser_class()
        parser = parser_class()
        return parser.parse(args)
    
    def get_command_runner_class(self):
        """Returns the CommandRunner class within the top-level module.
        """
        mod = import_module(self.package)
        return mod.CommandRunner
    
    def run_command(self, options):
        """Run the command.
        
        Args:
            options: A Namespace-like object.
        
        Returns:
            Tuple (retcode, summary).
        """
        runner_class = self.get_command_runner_class()
        runner = runner_class(options)
        return runner.run()
    
    def get_report_generator_class(self):
        """Returns the ReportGenerator class within the reports module.
        """
        mod = import_module(self.report_module)
        return mod.ReportGenerator
    
    def generate_reports(self, summary, options):
        """Generate reports.
        
        Args:
            summary: The summary dict.
            options: Command-line options.
            report_file: The report file name/prefix.
            report_formats: A list of formats.
        """
        generator_class = self.get_report_generator_class()
        generator = generator_class(options)
        generator.generate_reports(summary)

COMMANDS = dict(
    (name, Command(name))
    for _, name, ispkg in walk_packages([os.path.dirname(__file__)])
    if ispkg)

def get_command(name):
    """Returns the Command object for the specified command name.
    """
    if name not in COMMANDS:
        raise ValueError("Invalid command: {}".format(name))
    return COMMANDS[name]

def iter_commands():
    """Iterate over commands, ordered by command name.
    """
    for name in sorted(COMMANDS.keys()):
        yield COMMANDS[name]

def execute_cli(args=()):
    """Executes the Atropos command-line interface.
    
    The first argument is expected to be the command name. If not (i.e. args is
    empty or the first argument starts with a '-'), the 'trim' command is
    assumed. If the first argument is '-h' or '--help', the command-level help
    is printed.
        
    Args:
        args: Command-line arguments.
    
    Returns:
        The return code.
    """
    if len(args) == 0 or args[0] in ('-h', '--help'):
        print_subcommands()
        return 2
    
    config_args = None
    
    if args[0] == '--config':
        with open(args[1], 'rt') as config_file:
            config_args = list(
                token
                for line in config_file
                for token in line.rstrip().split())
        args = args[2:]
    
    def parse_command(args):
        if not args or args[0][0] == '-':
            return ('trim', args)
        else:
            return (args[0], args[1:])
    
    if len(args) == 0:
        command_name, args = parse_command(config_args)
    else:
        command_name, args = parse_command(args)
        if config_args:
            args = config_args + args
    
    try:
        command = get_command(command_name)
        retcode, summary = command.execute(args)
        if 'exception' in summary:
            logging.getLogger().error(
                "Error executing command %s", command_name, 
                exc_info=summary['exception']['details'])
        return retcode
    except Exception as err:
        logging.getLogger().error(
            "Error executing command: %s", command_name, exc_info=err)
        return 2

def print_subcommands():
    """Prints usage message listing the available subcommands.
    """
    print(__doc__.format(__version__, "\n".join(
        command.get_help() for command in iter_commands())))
