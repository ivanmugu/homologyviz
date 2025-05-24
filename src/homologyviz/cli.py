"""
Command-line interface (CLI) utilities for HomologyViz.

This module provides functions to parse command-line arguments and check system
dependencies (e.g., presence of BLASTn). It supports the GUI launcher and can be
extended for future command-line features.

Notes
-----
- This file is part of HomologyViz
- BSD 3-Clause License
- Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

import argparse
from argparse import Namespace
import pkg_resources
import sys

from homologyviz.miscellaneous import is_blastn_installed

# TODO: test it on Windows, and Ubuntu


def parse_command_line_input() -> Namespace:
    """
    Parse command-line arguments and validate the environment.

    This function sets up the command-line interface (CLI) for HomologyViz using
    `argparse`. It defines helper flags for displaying the help message and the program
    version, then parses user-provided arguments.

    Additionally, it checks if `blastn` is installed locally. If `blastn` is not found,
    the program exits early to prevent runtime errors.

    Returns
    -------
    argparse.Namespace
        A namespace object containing parsed command-line arguments.
        Typically empty unless extended in the future.

    Notes
    -----
    - `--help` or `-h`: Show usage information.
    - `--version` or `-v`: Show the installed version of the CLI.
    - Exits the program if `blastn` is not installed locally.
    """
    # Create parser.
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="homologyviz",
        formatter_class=argparse.RawTextHelpFormatter,
        description=("Make a graphical representation of BLASTn alignments."),
    )
    # Make argument groups.
    helper = parser.add_argument_group("Help")
    helper.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit."
    )
    prog_version = pkg_resources.get_distribution("msplotly").version
    helper.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {prog_version}",
        help="Show program's version number and exit",
    )

    # Parse command line arguments
    command_line_info = parser.parse_args()

    # check if user has BLASTn installed locally.
    if not is_blastn_installed():
        sys.exit()

    return command_line_info


if __name__ == "__main__":
    parse_command_line_input()
