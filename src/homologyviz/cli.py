"""Parse command line arguments provided by user.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import argparse
from argparse import Namespace
import pkg_resources
import subprocess
import sys
import platform

# TODO: test it on Windows


def parse_command_line_input() -> Namespace:
    """Parse command line arguments provided by user."""
    # Create parser.
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="msplotly",
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


def get_os() -> str:
    """Detect the current operating system."""
    os_name = platform.system()
    if os_name == "Windows":
        return "Windows"
    elif os_name == "Darwin":
        return "macOS"
    elif os_name == "Linux":
        return "Linux"
    else:
        return "Unknown"


def is_blastn_installed() -> bool:
    """Check if BLASTn is installed and accessible from the PATH."""
    current_os = get_os()
    shell_flag = current_os == "Windows"  # Set shell=True only for Windows

    # Attempt to run `blastn` and capture the output
    try:
        result = subprocess.run(
            ["blastn", "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=shell_flag,
        )
        if result.returncode == 0:
            print(f"BLASTn version: {result.stdout.strip()}")
            return True
        else:
            print(f"BLASTn found but returned an error: {result.stderr.strip()}")
    except FileNotFoundError:
        print("BLASTn is not installed or not found in the system PATH.")
    except Exception as e:
        print(f"An unexpected error occurred while checking for BLASTn: {e}")

    return False


if __name__ == "__main__":
    parse_command_line_input()
