"""
General-purpose utility functions used throughout the HomologyViz application.

This module provides helper functions for common tasks such as file deletion,
directory cleanup, and locating package resources. These are used internally by
multiple components (e.g., BLAST preparation, temporary file handling).

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

import os
import importlib.resources as resources
from pathlib import Path
import shutil
import math
import platform
import subprocess
import tempfile

from homologyviz.logger import get_logger

logger = get_logger(__name__)


def get_package_path(package: str = "homologyviz") -> Path:
    """
    Return the filesystem path to the root directory of the specified package.

    Useful in `src/`-layout projects for locating bundled resources (e.g., templates,
    static files) at runtime. This uses Python's `importlib.resources` to safely access
    installed package data in a cross-platform way.

    Parameters
    ----------
    package : str, default="homologyviz"
        The name of the package whose base path is being retrieved.

    Returns
    -------
    path : pathlib.Path
        Filesystem path to the package directory.
    """
    return resources.files(f"{package}")


def get_os() -> str:
    """
    Detect the current operating system.

    Returns
    -------
    str
        A string identifying the OS: "Windows", "macOS", "Linux", or "Unknown".
    """
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
    """
    Check whether the BLASTn program is installed and accessible from the system PATH.

    This function attempts to run `blastn -version` and returns True if the command
    executes successfully, indicating that BLASTn is available.

    Returns
    -------
    bool
        True if BLASTn is installed and accessible, False otherwise.
    """
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


def delete_files(documents: list[Path]) -> None:
    """
    Delete a list of files from the filesystem.

    Iterates through the provided list of file paths and attempts to delete each one.
    If a file does not exist, a message is printed and the function continues without
    raising an error.

    Parameters
    ----------
    documents : list
        List of file paths (as strings or Path-like objects) to be deleted.

    Returns
    -------
    None
    """
    for document in documents:
        if os.path.exists(document):
            os.remove(document)
            logger.info(f"Deleted file: {document}")
        else:
            logger.warning(f"File not found: {document}")


def clean_directory(directory_path: Path) -> None:
    """
    Recursively delete all files and subdirectories from the specified directory.

    This function removes all contents of the directory, including nested files and non-
    empty subdirectories. The target directory itself is not deleted, only its contents.

    Parameters
    ----------
    directory_path : pathlib.Path
        Path to the directory to be cleaned.

    Returns
    -------
    None
    """
    if not directory_path.exists():
        return

    for item in directory_path.iterdir():
        if item.is_file() or item.is_symlink():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item)


def round_up_to_nearest_significant_digit(number: float) -> int:
    """
    Round a number up to the nearest multiple of its most significant digit.

    This function rounds a given number up to the nearest "clean" value based on its most
    significant digit. It is useful for generating scale bar or axis limits in plots that
    are easy to interpret in plots.

    This function is used in HomologyViz to define the scale bar values.

    Examples
    --------
    >>> round_up_to_nearest_significant_digit(142)
    200
    >>> round_up_to_nearest_significant_digit(89)
    90
    >>> round_up_to_nearest_significant_digit(5)
    5

    Parameters
    ----------
    number : float
        The number to round up.

    Returns
    -------
    int
        The input number rounded up to the nearest multiple of its most significant digit.
    """
    if number <= 0:
        raise ValueError("Input must be a positive number grater than zero.")
    # Determine the nearest power of ten (e.g., 1000, 100, 10, etc.)
    power_of_ten = 10 ** math.floor(math.log10(number))
    # Round up to the next multiple of that power
    return math.ceil(number / power_of_ten) * power_of_ten


if __name__ == "__main__":
    # test_path = Path("/Users/msp/Documents/testing_function")
    # clean_directory(test_path)
    print(get_package_path())
    with tempfile.TemporaryDirectory() as tmp_dir:
        file1 = Path(tmp_dir) / "existing.txt"
        file1.write_text("dummy")
        missing = Path("missing.txt")
        delete_files([file1, missing])
