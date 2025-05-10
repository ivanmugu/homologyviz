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
        else:
            print(f"File {document} does not exist")


def clean_directory(directory_path: Path) -> None:
    """
    Delete all files and empty subdirectories from the specified directory.

    If the directory is empty, the function does nothing. If files or empty folders
    are present, they are deleted. Subdirectories must be empty; otherwise,
    `rmdir()` will raise an error.

    Parameters
    ----------
    directory_path : pathlib.Path
        Path to the directory to be cleaned.

    Returns
    -------
    None

    # TODO: Add support for recursively deleting non-empty subdirectories.
    """
    if not any(directory_path.iterdir()):
        return
    else:
        for item in directory_path.glob("*"):
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                item.rmdir()


if __name__ == "__main__":
    print(get_package_path())
