"""Miscellaneous functions.

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Ivan Munoz Gutierrez
"""

import os
import importlib.resources as resources
from pathlib import Path
import base64
import binascii


def testing_base64(content: str):
    try:
        base64.b64decode(content)
        print("good base64")
    except binascii.Error as e:
        print(e)
        print("Invalid base64 data:", e)


def get_package_path(package: str = "msplotly") -> Path:
    """Get path to package directory in a src-layout"""
    return resources.files(f"{package}")


def delete_files(documents: list) -> None:
    """Delete the files from `documents` list."""
    for document in documents:
        if os.path.exists(document):
            os.remove(document)
        else:
            print(f"File {document} does not exist")


def clean_directory(directory_path: Path) -> None:
    """If directory is not empty, delete all files"""
    if not any(directory_path.iterdir()):
        return
    else:
        for item in directory_path.glob("*"):
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                item.rmdir()


if __name__ == "__main__":
    # print(get_package_path())
    content_good = "data:text/plain;base64,SGVsbG8sIHdvcmxkIQ=="
    content_bad = "this_is_not_base64"
    testing_base64(content_good)
