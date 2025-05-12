"""
This module provides a reusable function for setting up and retrieving
a logger instance with consistent formatting and behavior across the application.

Usage
-----
Import and call `get_logger(__name__)` in any module to obtain a logger
that logs to the console with timestamps, log levels, and module names.

Example
-------
from logger import get_logger

logger = get_logger(__name__)
logger.info("This is an informational message.")

License
-------
This file is part of HomologyViz
BSD 3-Clause License
Copyright (c) 2024, Iván Muñoz Gutiérrez
"""

import logging


def get_logger(name: str) -> logging.Logger:
    """
    Retrieve a configured logger instance for the given module name.

    This function ensures that each logger is only configured once,
    avoiding duplicate handlers. The logger is set to INFO level by default
    and outputs to the console with a standard format that includes the
    timestamp, module name, log level, and message.

    Parameters
    ----------
    name : str
        The name of the logger, typically passed as `__name__` from the calling module.

    Returns
    -------
    logging.Logger
        A logger instance configured for console output.
    """
    logger = logging.Logger(name)

    if not logger.handlers:
        logger.setLevel(logging.INFO)

        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)

        logger.addHandler(handler)
        logger.propagate = False

    return logger
