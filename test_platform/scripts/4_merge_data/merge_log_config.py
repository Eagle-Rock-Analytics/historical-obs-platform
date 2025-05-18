"""
merge_log_config.py

Provides a function to configure a timestamped logger that writes to a file
and optionally to the console.

This logger is intended for QAQC routines, data processing workflows, or general
debugging where timestamped and structured logs are helpful.

Logging Behavior
----------------
The logger is set to level DEBUG, so all log messages are captured. Handlers
(file, console) can filter what gets recorded using their own level settings:

- File handler logs DEBUG and above.
- Console handler logs DEBUG and above if `verbose=True`.

Log Format
----------
Each log entry is formatted as:

    %(asctime)s - %(levelname)s - %(message)s

Example
-------
>>> from merge_log_config import init_logger
>>> logger = init_logger("logger.log", logs_dir="qaqc_logs", verbose=True)
>>> logger.info("Logger initialized")
"""

import logging
import os
from datetime import datetime


def _configure_logger(log_file: str, verbose: bool) -> logging.Logger:
    """
    Configures a logger to write to a file, with optional console output.

    Parameters
    ----------
    log_file : str
        Full path to the log file. Will overwrite existing file if present.
    verbose : bool
        If True, logs are also printed to the console.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    # Create or retrieve a logger instance named "sharedLogger"
    logger = logging.getLogger("sharedLogger")

    # If the logger is already configured with handlers, return it as-is
    if logger.hasHandlers():
        return logger

    # Set the logger to capture all messages at DEBUG level and above
    logger.setLevel(logging.DEBUG)

    # Create a file handler in overwrite mode
    file_handler = logging.FileHandler(log_file, mode="w")  # Overwrite if file exists
    file_handler.setLevel(logging.DEBUG)

    # Define log message format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)

    # Add the file handler to the logger
    logger.addHandler(file_handler)

    # Optionally, also add console logging if verbose is True
    if verbose:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    return logger


def init_logger(
    log_fname: str, logs_dir: str = "qaqc_logs", verbose: bool = True
) -> logging.Logger:
    """
    Initializes a timestamped logger with optional console output.

    Parameters
    ----------
    log_fname : str
        Log filename (e.g., station name, script name).
    logs_dir : str, optional
        Directory to store log files. Defaults to "qaqc_logs".
    verbose : bool, optional
        If True, also logs to the console. Defaults to True.

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """

    # Create the logs directory if it does not exist
    if not os.path.isdir(logs_dir):
        os.makedirs(logs_dir, exist_ok=True)
        print(f"Created directory: {logs_dir}")

    # Construct the full path to the log file
    log_path = os.path.join(logs_dir, log_fname)

    # Configure and return the logger
    return _configure_logger(log_path, verbose=verbose)
