"""
merge_log_config.py

Logging utility functions for QAQC workflows and data processing scripts.

Provides functions to:
- Configure a timestamped file-based logger (`setup_logger`)
- Upload a log file to S3 (`upload_log_to_s3`)

These functions are intended to be imported into other scripts. Internal-use
functions (prefixed with underscores) handle logger configuration logic.

Logging Behavior
----------------
The logger captures all log levels (`DEBUG` and above). Handlers filter what
actually gets recorded:

- File handler logs DEBUG and above.
- Console handler logs DEBUG and above if `verbose=True`.

Log Format
----------
Each log entry is formatted as:

    %(asctime)s - %(levelname)s - %(message)s

Usage Example
-------------
>>> from merge_log_config import setup_logger
>>> logger, log_path = setup_logger("station123", logs_dir="logs", verbose=True)
>>> logger.info("Logger initialized.")
"""

import logging
import os
import inspect
from datetime import datetime

import boto3


def setup_logger(
    station: str, logs_dir: str = "merge_logs", verbose: bool = True
) -> tuple[logging.Logger, str]:
    """
    Initialize a timestamped logger for the current station.

    Parameters
    ----------
    station : str
        Station identifier for log file naming.
    logs_dir : str
        Directory to save logs. Defaults to "merge_logs".
    verbose : bool
        If True, also logs to console.

    Returns
    -------
    logger : logging.Logger
        Configured logger instance.
    log_filepath : str
        Full path to the created log file.
    """
    timestamp = datetime.now().strftime("%Y%m%d%H%M")
    log_filename = f"merge_{station}.{timestamp}.log"
    log_filepath = f"{logs_dir}/{log_filename}"

    logger = _init_logger(log_filename, logs_dir=logs_dir, verbose=verbose)
    logger.info(f"Starting merge script for station: {station}")

    return logger, log_filepath


def upload_log_to_s3(
    bucket_name: str,
    merge_dir: str,
    network: str,
    logfile_fpath: str,
    logger: logging.Logger,
) -> None:
    """
    Upload a local log file to an S3 bucket and log the operation.

    Parameters
    ----------
    bucket_name : str
        Name of the S3 bucket.
    merge_dir : str
        Prefix (directory path) inside the bucket to store the log file.
    network: str
        Name of network. Used to construct upload path for logfile
    logfile_fpath : str
        Path to local logfile
    logger : logging.Logger
        Logger instance for logging actions and errors.

    Raises
    ------
    Exception
        If upload to S3 fails, the original exception is raised.
    """

    # Construct full S3 URI for logging
    s3_key = f"{merge_dir}/{network}/{logfile_fpath}"
    logfile_s3_uri = f"s3://{bucket_name}/{s3_key}"

    # Upload to s3
    s3 = boto3.resource("s3")
    try:
        s3.Bucket(bucket_name).upload_file(logfile_fpath, s3_key)
        logger.info(f"Saved log file to {logfile_s3_uri}")
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to upload log file to S3: {e}"
        )
        raise

    _close_logger(logger)


def _close_logger(logger: logging.Logger) -> None:
    """
    Closes all handlers associated with the given logger.

    Parameters
    ----------
    logger : logging.Logger
        The logger instance to close.
    """
    for handler in logger.handlers[:]:  # Copy the list to avoid modification issues
        handler.close()
        logger.removeHandler(handler)


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


def _init_logger(
    log_fname: str, logs_dir: str = "merge_logs", verbose: bool = True
) -> logging.Logger:
    """
    Initializes a timestamped logger with optional console output.

    Parameters
    ----------
    log_fname : str
        Log filename (e.g., station name, script name).
    logs_dir : str, optional
        Directory to store log files. Defaults to "merge_logs".
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
