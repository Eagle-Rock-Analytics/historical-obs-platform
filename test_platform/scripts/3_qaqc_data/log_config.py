"""
This is a script setting up Stage 3: QA/QC runtime logging for error tracing and timing. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

import logging
from mpi4py import MPI
import os


# Configure the logger
def setup_logger(log_file=f"{os.getcwd()}/default_qaqc_log.log", verbose=False):
    """Configures logger for more efficient tracing of QAQC processes and errors.

    Parameters
    ----------
    log_file : str
        Path to QAQC log file
    verbose : boolean, optional
        Prints script progress to local terminal. Default is False.

    Returns
    -------
    logger : [?]
        Information to be sent to logger / log file
    """

    # Retrieve the existing logger (or create a new one if it doesn't exist)
    logger = logging.getLogger("sharedLogger")
    logger.setLevel(logging.DEBUG)

    ## Check if the logger already has handlers to avoid duplication

    # Create a file handler that logs to a file
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)

    # Create a formatter and attach it to the handler
    formatter = logging.Formatter(
        "%(asctime)s - %(rank)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)

    # Add the file handler to the logger
    logger.addHandler(file_handler)

    # Optionally, add a stream handler to log to console
    if verbose:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # Add a custom filter to include the MPI rank in the log messages
    class RankFilter(logging.Filter):
        def filter(self, record):
            # Attach the MPI rank to the log record
            record.rank = MPI.COMM_WORLD.Get_rank()
            record.rank = "(Rank {}/{})".format(
                MPI.COMM_WORLD.Get_rank() + 1, MPI.COMM_WORLD.Get_size()
            )
            return True

    # Add the filter to the logger
    logger.addFilter(RankFilter())

    return logger


def remove_file_handler_by_filename(logger, filename):
    """Remove a specific FileHandler from the logger by matching the filename.

    Parameters
    ----------
    logger : [?]
        Input logger handler
    filename : str
        Filename of logger handler

    Returns
    -------
    None
        This function does not return a value

    """
    ## NEEDS DOCUMENTATION IMPROVEMENT

    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            # Check if the handler's baseFilename matches the specified filename
            if handler.baseFilename == filename:
                logger.removeHandler(handler)
                handler.close()  # Close the handler after removing it
                # print(f"Removed FileHandler for {filename}")
                break  # Exit after removing the first matching handler
    else:
        print(f"No FileHandler found for {filename}")

    return None


# Call this to initialize the logger
logger = setup_logger()

# Now, let's remove the handler for default log file to avoid duplication in printing to console (verbose=True)
remove_file_handler_by_filename(logger, f"{os.getcwd()}/default_qaqc_log.log")
