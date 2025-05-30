"""
MERGE_run_for_single_station.py

Inputs:
-------
- Station ID (string) corresponding to a cleaned file in either zarr or netCDF format.

Outputs:
--------
- QA/QC-processed data for the specified station, including priority variables for all timestamps.
  The data is exported as a zarr file.

Example usage:
--------------
python MERGE_run_for_single_station.py --station="LOXWFO_CBGC1"

"""

import argparse
from MERGE_pipeline import run_merge_one_station


def main():
    """
    This function is designed to create an argument parser, define the arguments for the script, parse the argument, and then run the QAQC pipeline for the specified station.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """
    # Create argument parser
    parser = argparse.ArgumentParser(
        prog="MERGE_run_for_single_station",  # Program name
        description="""This script runs the full merge pipeline for a single weather station dataset.
                        It processes cleaned QA/QC output and applies a sequence of standardization,
                        homogenization, and export operations to generate a finalized station file. It is designed to be 
                       network-independent, capable of processing data for individual stations across various networks.""",
    )

    # Define arguments for the script
    parser.add_argument(
        "-s",
        "--station",
        required=True,
        help="The Station ID to process (required). This corresponds to a cleaned file in either zarr or netCDF format.",
        type=str,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        help="Enable verbose output for debugging and logging (default: False).",
        type=bool,
    )

    # Parse arguments
    args = parser.parse_args()

    # Run the QAQC pipeline for the specified station
    run_merge_one_station(
        station=args.station,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
