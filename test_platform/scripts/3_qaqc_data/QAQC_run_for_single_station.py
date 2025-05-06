"""
QAQC_run_for_single_station.py

This script applies QA/QC protocols to cleaned weather station data for ingestion into the Historical Observations Platform.
It is designed to work on individual stations, regardless of network.

Overview:
---------
1. Remove duplicate station entries
2. Handle variables with inconsistent reporting intervals or changing frequencies (e.g., convert to hourly data)
3. Apply various QA/QC tests:
   - Consistency checks
   - Missing data and gap detection
   - Climatological outlier detection
   - Cross-variable sanity checks
4. (Optional) Perform case study analysis for accuracy
   [Note: This can be separated into a standalone script if needed]

Inputs:
-------
- Station ID (string) corresponding to a cleaned file in either zarr or netCDF format.

Outputs:
--------
- QA/QC-processed data for the specified station, including priority variables for all timestamps.
  The data is exported as a zarr file.

Example usage:
--------------
python QAQC_run_for_single_station.py --station="LOXWFO_CBGC1"

"""

import argparse
from QAQC_pipeline import run_qaqc_one_station


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
        prog="QAQC_run_for_single_station",  # Program name
        description="""This script applies quality assurance and quality control (QA/QC) protocols to cleaned weather 
                       station data for ingestion into the Historical Observations Platform (HOP). It is designed to be 
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
        "-r",
        "--rad_scheme",
        default="remove_zeros",
        help="Radiation handling scheme for frequent values check. Refer to 'qaqc_frequent_values' for options. Default is 'remove_zeros'.",
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
    run_qaqc_one_station(
        station=args.station,
        verbose=args.verbose,
        rad_scheme=args.rad_scheme,
    )


if __name__ == "__main__":
    main()
