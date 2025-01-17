"""
This script is a template structure for data merging across all available data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Reads in all qa/qc processed data files
(2) Merges data into a single .nc file
(3) Writes necessary information regarding data product, focusing on flexible usage
Inputs: All QA/QC-processed data for each network
Outputs: Final data product as .nc file (or .zarr?)

Example how to run in the command line: 
python ALLNETWORKS_merge.py -n="ASOSAWOS"

"""

# Environment set-up
# Import libraries
import os
import tempfile
import argparse

# =================================================================================================
# Main Function
if __name__ == "__main__":
    # Create parser
    parser = argparse.ArgumentParser(
        prog="ALLNETWORKS_merge",
        description="""This script performs the final merge protocols for cleaned and quality controlled station data for ingestion into the Historical 
                       Observations Platform, and is independent of network.""",
        epilog="""Possible stations:
                  [ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, 
                   HOLFUY, HPWREN, LOXWFO, MAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, OtherISD, 
                   RAWS, SGXWFO, SHASAVAL, VALLEYWATER, VCAPCD, MARITIME, NDBC, SCAN, SNOTEL]
               """,
    )

    # Define arguments for the program
    parser.add_argument(
        "-n",
        "--network",
        default="TRAINING",
        help="Network name (default to 'TRAINING').",
        type=str,
    )
    parser.add_argument(
        "-l",
        "--local",
        default=False,
        help="Save files and plots locally (default to False).",
        type=bool,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        help="Print statements throughout script (default to False).",
        type=bool,
    )
    parser.add_argument(
        "-s",
        "--sample",
        default="all",
        help="How many stations to run (default to 'all'l).",
        type=str,
    )

    # Parse arguments
    args = parser.parse_args()
    network = args.network
    if network.lower() == "training":
        network = "TRAINING"
    verbose = args.verbose
    local = args.local
    sample = args.sample

    # Set zarr argument
    zarrified_networks = ["VALLEYWATER"]  # Any networks with zarrified data
    if network in zarrified_networks:
        zarr = True
    else:
        zarr = False

    # Set paths to data in AE bucket
    if network == "NETWORK":
        rawdir, cleandir, qaqcdir, mergedir = None, None, None, None
    else:
        rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
    whole_station_qaqc(
        network,
        cleandir,
        qaqcdir,
        rad_scheme,
        verbose=verbose,
        local=local,
        sample=sample,
        zarr=zarr,
    )