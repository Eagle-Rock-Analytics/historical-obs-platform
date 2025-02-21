"""
This script is a template structure for data merging across all available data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Reads in all qa/qc processed data files
(2) Merges data into a single station file
(3) Writes necessary information regarding data product, focusing on flexible usage
Inputs: All QA/QC-processed data for each network
Outputs: Final data product as .zarr (or.nc file)

Example how to run entire network in the command line: 
python3 ALLNETWORKS_merge.py -n="ASOSAWOS"
Example how to run a sample of a network in the command line with verbose print statements:
python3 ALLNETWORKS_merge.py -n="ASOSAWOS" -s 4 -v True
"""

# Import libraries
import argparse

# Import merge pipeline
try:
    from MERGE_pipeline import *
except Exception as e:
    print("Error importing MERGE_pipeline.py: {}".format(e))

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
    sample = args.sample

    # Set zarr argument -- this may no longer be true as a part of QAQC
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
    run_merge_pipeline(
        network,
        qaqcdir,
        mergedir,
        verbose=verbose,
        sample=sample,
        zarr=zarr,
    )
