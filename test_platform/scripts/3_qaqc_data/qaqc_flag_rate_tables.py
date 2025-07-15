"""
qaqc_success_report_tables.py

Creates QAQC flag counts csv files per network from the corresponding eraqc_counts_timestep files that were
generated as a part of the final processing step for stations within the Historical Data Pipeline.
These tables are used to then generate statistics for the QAQC success report.

This is carried out in two steps:

1. Generate the per-network QAQC flag count tables, at native and hourly timesteps
2. Generates one flag count table that sums all per-network tables, at native and hourly timesteps

Functions
---------
- _pairwise_sum(): helper function that merges two input flag tables, used by network_sum_flag_counts() and total_sum_flag_counts().
- _network_format_table:
- _total_format_table:
- network_sum_flag_counts(): sums all station flag count tables for a given network, creating one flag count table for that network
- generate_station_tables(): runs network_sum_flag_counts() for every network
- total_sum_flag_counts(): sums all network flag count tables, creating one final flag count table

Intended Use
------------
Run this script to generate all flag sum tables at the network and total levels.

Run "python qaqc_generate_flag_sum.py"

"""

import time
import boto3
import numpy as np
import pandas as pd
import xarray as xr

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Set relative paths to other folders and objects in repository.
BUCKET_NAME = "wecc-historical-wx"
QAQC_DIR = "3_qaqc_wx"
MERGE_DIR = "4_merge_wx"
stations_csv_path = f"s3://{BUCKET_NAME}/{QAQC_DIR}/all_network_stationlist_qaqc.csv"





def main():
    """
    Run this script.

    Parameters
    ----------
    timestep: str
        if set to 'hourly', merge all hourly QAQC flag count tables
        if set to 'native', merge all native timestep QAQC flag count tables

    Returns
    -------
    None

    """

    # step 1: generate flag sum tables for each network
    print(
        "Starting flag sum table generation per network -- anticipated time to complete: 1 hour"
    )

    print("Generating native timestep tables...")
    generate_station_tables("native")
    print("Generating hourly timestep tables...")
    generate_station_tables("hourly")

    # step 2: generate total flag sum table
    print("Starting total sum flag counts...")

    print("Generating native timestep tables...")
    total_sum_flag_counts("native")

    print("Generating hourly timestep tables...")
    total_sum_flag_counts("hourly")

    return None


if __name__ == "__main__":
    main()
