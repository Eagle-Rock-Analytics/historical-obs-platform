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


def _pairwise_rate(
    flag_df: pd.DataFrame, running_rate_df, station_name: str
) -> pd.DataFrame:
    """
    Generates flag rates dataframe for input flag counts dataframe and then adds it to the running flag rate dataframe.
    Helper function for network_rate_tables() and station_rate_table().

    Parameters
    ----------
    flag_df: pd.DataFrame
        flag rates dataframe for next station
    running_rate_df: pd.DataFrame
        dataframe of previously added station flag rates

    Returns
    -------
    rates_df_merged: pd.DataFrame

    """
    # Make the eraqc_flag_values column the index
    flag_df = flag_df.set_index("eraqc_flag_values")

    # Count up the flagged observations - so counts in all but the "no_flag" and "total_obs_count" rows
    subset = flag_df[~flag_df.index.isin(["no_flag", "total_obs_count"])]
    totals = subset.sum(numeric_only=True)
    flag_df.loc["total_flag"] = pd.Series(totals)

    # And then use those total to calculate the per-variable flag rates
    frac = flag_df.loc["total_flag"] / flag_df.loc["total_obs_count"]
    flag_df.loc["flag_rate"] = pd.Series(frac)

    # Keep only the rate
    rates_df = flag_df.loc[["flag_rate"]]
    rates_df = rates_df.rename(index={"flag_rate": station_name})

    rates_df = rates_df.reset_index()

    # Finally, append column of total observation count
    flag_df = flag_df.reset_index()
    total_obs = flag_df[flag_df["eraqc_flag_values"] == "total_obs_count"].iloc[0, 1]
    rates_df["total_obs_count"] = total_obs

    if len(running_rate_df) == 0:
        return rates_df

    else:
        rates_df_merged = pd.merge(rates_df, running_rate_df, how="outer")
        return rates_df_merged


def network_rate_tables(timestep: str) -> None:
    """
    Generates a table of flag rates per network and uploads it to AWS.

    Parameters
    ----------
    timestep: str
        if set to 'hourly', merge all hourly QAQC flag count tables
        if set to 'native', merge all native timestep QAQC flag count tables

    Returns
    -------
    None

    """
    # Record start time
    start_time = time.time()

    ## Setup

    # Only run for a valid "timestep" input
    if timestep not in ("hourly", "native"):
        print("invalid timestep: ", timestep)
        return None

    # The function iteratively adds in flag counts to this dataframe
    flag_rate_df = []

    # Point to folder containing station flag count CSVs
    flags_prefix = f"{MERGE_DIR}/per_network_flag_counts_{timestep}_timestep"

    ## Merge flag counts

    # Loop through all CSVs are the given level
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=flags_prefix):
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item.key)
        flags = pd.read_csv(obj["Body"])
        station_name = item.key.split(flags_prefix + "/")[1].split("_flag")[0]
        # the CSV is empty
        if flags.empty:
            continue
        # the CSV is not empty
        else:
            # Remove "QAQC_function" and "Flag_meaning" columns - we don't need these
            flags = flags.drop(["QAQC_function", "Flag_meaning"], axis=1)

            # Send current dataframe and dataframe of previously generated rates to helper function
            flag_rate_df = _pairwise_rate(flags, flag_rate_df, station_name)

    # Change "eraqc_flag_values" to "networks"
    flag_rate_df = flag_rate_df.rename(columns={"eraqc_flag_values": "networks"})

    ## Send final flag rates file to AWS as CSV
    csv_s3_filepath = (
        f"s3://wecc-historical-wx/4_merge_wx/network_{timestep}_flag_rates.csv"
    )

    print(f"Sending {timestep} timestep network flag rates CSV to: {csv_s3_filepath}")
    flag_rate_df.to_csv(csv_s3_filepath, index=False)

    ## Output time elapsed
    end_time = time.time()
    time_elapsed = (end_time - start_time) / 60
    print(f"{time_elapsed} minutes")

    return None


def station_rate_tables(timestep: str) -> None:
    """
    Generates a table of flag rates per station.

    Parameters
    ----------
    timestep: str
        if set to 'hourly', merge all hourly QAQC flag count tables
        if set to 'native', merge all native timestep QAQC flag count tables
    Returns
    -------
    None

    """
    # Record start time
    start_time = time.time()

    ## Setup

    # Only run for a valid "timestep" input
    if timestep not in ("hourly", "native"):
        print("invalid timestep: ", timestep)
        return None

    # List of networks to iterate over
    station_list = pd.read_csv(stations_csv_path)
    network_list = station_list["network"].unique()

    network_list = ["ASOSAWOS", "SNOTEL"]

    # The function iteratively adds in flag counts to this dataframe
    flag_rate_df = []

    for network in network_list:
        # Point to folder containing station flag count CSVs
        flags_prefix = f"{MERGE_DIR}/{network}/eraqc_counts_{timestep}_timestep"

        ## Merge flag counts

        # Loop through all CSVs at the given level
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=flags_prefix):
            obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item.key)
            flags = pd.read_csv(obj["Body"])
            station_name = item.key.split(flags_prefix + "/")[1].split("_flag")[0]
            # the CSV is empty
            if flags.empty:
                continue
            # the CSV is not empty
            else:
                # Send current dataframe and dataframe of previously generated rates to helper function
                flag_rate_df = _pairwise_rate(flags, flag_rate_df, station_name)

    # Change "eraqc_flag_values" to "stations"
    flag_rate_df = flag_rate_df.rename(columns={"eraqc_flag_values": "networks"})

    ## Send final flag rates file to AWS as CSV
    csv_s3_filepath = (
        f"s3://wecc-historical-wx/4_merge_wx/station_{timestep}_flag_rates.csv"
    )

    print(f"Sending {timestep} timestep station flag rates CSV to: {csv_s3_filepath}")
    flag_rate_df.to_csv(csv_s3_filepath, index=False)

    ## Output time elapsed
    end_time = time.time()
    time_elapsed = (end_time - start_time) / 60
    print(f"{time_elapsed} minutes")

    return None


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

    # step 1: generate per network flag rate table
    print(
        "Starting flag rate table generation per network -- anticipated time to complete: 3 sec"
    )

    print("Generating native timestep network rates table...")
    network_rate_tables("native")

    print("Generating hourly timestep network rates table...")
    network_rate_tables("hourly")

    # step 2: generate per station flag rate table
    print(
        "Starting flag rate table generation per station -- anticipated time to complete: 1 hr"
    )

    print("Generating native timestep station rates table...")
    station_rate_tables("native")

    print("Generating hourly timestep station rates table...")
    station_rate_tables("hourly")

    return None


if __name__ == "__main__":
    main()
