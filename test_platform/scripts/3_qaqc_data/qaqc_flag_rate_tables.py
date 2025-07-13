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
    flag_df_1: pd.DataFrame, flag_df_2: pd.DataFrame, station_name: str
) -> pd.DataFrame:
    """
    Sums two input flag count dataframes. This is a helper function for sum_flag_counts().

    Parameters
    ----------
    flag_df_1: pd.DataFrame
        dataframe of previously summed station flag counts
    flag_df_2: pd.DataFrame
        flag counts dataframes for next station

    Returns
    -------
    summed_df: pd.DataFrame

    """
    flag_df_1 = flag_df_1.set_index("eraqc_flag_values")
    subset = flag_df_1[~flag_df_1.index.isin(["no_flag", "total_obs_count"])]

    totals = subset.sum(numeric_only=True)
    flag_df_1.loc["total_flag"] = pd.Series(totals)

    frac = flag_df_1.loc["total_flag"] / flag_df_1.loc["total_obs_count"]
    flag_df_1.loc["frac"] = pd.Series(frac)

    rates_df = flag_df_1.loc[["frac"]]
    rates_df = rates_df.rename(index={"frac": station_name})

    rates_df = rates_df.reset_index()

    # append column of total observation count
    flag_df_1 = flag_df_1.reset_index()
    total_obs = flag_df_1[flag_df_1["eraqc_flag_values"] == "total_obs_count"].iloc[
        0, 1
    ]
    rates_df["total_obs_count"] = total_obs

    if len(flag_df_2) == 0:
        return rates_df

    else:
        rates_df_merge = pd.merge(rates_df, flag_df_2, how="outer")
        return rates_df_merge


def station_rate_tables() -> None:
    """
    Generates flag rates tables at either the station or network level.

    Parameters
    ----------
    flag_df_1: pd.DataFrame
        dataframe of previously summed station flag counts
    flag_df_2: pd.DataFrame
        flag counts dataframes for next station

    Returns
    -------
    summed_df: pd.DataFrame
    """
    network_list = ["VCAPCD", "CDEC"]

    # the function iteratively adds in flag counts to this dataframe
    flag_rate_df = []

    for network in network_list:
        # point to folder containing station flag count CSVs
        flags_prefix = f"{MERGE_DIR}/{network}/eraqc_counts_native_timestep"

        ## Merge flag counts

        # loop through all CSVs are the given level
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=flags_prefix):
            obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item.key)
            flags = pd.read_csv(obj["Body"])
            station_name = item.key.split(flags_prefix + "/")[1].split("_flag")[0]
            # the CSV is empty
            if flags.empty:
                continue
            # the CSV is not empty
            else:
                # send current dataframe and dataframe of previously summed counts to helper function
                flag_rate_df = _pairwise_rate(flags, flag_rate_df, station_name)

        ## Send final flag rates file to AWS as CSV
        csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/station_flag_rates.csv"
        # flag_rate_df.to_csv(csv_s3_filepath, index=False)
        print(f"Sending station flag rates CSV to: {csv_s3_filepath}")

    return None


def network_rate_tables() -> None:
    """
    Generates flag rates tables at either the station or network level.

    Parameters
    ----------
    flag_df_1: pd.DataFrame
        dataframe of previously summed station flag counts
    flag_df_2: pd.DataFrame
        flag counts dataframes for next station

    Returns
    -------
    summed_df: pd.DataFrame
    """

    # the function iteratively adds in flag counts to this dataframe
    flag_rate_df = []

    # point to folder containing station flag count CSVs
    flags_prefix = f"{MERGE_DIR}/per_network_flag_counts_native_timestep"

    ## Merge flag counts

    # loop through all CSVs are the given level
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=flags_prefix):
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item.key)
        flags = pd.read_csv(obj["Body"])
        station_name = item.key.split(flags_prefix + "/")[1].split("_flag")[0]
        # the CSV is empty
        if flags.empty:
            continue
        # the CSV is not empty
        else:
            # send current dataframe and dataframe of previously summed counts to helper function
            flag_rate_df = _pairwise_rate(flags, flag_rate_df, station_name)

    ## Send final flag rates file to AWS as CSV
    csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/network_flag_rates.csv"
    # flag_rate_df.to_csv(csv_s3_filepath, index=False)
    print(f"Sending station flag rates CSV to: {csv_s3_filepath}")

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
