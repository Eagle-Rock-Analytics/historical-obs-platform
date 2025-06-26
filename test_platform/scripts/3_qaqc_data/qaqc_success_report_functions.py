"""
qaqc_success_report_functions.py

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
Import into qaqc_flag_counts_sum.ipynb to generate information for the QAQC success report.
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
bucket_name = "wecc-historical-wx"
stations_csv_path = f"s3://{bucket_name}/2_clean_wx/temp_clean_all_station_list.csv"
qaqc_dir = "3_qaqc_wx"
merge_dir = "4_merge_wx"


# -----------------------------------------------------------------------------
def _pairwise_sum(flag_df_1, flag_df_2) -> pd.DataFrame:
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
    if len(flag_df_1) == 0:
        return flag_df_2
    else:
        total_df = pd.concat([flag_df_1, flag_df_2])

        summed_df = total_df.groupby("eraqc_flag_values", as_index=False).sum(
            numeric_only=True
        )
        return summed_df


# -----------------------------------------------------------------------------
def _network_format_table(
    summed_counts: pd.DataFrame, flag_table: pd.DataFrame
) -> pd.DataFrame:
    """
    A helper function that formats the network-level counts tables

    Parameters
    ----------
    summed_counts: pd.DataFrame
        dataframe of summed station flag counts
    flag_table: pd.DataFrame
        flag counts dataframes for next station

    Returns
    -------
    final_format: pd.DataFrame

    """
    ## Format flag meanings df
    flag_table = flag_table.rename(columns={"Flag_value": "eraqc_flag_values"})

    ## Format summed counts df

    # remove the ".0" from the flag values
    summed_counts["eraqc_flag_values"] = summed_counts["eraqc_flag_values"].str.replace(
        ".0", "", regex=True
    )

    # convert flag value strings to integers
    summed_counts["eraqc_flag_values"] = summed_counts["eraqc_flag_values"].apply(
        lambda x: int(x) if x not in ["no_flag", "total_obs_count"] else x
    )

    ## Merge the the counts and flag meanings dataframes
    merged_dfs = summed_counts.merge(
        flag_table, on="eraqc_flag_values", how="outer"
    ).fillna(0)

    ## Format final dataframe

    # order by flag value, in descending numerical order
    final_format = (
        merged_dfs.groupby(
            merged_dfs.eraqc_flag_values.apply(type) != str, group_keys=True
        )
        .apply(lambda g: g.sort_values("eraqc_flag_values"))
        .reset_index(drop=True)
    )

    # move string flag value entries to the bottom
    final_format = final_format.loc[
        pd.to_numeric(final_format["eraqc_flag_values"], errors="coerce")
        .sort_values()
        .index
    ]

    # convert all counts to integers
    final_format = final_format.applymap(
        lambda x: int(x) if not isinstance(x, str) else x
    )

    return final_format


# -----------------------------------------------------------------------------
def _total_format_table(summed_counts: pd.DataFrame) -> pd.DataFrame:
    """
    A helper function that formats the final, 'total' counts table

    Parameters
    ----------
    summed_counts: pd.DataFrame
        dataframe of summed station flag counts

    Returns
    -------
    final_format: pd.DataFrame

    """

    # convert flag value strings to integers

    summed_counts = summed_counts.applymap(
        lambda x: int(x) if not isinstance(x, str) else x
    )

    ## Format final dataframe

    # order by flag value, in descending numerical order
    final_format = (
        summed_counts.groupby(
            summed_counts.eraqc_flag_values.apply(type) != str, group_keys=True
        )
        .apply(lambda g: g.sort_values("eraqc_flag_values"))
        .reset_index(drop=True)
    )

    # move string flag value entries to the bottom
    final_format = final_format.loc[
        pd.to_numeric(final_format["eraqc_flag_values"], errors="coerce")
        .sort_values()
        .index
    ]

    return final_format


# -----------------------------------------------------------------------------
def network_sum_flag_counts(network: str, timestep: str) -> None:
    """
    Sums all station QAQC flag counts in a network for a given timestep (hourly or native) and sends to AWS.
    These counts are used to generate statistics for the QAQC success report.

    Parameters
    ----------
    network: str
        network name
    timestep: str
        if set to 'hourly', merge all hourly QAQC flag count tables
        if set to 'native', merge all native timestep QAQC flag count tables

    Returns
    -------
    None

    """
    ## Setup

    # read in flag meanings CSV

    flag_meanings = pd.read_csv("era_qaqc_flag_meanings.csv")

    # only run for a valid "timestep" input
    if timestep not in ("hourly", "native"):
        print("invalid timestep: ", timestep)
        return None

    # the function iteratively adds in flag counts to this dataframe
    summed_counts_df = []

    # point to folder containing station flag count CSVs
    flags_prefix = f"{merge_dir}/{network}/eraqc_counts_{timestep}_timestep"

    ## Merge flag counts

    # loop through all CSVs are the given level
    for item in s3.Bucket(bucket_name).objects.filter(Prefix=flags_prefix):
        obj = s3_cl.get_object(Bucket=bucket_name, Key=item.key)
        flags = pd.read_csv(obj["Body"])
        # the CSV is empty
        if flags.empty:
            continue
        # the CSV is not empty
        else:
            # send current dataframe and dataframe of previously summed counts to helper function
            summed_counts_df = _pairwise_sum(summed_counts_df, flags)

    counts_final = _network_format_table(summed_counts_df, flag_meanings)

    ## Send final counts file to AWS as CSV

    csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/per_network_flag_counts_{timestep}_timestep/{network}_flag_counts_{timestep}_timestep.csv"
    counts_final.to_csv(csv_s3_filepath, index=False)
    print(f"Sending summed counts dataframe for {network} to: {csv_s3_filepath}")

    return None


# -----------------------------------------------------------------------------
def total_sum_flag_counts(timestep: str) -> None:
    """
    Sums all network-level QAQC flag counts for a given timestep (hourly or native) and sends to AWS.
    These counts are used to generate statistics for the QAQC success report.

    Parameters
    ----------
    timestep: str
        if set to 'hourly', merge all hourly QAQC flag count tables
        if set to 'native', merge all native timestep QAQC flag count tables

    Returns
    -------
    None

    """
    ## Setup

    # only run for a valid "timestep" input
    if timestep not in ("hourly", "native"):
        print("invalid timestep: ", timestep)
        return None

    # the function iteratively adds in flag counts to this dataframe
    summed_counts_df = []

    # point to folder containing network-level flag count CSVs
    flags_prefix = f"{merge_dir}/per_network_flag_counts_{timestep}_timestep"

    ## Merge flag counts

    # loop through all networks CSVs
    for item in s3.Bucket(bucket_name).objects.filter(Prefix=flags_prefix):
        obj = s3_cl.get_object(Bucket=bucket_name, Key=item.key)
        flags = pd.read_csv(obj["Body"])
        # the CSV is empty
        if flags.empty:
            continue
        # the CSV is not empty
        else:
            # send current dataframe and dataframe of previously summed counts to helper function
            print(f"summing for {item.key}")
            summed_counts_df = _pairwise_sum(summed_counts_df, flags)

    # format final table
    final_table = _total_format_table(summed_counts_df)

    ## Send final counts file to AWS as CSV
    if len(summed_counts_df) == 0:
        return None
    else:
        csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/total_flag_counts_{timestep}_timestep.csv"
        final_table.to_csv(csv_s3_filepath, index=False)
        print(f"Sending final summed counts dataframe for to: {csv_s3_filepath}")

        return None


# -----------------------------------------------------------------------------
def generate_station_tables(timestep: str) -> None:
    """
    Runs network_sum_flag_counts() for every network.

    Parameters
    ----------
    timestep: str
        if set to 'hourly', merge all hourly QAQC flag count tables
        if set to 'native', merge all native timestep QAQC flag count tables

    Returns
    -------
    None

    """
    # record start time
    start_time = time.time()

    # only run for a valid "timestep" input
    if timestep not in ("hourly", "native"):
        print("invalid timestep: ", timestep)
        return None

    station_list = pd.read_csv(stations_csv_path)
    network_list = station_list["network"].unique()

    for network in network_list:
        network_sum_flag_counts(network, timestep)

    # record end time
    end_time = time.time()

    # output time elapsed
    time_elapsed = (end_time - start_time) / 60
    print(f"{time_elapsed} minutes")

    return None


if __name__ == "__main__":
    main()
