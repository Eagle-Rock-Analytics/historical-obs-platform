"""
qaqc_concatenate_stations.py

This script identifies weather stations with identical locations within specified networks,
flags groups of stations for concatenation, and merges their data with special handling for:

- Groups containing more than two stations
- Stations with overlapping time periods (keeping newer data in overlaps)

After concatenation, the original input datasets are moved or renamed in AWS S3
to maintain data provenance without deletion.

Key steps
---------
1. Identify stations to concatenate based on location
2. Concatenate station data with overlap handling
3. Manage original datasets post-concatenation in AWS

Functions
---------
- main: Orchestrates the full workflow: identifies candidate stations, performs concatenation, and exports results.
- concatenation_check: Flags stations at identical lat/lon locations and assigns group IDs.
- apply_concat_check: Applies `concatenation_check()` to each network and uploads results to AWS.
- concatenate_stations: Coordinates concatenation for each group of stations and handles export.
- _df_concat: Concatenates two station datasets, handling overlaps if present.
- _overlap_concat: Helper used by `_df_concat()` to resolve overlapping time periods.
- _more_than_2: Iteratively merges more than two stations in a group via pairwise concatenation.
- _concat_export_help: Prepares and formats concatenated dataset for export, including metadata.
- concat_export: Exports final dataset to S3.
- _rename_file: Renames original station datasets in S3 to mark them as deprecated.

Intended Use
------------
This script only should be run once to concatenate the identified stations into a single record, rather than split across many stations.
This script will fail on any additional runs because the input stations will be removed from the QAQC directory, causing failure.
Regeneration of the original input QAQC'd stations is required for a full re-run.
"""

import datetime
import boto3
import pandas as pd
import xarray as xr
from io import StringIO
from time import time
from QAQC_pipeline import qaqc_ds_to_df

# Silence warnings
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# A list of networks to be checked for concatenation
target_networks = ["ASOSAWOS", "MARITIME"]

# AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")

# AWS buckets
BUCKET_NAME = "wecc-historical-wx"
QAQC_DIR = "3_qaqc_wx/"


def main():

    print("Starting script qaqc_concatenate_stations.py")
    t0 = time()  # start time

    # Identify MARITIME and ASOSAWOS stations with identical lat/lon,
    # assign a unique subset ID to each group, and upload a DataFrame
    # of ERA-IDs and subset IDs to AWS for each network.
    apply_concat_check(target_networks)

    # Concatenates stations for a given network using helper functions,
    # then exports the final result to s3 and returns the list of concatenated ERA-IDs.
    for network_name in target_networks:
        final_concat_list = concatenate_stations(network_name)

    # Print elapsed time
    ttime = time() - t0
    print(f"Script complete. Elapsed time: {ttime:.2f} s.\n")


# Helper functions
def concatenation_check(station_list: list) -> pd.DataFrame:
    """
    This function flags stations that need to be concatenated.

    Rules
    ------
    1.) Stations are flagged if they have identical latitudes and longitudes

    Parameters
    ------
    station_list: list of str
        list of station information

    Returns
    -------
    if success: returns input station list with a flag column assigning an integer to each group of repeat latitudes and longitudes
    if failure: None
    """
    # Flag stations with identical latitudes and longitudes, then assign each group a unique integer
    # List of possible variable names for longitudes and latitudes
    lat_lon_list = [
        "LAT",
        "LON",
        "latitude",
        "longitude",
        "LATITUDE",
        "LONGITUDE",
        "lat",
        "lon",
    ]
    # Extract the latitude and longitude variable names from the input dataframe
    lat_lon_cols = [col for col in station_list.columns if col in lat_lon_list]

    # Generate column flagging duplicate latitudes and longitudes
    station_list["concat_subset"] = station_list.duplicated(
        subset=lat_lon_cols, keep=False
    )
    # within each group of identical latitudes and longitudes, assign a unique integer
    station_list["concat_subset"] = (
        station_list[station_list["concat_subset"] == True]
        .groupby(lat_lon_cols)
        .ngroup()
    )

    # Order station list by flag
    concat_station_list = station_list.sort_values("concat_subset")

    # Keep only flagged stations
    concat_station_list = concat_station_list[
        ~concat_station_list["concat_subset"].isna()
    ]

    # Convert flags to integers - this is necessary for the final concatenation step
    concat_station_list["concat_subset"] = concat_station_list["concat_subset"].astype(
        "int32"
    )
    # Now keep only the ERA-ID and flag column
    era_id_list = ["ERA-ID", "era-id"]
    era_id_col = [col for col in station_list.columns if col in era_id_list]
    concat_station_list = concat_station_list[era_id_col + ["concat_subset"]]

    # Standardize ERA id to "ERA-ID" (this is specific to Valleywater stations)
    if "era-id" in era_id_col:
        concat_station_list.rename(columns={"era-id": "ERA-ID"}, inplace=True)

    return concat_station_list


def apply_concat_check(
    station_names_list: list,
):
    """
    This function applies the conatenation check to a list of target stations.
    It then upload a csv containing the ERA IDs and concatenation subset ID for
    all identified stations in a network.

    Parameters
    ----------
    station_names_list : list of str
        List of the target station names.

    Returns
    -------
    if success: uploads list of stations to be concatenated to AWS
    if failure: None
    """

    final_list = pd.DataFrame([])
    for station in station_names_list:

        # Need to use the QAQC stationlist to implement the QAQC binary Y/N flag
        # Import station list of target station
        key = f"3_qaqc_wx/{station}/stationlist_{station}_qaqc.csv"
        station_list = pd.read_csv(f"s3://{BUCKET_NAME}/{key}")

        # subset for only stations that passed QA/QC
        station_list_qc = station_list.loc[station_list["QAQC"] == "Y"]

        # Apply concatenation check
        concat_list = concatenation_check(station_list_qc)

        # Rename the flags for each subset to <station>_<subset number>
        concat_list["concat_subset"] = (
            station + "_" + concat_list["concat_subset"].astype(str)
        )

        # Append to final list of stations to concatenate
        final_list = pd.concat([final_list, concat_list])

        # Upload to QAQC directory in AWS
        new_buffer = StringIO()
        final_list.to_csv(new_buffer, index=False)
        content = new_buffer.getvalue()

        # the csv is stored in each station folder within 3_qaqc_wx
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=QAQC_DIR + station + f"/concat_list_{station}.csv",
        )

    return None


def _overlap_concat(df_new: pd.DataFrame, df_old: pd.DataFrame) -> pd.DataFrame:
    """
    Handles the cases in which there is overlap between the two input stations

    Rules
    ------
    1.) concatenation: keep the newer station data in the time range in which both stations overlap

    Parameters
    ----------
    df_new: pd.DataFrame
        weather station network
    df_old: pd.DataFrame
        weather station network

    Returns
    -------
    if success: returns pd.DataFrame
    if failure: None
    """

    # identify where there is overlap in timestamps, and keep from newer station data
    df_overlap = df_new[df_new["time"].isin(df_old["time"])]
    print(f"Length of overlapping period: {len(df_overlap)}")

    # Split datframes into subsets
    # Remove data in time overlap between old and new
    df_old_cleaned = df_old[~df_old["time"].isin(df_overlap["time"])]
    df_new_cleaned = df_new[~df_new["time"].isin(df_overlap["time"])]

    # Concatenate subsets
    df_concat = pd.concat([df_old_cleaned, df_overlap, df_new_cleaned])

    return df_concat


def _df_concat(
    df_1: pd.DataFrame,
    df_2: pd.DataFrame,
    attrs_1: dict,
    attrs_2: dict,
    var_attrs_1: dict,
    var_attrs_2: dict,
) -> tuple[pd.DataFrame, str, str, str, dict]:
    """
    Performs concatenation of input datasets, handling two cases
        1.) temporal overlap between the datasets
        2.) no temporal overlap

    Rules
    ------
    1.) concatenation: keep the newer station data in the time range in which both stations overlap

    Parameters
    ----------
    df_1: pd.DataFrame
        station data
    df_2: pd.DataFrame
        dtation data
    attrs_1: list of str
        attributes of df_1
    attrs_2: list of str
        attributes of df_2
    var_attrs_1: list of str
        variable attributes of ds_1
    var_attrs_2: list of str
        variable attributes of ds_2

    Returns
    -------
    if success:
        df_concat: pd.DataFrame
        stn_n_to_keep: str
        stn_n_to_drop: str
        attrs_new: dict
        var_attrs_new: dict
    if failure: None
    """

    # determine which dataset is older
    if df_1["time"].max() < df_2["time"].max():
        # if df_1 has an earlier end tiem than df_2, then d_2 is newer
        # we also grab the name of the newer station in this step, for use later
        df_new = df_2
        attrs_new = attrs_2
        var_attrs_new = var_attrs_2
        df_old = df_1

    else:
        df_new = df_1
        attrs_new = attrs_1
        var_attrs_new = var_attrs_1
        df_old = df_2

    stn_n_to_keep = df_new["station"].unique()[0]
    stn_n_to_drop = df_old["station"].unique()[0]
    print(f"\nStation will be concatenated and saved as: {stn_n_to_keep}")

    # now set things up to determine if there is temporal overlap between df_new and df_old
    df_overlap = df_new[df_new["time"].isin(df_old["time"])]

    # If there is no overlap between the two time series, just concatenate
    if len(df_overlap) == 0:
        print("No overlap!")
        df_concat = pd.merge(df_old, df_new, how="outer")

    # If overlap exists, split into subsets and concatenate
    else:
        print("There is overlap")
        df_concat = _overlap_concat(df_old, df_new)

    # Reset station name to be the newer station
    df_concat["station"] = stn_n_to_keep

    return df_concat, stn_n_to_keep, stn_n_to_drop, attrs_new, var_attrs_new


def _more_than_2(
    network_name: str, stns_to_pair: pd.DataFrame
) -> tuple[pd.DataFrame, dict, dict, list]:
    """
    Performs pairwise concatenation on subsets of more than two stations flagged for concatenation

    Rules
    ------
    1.) concatenation: keep the newer station data in the time range in which both stations overlap

    Parameters
    ----------
    network_name: string
        weather station network
    stns_to_pair: pd.DataFrame
        dataframe of the input station names

    Returns
    -------
    if success:
        df_concat: pd.DataFrame
        station_names: dict
        attrs_new: dict
        datasets: list of xr.Dataset
    if failure: None
    """

    print(f"Concatenating the following stations: {stns_to_pair}")

    # Load datasets into a list
    datasets = [
        xr.open_zarr(
            f"s3://wecc-historical-wx/3_qaqc_wx/{network_name}/{stn}.zarr",
            consolidated=True,
        )
        for stn in stns_to_pair["ERA-ID"]
    ]

    # Sort datasets by their max 'time'
    datasets_sorted = sorted(datasets, key=lambda ds: ds["time"].max())

    # Store station names, in order from oldest to newest
    names = [ds.coords["station"].values[0] for ds in datasets_sorted]

    print(f"Newest station: {names[-1]}")

    # Setup for the while loop
    ds_1 = datasets_sorted[0]
    df_1, MultiIndex_1, attrs_1, var_attrs_1, era_qc_vars_1 = qaqc_ds_to_df(ds_1)
    i = 0
    end = len(datasets_sorted) - 1

    while i < end:

        print("iteration:", i)

        ds_2 = datasets_sorted[i + 1]
        df_2, MultiIndex_2, attrs_2, var_attrs_2, era_qc_vars_2 = qaqc_ds_to_df(ds_2)

        # Send to helper function for concatenation
        df_concat, stn_n_to_keep, stn_n_to_drop, attrs_new, var_attrs_new = _df_concat(
            df_1, df_2, attrs_1, attrs_2, var_attrs_1, var_attrs_2
        )

        df_1 = df_concat
        attrs_1 = attrs_new
        var_attrs_1 = var_attrs_new

        i += 1

    # Construct station names list, for updating attributes
    newest_station = names[-1]  # Get last station name from station name list
    older_stations = ", ".join(
        names[:-1]
    )  # Create a string containing all older station names
    station_names = {"station_name_new": newest_station, "old_stations": older_stations}

    print("Progressive concatenation for 2+ stations is complete.")

    new_column = [newest_station] * len(df_concat)
    df_concat["station"] = new_column

    return df_concat, station_names, attrs_new, var_attrs_new, datasets


def _concat_export_help(
    df_concat: pd.DataFrame,
    final_concat_list: list[str],
    network_name: str,
    attrs_new: dict,
    var_attrs_new: dict,
    station_names: dict,
) -> tuple[pd.DataFrame, str]:
    """
    Prepares the final concatenated dataset for export by
    - updating the attributes and
    - converting one of the mulit-index levels to the correct datatype
    then exports the dataset to AWS

    Rules
    ------
    1.) retains the name of the newest station

    Parameters
    ----------
    df_concat: pd.DataFrame
        dataframe of concatenated dataframes
    final_concat_list: list[str]
        list of stations that have been concatenated
    network_name: str
        weather station network
    attrs_new: dict
        attributes of newer dataframe that was input to concatenation
    var_attrs_new: dict
        attributes of newer dataframe that was input to concatenation
    station_names: dict
        library of station names, including the single new station name and a string of all the older station names

    Returns
    -------
    if successful, exports dataset of concatenated dataframes to AWS
    if failure, returns None
    """

    # Prepare concatenated dataset for export
    # Final sort on time
    df_concat = df_concat.sort_values(by="time")

    # Delete unnecessary columns and set index
    df_concat = df_concat.drop(["hour", "day", "month", "year", "date"], axis=1)
    df_to_export = df_concat.set_index(["station", "time"])

    # Convert concatenated dataframe to dataset
    ds_concat = df_to_export.to_xarray()

    # Convert datatype of station coordinate
    ds_concat.coords["station"] = ds_concat.coords["station"].astype("<U20")

    # Include past attributes
    for i in attrs_new:
        ds_concat.attrs[i] = attrs_new[i]

    # for var, value in var_attrs_new:
    #     ds_concat[var] = ds_concat[var].assign_attrs(value)

    # Update 'history' attribute
    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")
    ds_concat.attrs["history"] = (
        ds_concat.attrs["history"] + f" \nstation_matching.ipynb run on {timestamp} UTC"
    )

    # Update 'comment' attribute
    ds_concat.attrs["comment"] = (
        "Intermediary data product. This data has been subjected to cleaning, QA/QC, but may not have been standardized."
    )

    # Extract old and new station names from name dictionary
    station_name_new = station_names["station_name_new"]
    station_name_old = station_names["old_stations"]

    # Add new qaqc_files_merged attribute
    ds_concat.attrs["qaqc_files_merged"] = (
        f"{station_name_old}, {station_name_new} merged. Overlap retained from newer station data."
    )

    return ds_concat, station_name_new


def concat_export(ds: xr.Dataset, network_name: str, station_name_new: str):
    """
    Exports the concatenated dataset to s3.

    Parameters
    ----------
    ds : xr.Dataset
        qaqc'd and concatenated dataset
    network_name : str
        name of network
    station_name_new : str
        new name of the station to export

    Returns
    -------
    None
    """

    # Export
    export_url = (
        f"s3://wecc-historical-wx/3_qaqc_wx/{network_name}/{station_name_new}.zarr"
    )
    print(f"Exporting concatenated dataset... {export_url}")

    ## TURN OFF FOR TESTING MODE -- WHEN READY FOR EXPORT
    ds.to_zarr(export_url, mode="w")

    return None


def _rename_file(
    stn: xr.Dataset,
    network: str,
):
    """
    Renames a given file in AWS by copying it over into the new name and then deleting the old file

    Parameters
    ----------
    ds: xr.Dataset
        qaqc'd dataset to rename in s3
    network: str
        weather station network name

    Returns
    -------
    if success: None
    if failure: None
    """

    try:
        # pull station name
        old_name = stn.station.values[0]

        old_url = f"s3://{BUCKET_NAME}/3_qaqc_wx/{network}/{old_name}.zarr"
        print(f"Original file name: {old_url}")

        # build new name with _c identifier
        new_name = f"{old_name}_c"
        new_url = f"s3://{BUCKET_NAME}/3_qaqc_wx/{network}/{new_name}.zarr"

        # Rename using input and save to new path
        stn.to_zarr(new_url, mode="w")

        # Delete older version of the file using just old stn, no s3 bucket portion
        ## TURN OFF FOR TESTING MODE -- WHEN READY FOR EXPORT
        for obj in s3.Bucket(BUCKET_NAME).objects.filter(
            Prefix=f"3_qaqc_wx/{network}/{old_name}.zarr/"
        ):
            obj.delete()

        print(
            f"File {old_name} renamed to {new_name}, exported, and the redundant older file deleted from s3. "
        )

    except Exception as e:
        print(f"Error renaming file: {e}")


def concatenate_stations(network_name: str) -> None:
    """
    Coordinates the concatenation of input datasets and exports the final concatenated dataset.
    Also returns a list of the ERA-IDs of all stations that are concatenated.

    Parameters
    ----------
    network_name: string
        weather station network

    Returns
    -------
    final_concat_list: list
        List of ERA-IDs of all stations that are concatenated

    Notes
    -----
    Uses the following helper functions
        _df_concat(): concatenates two dataframes
        _overlap_concat(): used by _df_concat() to concatenate two stations with overlapping time ranges
        _more_than_2(): handles subsets with more than two stations, passing pairs to _df_concat() iteratively
        _concat_export_help(): formats and exports concatenated dataframe

    """
    # Initiate empty list, to which we will iteratively add the ERA-IDs of stations that are concatenated
    final_concat_list = []

    # Read in full concat station list
    print(network_name)
    concat_list = pd.read_csv(
        f"s3://{BUCKET_NAME}/3_qaqc_wx/{network_name}/concat_list_{network_name}.csv"
    )

    # Identify stns within designated network
    concat_by_network = concat_list.loc[
        concat_list.concat_subset.str.contains(network_name)
    ]

    unique_pair_names = concat_by_network.concat_subset.unique()
    # For MARITIME, remove these stations becuase they're actually separate stations
    if network_name == "MARITIME":
        unique_pair_names = unique_pair_names[1:]
    else:
        pass

    print(
        f"There are {len(concat_by_network)} stations to be concatenated into {len(unique_pair_names)} station pairs within {network_name}..."
    )
    print(unique_pair_names)

    # Set up pairs
    for pair in unique_pair_names:
        print("\n", pair)
        # Pull out stations corresponding to pair name
        stns_to_pair = concat_by_network.loc[concat_by_network.concat_subset == pair]

        if len(stns_to_pair) == 2:  # 2 stations to concat together
            print(stns_to_pair)

            # Import this subset of datasets and convert to dataframe
            stnpair0 = stns_to_pair.iloc[0]["ERA-ID"]
            stnpair1 = stns_to_pair.iloc[1]["ERA-ID"]
            url_1 = f"s3://{BUCKET_NAME}/3_qaqc_wx/{network_name}/{stnpair0}.zarr"
            url_2 = f"s3://{BUCKET_NAME}/3_qaqc_wx/{network_name}/{stnpair1}.zarr"

            print(f"Retrieving.... {url_1}")
            print(f"Retrieving.... {url_2}")
            ds_1 = xr.open_zarr(url_1)
            ds_2 = xr.open_zarr(url_2)

            # Add stns to ds_list for export control
            ds_list_to_export = [ds_1, ds_2]

            # Convert to dataframes with corresponding information
            df_1, MultiIndex_1, attrs_1, var_attrs_1, era_qc_vars_1 = qaqc_ds_to_df(
                ds_1
            )
            df_2, MultiIndex_2, attrs_2, var_attrs_2, era_qc_vars_2 = qaqc_ds_to_df(
                ds_2
            )

            # Send to helper function for concatenation
            df_concat, stn_n_to_keep, stn_n_to_drop, attrs_new, var_attrs_new = (
                _df_concat(df_1, df_2, attrs_1, attrs_2, var_attrs_1, var_attrs_2)
            )
            print(f"Length of input dataframes: {len(df_1)}, {len(df_2)}")

            # Construct dictionary of old and new station names
            station_names = {
                "station_name_new": stn_n_to_keep,
                "old_stations": stn_n_to_drop,
            }

        else:
            # If there are more than 2 stations in the given subset, pass to _more_than_2()
            print(f"More than 2 stations within a subset: {len(stns_to_pair)}.")
            df_concat, station_names, attrs_new, var_attrs_new, ds_list_to_export = (
                _more_than_2(network_name, stns_to_pair)
            )

        print(f"Length of new dataframe: {len(df_concat)}")

        # Add concatenated station names to station name list
        final_concat_list.extend(stns_to_pair["ERA-ID"].tolist())

        # Send concatenated dataframe to helper function to package for export
        ds_concat, station_name_new = _concat_export_help(
            df_concat,
            final_concat_list,
            network_name,
            attrs_new,
            var_attrs_new,
            station_names,
        )

        # Rename input files
        print(f"Renaming {len(stns_to_pair)} input stations...")
        for stn in ds_list_to_export:
            _rename_file(stn, network_name)

        # Export final dataset
        concat_export(ds_concat, network_name, station_name_new)

    return final_concat_list


if __name__ == "__main__":
    main()
