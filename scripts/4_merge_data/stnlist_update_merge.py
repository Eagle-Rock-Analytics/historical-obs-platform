"""
stnlist_update_merge.py

This script iterates through a specified network and checks to see what stations have been successfully hourly standardized and merged,
updating the station list in the 1_raw_wx folder to reflect station availability. Error.csvs in the cleaned bucket are also parsed,
with relevant errors added to the corresponding stations if station files do not pass merge, or if the errors occur during or after the merge process.

Note that because errors.csv are parsed, very old errors.csv may want to be removed manually from AWS or thresholded below
(removing those produced during code testing). 

Functions
---------
- get_station_list: Retrieves specific network stationlist from QAQC bucket
- get_zarr_last_mod: Identifies the last modified date from a zarr 
- get_merge_stations: Retrieves list of all stations that pass the merge process
- fix_start_end_dates: Fixes two kinds of incorrect date encoding listed in the network stationlists. 
- parse_error_csv: Retrieves all processing error files for a network
- merge_qa: Update station list and save to AWS, adding merge status, time of merge pass and any relevant errors.

Intended Use
------------
Run this script after merge has been completed for a network (via pcluster run) to update the network stationlist for merge success rate tracking.
"""

import os
import sys
import pandas as pd
import xarray as xr
from io import BytesIO, StringIO
import numpy as np
import boto3
import s3fs

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from paths import BUCKET_NAME, QAQC_WX, MERGE_WX

# Set environment variables
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")


def get_station_list(network: str) -> pd.DataFrame:
    """
    Given a network name, return a pandas dataframe containing the network's station list from the QAQC bucket.
    Intentionally grabbing the QAQC version of the stationlist to retain information about whether a station was QAQC going into merge process.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    station_list : pd.DataFrame
        station list of all stations within a network
    """
    station_list = pd.read_csv(
        f"s3://wecc-historical-wx/{QAQC_WX}/{network}/stationlist_{network}_qaqc.csv"
    )
    return station_list


def get_zarr_last_mod(fn: str) -> str:
    """Identifies the "last_modified" date within a .zarr file.

    Parameters
    ----------
    fn : str
        filename of a station

    Returns
    -------
    last_mod : str
        last modified date from within zarr
    """

    # Grab only path name without extension or bucket name
    path_no_ext = fn.split(".")[0]
    path_no_bucket = path_no_ext.split(BUCKET_NAME)[-1][1:]

    # idenitfy last_modified date from metadata date within .zarr
    mod_list = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=path_no_bucket):
        mod_list.append(str(item.last_modified))

    # return most recent datetime value
    last_mod = max(mod_list)
    return last_mod


def get_merge_stations(network: str) -> pd.DataFrame:
    """
    Retrieves a list of all stations in a given network that have passed the merge process.

    This function scans the `3_qaqc_wx` folder in the AWS S3 bucket for the specified network,
    identifying which stations have successfully completed based on the presence of `.zarr` files.
    It returns a DataFrame containing the station ID, a placeholder for the timestamp when merge occurred,
    and a binary indicator ('Y/N') showing whether the `.zarr` file exists for that station.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    pd.DataFrame
        A DataFrame with three columns:
        - 'ID': Station identifier
        - 'Time_Merge': QA/QC timestamp
        - 'merged': 'Y' if a .zarr file exists (passed merge), 'N' otherwise
    """
    df = {"ID": [], "Time_Merge": [], "merged": []}  # Initialize results dictionary

    # Construct the S3 path prefix for the network inside the QAQC folder
    parent_s3_path = f"{BUCKET_NAME}/{MERGE_WX}/{network}"

    # Use s3fs to list all items under this path
    s3_fs = s3fs.S3FileSystem(anon=False)
    all_paths = s3_fs.ls(parent_s3_path)

    # Filter to only .zarr folders
    zarr_folders = [f"{path}" for path in all_paths if path.endswith(".zarr")]

    for item in zarr_folders:
        # Extract the station ID from the folder name, which is usually the last part of the path
        station_id = item.split("/")[-1].split(".")[-2].replace(".zarr", "")
        df["ID"].append(station_id)
        df["merged"].append("Y")  # merge passed

        # Retrieving "last_modified" timestamp on a .zarr requires going into group
        time_mod = get_zarr_last_mod(item)
        df["Time_Merge"].append(time_mod)

    return pd.DataFrame(df)


def fix_start_end_dates(network: str, stations: pd.DataFrame) -> pd.DataFrame:
    """
    Fixes two kinds of incorrect date encoding listed in the network stationlists.
    Kind 1: Missing start and end dates listed in the stationlist for specific networks.
    Issue originated from source network (not provided via raw station list).
    Impacted networks: MARITIME, NDBC, CW3E
    #! Handful of individual stations in RAWS, HADS, CWOP (handle separately)

    Kind 2: End date incorrectly encoded as 2100-12-31. Issue originated from source network (known issue).
    Impacted networks: SCAN, SNOTEL

    Parameters
    ----------
    stations : pd.DataFrame
        original stationlist

    Returns
    -------
    fixed_stns : pd.DataFrame
        stationlist with corrected start and end dates

    Notes
    -----
    To identify the correct start/end date, the station file has to be opened. For the "Kind 2" error,
    we are also "resetting" the start date for ease of computation. Dates should be identical.
    #! Nice to have: build in check for identical dates
    """

    # networks with incorrect date encoding in station list
    wrong_date = ["MARITIME", "NDBC", "CW3E", "SCAN", "SNOTEL"]

    if network in wrong_date:
        print(
            "Network has known issue with start/end date coverage. This may take some time to correct."
        )
        for id in stations["ERA-ID"]:
            print(f"Checking start/end date encoding for {id}...")

            # identify correct start/end date from station timestamps
            try:
                ds = xr.open_zarr(
                    f"s3://{BUCKET_NAME}/{MERGE_WX}/{network}/{id}.zarr",
                    consolidated=False,
                )
            except Exception as e:
                continue

            correct_start = str(ds.time.values[0])
            correct_end = str(ds.time.values[-1])
            # save memory
            ds.close()

            # set correct start/end date time to station list
            try:
                stations.loc[stations["ERA-ID"] == id, "start-date"] = correct_start
                stations.loc[stations["ERA-ID"] == id, "end-date"] = correct_end

            except Exception as e:
                print(f"Issue setting correct start / end dates for {id}... {e}")

    else:
        print("Network has no known start/end date issues.")

    # put these columns at specific index (so they're not at the end)
    try:
        start = stations.pop("start-date")
        end = stations.pop("end-date")
        stations.insert(4, "start-date", start)
        stations.insert(5, "end-date", end)
    except:
        print(f"Issue resetting index on start/end")

    return stations


def parse_error_csv(network: str) -> pd.DataFrame:
    """
    Given a network name, return a pandas dataframe containing all errors reported for the network in the merge stage.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    errordf : pd.Dataframe
        dataframe of all errors produced during QAQC
    """

    # List to store all non-empty error DataFrames found
    errordf = []

    # Define the path prefix in the S3 bucket for error CSVs
    errors_prefix = f"{MERGE_WX}/{network}/merge_errs/errors"

    # Loop over all objects in the specified S3 prefix
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=errors_prefix):
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item.key)
        errors = pd.read_csv(obj["Body"])
        if errors.empty:  # If file empty
            continue
        else:
            errors = errors[["File", "Time", "Error"]]
            errordf.append(errors)  # Add to list of error records
    if not errordf:  # If no errors in cleaning
        return pd.DataFrame()
    else:
        errordf = pd.concat(errordf)

        # Drop duplicate error messages for the same file
        errordf = errordf.drop_duplicates(subset=["File", "Error"])

        # Remove general network-wide errors, keeping only station-specific ones
        errordf = errordf[errordf.File != "Whole network"]
        return errordf


def merge_qa(network: str):
    """
    Update station list and save to AWS, adding merge status, time of merge pass and any relevant errors.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    None
    """
    # Fixing capitalization issues
    if network == "otherisd":
        network = "OtherISD"

    # Call functions
    stations = get_station_list(network)  # grabs stationlist_qaqcd
    merge_ids = get_merge_stations(network)  # grabs stations that pass merge
    errors = parse_error_csv(network)  # grabs station error files

    if merge_ids.empty:
        print(
            "No merged files for this network. Please run the relevant merge script and try again."
        )
        exit()

    # Fix start/end date issues
    stations = fix_start_end_dates(network, stations)

    # Join qaqc'd columns to column list
    stations = stations.merge(merge_ids, left_on="ERA-ID", right_on="ID", how="outer")
    if "ID_y" in stations.columns:
        # Make binary merged column
        stations["merged"] = np.where(stations.ID_y.isna(), "N", "Y")
        # Drop ID column
        stations = stations.drop(["ID_x", "ID_y"], axis=1)
    else:
        # Make binary merged column
        stations["merged"] = np.where(stations.ID.isna(), "N", "Y")
        # Drop ID column
        stations = stations.drop("ID", axis=1)

    # Move Time_Merge to last
    s = stations.pop("Time_Merge")
    stations = pd.concat([stations, s], axis=1)

    # Add errors to column by station - only add error if error occurred at or after file clean, if file cleaned.
    stations["Errors_Merge"] = np.nan

    # Remove any NAs from ERA-ID
    stations = stations.loc[stations["ERA-ID"].notnull()]

    # Get list of station IDs
    ids = [id.split("_")[-1] for id in stations["ERA-ID"].tolist()]

    if errors.empty:  # If no errors, stop here.
        pass

    else:
        # Add relevant ID to errors csv
        errors["ID"] = np.nan
        errors.reset_index(inplace=True, drop=True)
        errors["Time"] = pd.to_datetime(errors["Time"], format="%Y%m%d%H%M", utc=True)

        for index, row in errors.iterrows():
            id = []
            id = [x for x in ids if x in row["File"]]
            if id:
                errors.loc[index, "ID"] = network + "_" + id[-1]

        # For each station
        for index, row in stations.iterrows():
            error_sta = errors.loc[errors.ID == row["ERA-ID"]]
            if error_sta.empty:
                # if no errors for station
                continue
            else:
                if not pd.isnull(row["Time_Merge"]):
                    # If file cleaned
                    # Only keep errors from merge at or after time of merge
                    error_sta = error_sta.loc[
                        (error_sta.Time >= row["Time_Merge"]) | (error_sta.Time.isna()),
                        :,
                    ]

                if len(error_sta) == 1:
                    stations.loc[index, "Errors_Merge"] = error_sta["Error"].values[0]
                elif len(error_sta) > 1:
                    values = [
                        f"{x.File}: {x.Error}" for index, x in error_sta.iterrows()
                    ]
                    value = " ".join(values)
                    stations.loc[index, "Errors_Merge"] = value

    # Print summary
    merge_Y = stations["merged"].value_counts()["Y"]
    try:
        merge_N = stations["merged"].value_counts()["N"]
    except:
        # for some reason for networks with no failed merges, this is causing problems
        merge_N = 0

    if "Y" in stations["merged"].values:
        if "N" not in stations["merged"].values:
            # order is important here, if no "N" is present in a merged network, it will bark without this
            print(
                f"Station list updated for {network} stations that pass merge. All stations pass: {merge_Y} stations."
            )
        else:
            print(
                f"Station list updated for {network} stations that pass merge. {merge_Y} stations passed merge, {merge_N} stations were not merged."
            )
    else:
        print(
            f"Station list updated for {network} stations. No stations successfully pass merge. {merge_N} stations do not yet pass merge."
        )

    # Save station file to cleaned bucket
    new_buffer = StringIO()
    stations.to_csv(new_buffer, index=False)
    content = new_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=f"{MERGE_WX}/{network}/stationlist_{network}_merge.csv",
    )


if __name__ == "__main__":
    import argparse

    ALL_NETWORKS = [
        "ASOSAWOS",
        "CAHYDRO",
        "CIMIS",
        "CW3E",
        "CDEC",
        "CNRFC",
        "CRN",
        "CWOP",
        "HADS",
        "HNXWFO",
        "HOLFUY",
        "HPWREN",
        "LOXWFO",
        "MAP",
        "MTRWFO",
        "NCAWOS",
        "NOS-NWLON",
        "NOS-PORTS",
        "otherisd",
        "RAWS",
        "SGXWFO",
        "SHASAVAL",
        "VCAPCD",
        "MARITIME",
        "NDBC",
        "SCAN",
        "SNOTEL",
        "VALLEYWATER",
    ]

    parser = argparse.ArgumentParser(
        description="Update station list with merge status for a given network."
    )
    parser.add_argument(
        "network",
        type=str,
        help=f"Network name (e.g. ASOSAWOS, CWOP, HADS). Valid options: {', '.join(ALL_NETWORKS)}",
    )
    args = parser.parse_args()
    merge_qa(args.network)
