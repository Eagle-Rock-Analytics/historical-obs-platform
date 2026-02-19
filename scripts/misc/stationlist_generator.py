"""
stationlist_generator.py

Generates the "all network" station list, using the "per network" station lists. Functionality to generate for each stage of development 
(pull, clean, qa/qc, and merge). To replace the "station list" generation in the figure notebooks. 

Functions
---------
- get_station_list_paths: Builds a dictionary of all stationlists within a sublevel of s3 bucket
- retrieve_and_concat_stnlists: Retrieves the stationlists using the paths dictionary, and concats together. 
- export_stationlist: Helper function to export final csv file
- generate_stationlist: Core processing function

Intended Use
-------------
Run this script after the corresponding "stnlist_update" script has been completed for all networks in the respective stage of development. 
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import boto3
import pandas as pd
import numpy as np
import datetime
from io import BytesIO, StringIO

from paths import BUCKET_NAME, RAW_WX, CLEAN_WX, QAQC_WX, MERGE_WX

PULL_DIR = f"{RAW_WX}/"
CLEAN_DIR = f"{CLEAN_WX}/"
QAQC_DIR = f"{QAQC_WX}/"
MERGE_DIR = f"{MERGE_WX}/"

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")

CLEANED_VARS = [
    "tas_nobs",
    "tdps_nobs",
    "tdps_derived_nobs",
    "ps_nobs",
    "ps_derived_nobs",
    "psl_nobs",
    "ps_altimeter_nobs",
    "pr_nobs",
    "pr_5min_nobs",
    "pr_1h_nobs",
    "pr_24h_nobs",
    "pr_localmid_nobs",
    "hurs_nobs",
    "sfcwind_nobs",
    "sfcwind_dir_nobs",
    "rsds_nobs",
    "total_nobs",
]


def get_station_list_paths(directory: str) -> dict:
    """Iterate through clean folder and get all station lists

    Parameters
    ----------
    directory : str
        s3 bucket corresponding to option

    Returns
    -------
    networks : dict
        dictionary of all input network paths

    Notes
    -----
    1. Currently doesn't have a method for dealing with more than one normally formatted station file
    """

    # Get list of folder prefixes
    response = s3_cl.list_objects_v2(
        Bucket=BUCKET_NAME, Prefix=directory, Delimiter="/"
    )

    # Set-up df of station list names
    networks = {"Network": [], "NetworkPath": [], "StationFile": []}

    for prefix in response["CommonPrefixes"]:  # For each folder path
        networkpath = prefix["Prefix"][:-1]
        networkname = networkpath.split("/")[-1]
        station_file = (
            s3.Bucket(BUCKET_NAME)
            .objects.filter(Prefix=networkpath + "/" + "stationlist_")
            .all()
        )

        # If 1 stationlist per network
        if len(list(station_file)) == 1:
            for item in station_file:
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(item.key)
                # If more than one file of this format found in folder, just take the most recent
                break

        # If no station lists
        elif len(list(station_file)) == 0:
            # List all files in folder
            files = s3.Bucket(BUCKET_NAME).objects.filter(Prefix=networkpath + "/")
            # More general search for 'station'
            file = [file for file in files if "station" in file.key]

            # Keep all files found here. These files may be different (e.g. ISD ASOS/AWOS vs ASOS/AWOS station lists)
            for item in file:
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(item.key)

        # If more than one identically formatted station list returned, take most recent
        elif len(list(station_file)) > 1:
            # sort station lists by last edit
            file_all = [obj.key for obj in station_file]

            # need to handle for subsetted CWOP stationfiles
            if networkname == "CWOP":
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                # Add path to first in list -- alphabetical
                networks["StationFile"].append(file_all[0])
            else:
                file = [
                    obj.key
                    for obj in sorted(
                        station_file, key=lambda x: x.last_modified, reverse=True
                    )
                ]
                # Add path to most recently changed file
                networks["Network"].append(networkname)
                networks["NetworkPath"].append(networkpath)
                networks["StationFile"].append(file[0])

    return networks


def retrieve_and_concat_stnlists(directory: str, option: str) -> pd.DataFrame:
    """
    Retrieves the stationlists using the paths dictionary, and concats together.

    Parameters
    ----------
    directory : str
        s3 bucket corresponding to option
    option : str
        which stage of development to generate a stationlist for

    Returns
    -------
    df_to_save : pd.DataFrame
    """

    # Set-up all network stationlist
    dffull = pd.DataFrame()

    # Retrieve dictionary of paths of all network stationlists
    networks = get_station_list_paths(directory)
    networks = pd.DataFrame(networks)

    # Remove all Valley Water and Marin County stations from list
    networks = networks[networks["Network"] != "VALLEYWATER"]
    networks = networks[networks["Network"] != "MARINCOUNTY"]
    networks = networks[networks["Network"] != "station_metadata"]
    print(f"{len(networks)} networks retrieved in {directory}...")

    # Check that no network has >1 station file and break code if it does.
    boolean = networks["Network"].is_unique
    total = 0

    if boolean is False:
        dupes = networks["Network"].duplicated()
        doubles = list(networks["Network"][dupes])
        print(
            f"Error: More than one station file found for the following networks: {doubles}. \
              Inspect folder, add duplicated files to remove_list and re-run function."
        )
        exit()

    for index, row in networks.iterrows():
        # Get data
        df = pd.DataFrame()
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=row["StationFile"])  # get file
        body = obj["Body"].read()
        try:
            temp = pd.read_csv(BytesIO(body), encoding="utf8")
        except:
            # there's an encoding error in the CIMIS pull stationlist
            temp = pd.read_excel(BytesIO(body), engine="openpyxl")

        try:
            networkname = row["Network"]
            print(f"Checking for {option} stations: {networkname}")
            total += len(temp)

            # Make column names to lower
            temp.columns = temp.columns.str.lower()

            # Delete index column
            remove = [col for col in temp.columns if "unnamed" in col]
            if remove is not None:
                temp = temp.drop(remove, axis=1)

            # name: station name, station_name, name, dcp location name --> use name as filter
            if option == "pull":
                if any("name" in str for str in temp.columns):
                    colname = [col for col in temp.columns if "name" in col]
                    if len(colname) > 1:
                        # If more than one col returned
                        removelist = set(["countyname"])
                        # Use sets to exclude partial matches (e.g. 'name' in 'countyname')
                        colname = list(set(colname) - removelist)
                        if len(colname) > 1:
                            print(
                                f"Too many options for station name columns. Add manually to removelist: {colname}"
                            )
                            break

                    df["name"] = temp[colname].values.reshape(
                        -1,
                    )

            else:
                if any("era-id" in str for str in temp.columns):
                    colname = [col for col in temp.columns if "era-id" in col]
                    df["era-id"] = temp[colname].values.reshape(-1)

            # latitude: lat or latitude
            if any("lat" in str for str in temp.columns):
                colname = [col for col in temp.columns if "lat" in col]
                if len(colname) > 1:  # If more than one col returned
                    removelist = set([])
                    colname = list(set(colname) - removelist)
                df["latitude"] = temp[colname].values.reshape(-1)
            else:
                df["latitude"] = np.nan

            # longitude: lat or latitude
            if any("lon" in str for str in temp.columns):
                colname = [col for col in temp.columns if "lon" in col]
                if len(colname) > 1:  # If more than one col returned
                    removelist = set([])
                    colname = list(set(colname) - removelist)
                df["longitude"] = temp[colname].values.reshape(-1)
            else:
                df["longitude"] = np.nan

            # elevation: elev or elevation (TO DO: convert to same unit!!)
            if any("elev" in str for str in temp.columns):
                colname = [col for col in temp.columns if "elev" in col]
                if len(colname) > 1:  # If more than one col returned
                    removelist = set(
                        ["elev(m)", "barometer_elev", "anemometer_elev"]
                    )  # remove sensor heights
                    colname = list(set(colname) - removelist)
                    if len(colname) > 1:
                        if "elev_dem" in colname:
                            colname.remove("elev_dem")
                df["elevation"] = temp[colname].values.reshape(-1)
            else:
                df["elevation"] = np.nan

            # start-date: search for start, begin or connect
            if any(y in x for x in temp.columns for y in ["begin", "start", "connect"]):
                colname = [
                    col
                    for col in temp.columns
                    if any(sub in col for sub in ["begin", "start", "connect"])
                ]
                if len(colname) > 1:  # If more than one col returned
                    removelist = set(
                        ["startdate", "begindate"]
                    )  # Add any items to be manually removed here.
                    colname = list(set(colname) - removelist)
                    if len(colname) > 1:
                        # If both start_time (parsed) and begin (not parsed) columns present, remove begin.
                        if "start_time" in colname:
                            if "begin" in colname:
                                colname.remove("begin")
                        if "disconnect" in colname:
                            colname.remove("disconnect")
                df["start-date"] = temp[colname].values.reshape(-1)
            else:  # If no start date provided
                df["start-date"] = np.nan

            # end-date: search for end or disconnect
            if any(y in x for x in temp.columns for y in ["end", "disconnect"]):
                colname = [
                    col
                    for col in temp.columns
                    if any(sub in col for sub in ["end", "disconnect"])
                ]
                if len(colname) > 1:  # If more than one col returned
                    removelist = set(
                        ["enddate"]
                    )  # Add any items to be manually removed here.
                    colname = list(set(colname) - removelist)
                    if len(colname) > 1:
                        # If both start_time (parsed) and begin (not parsed) columns present, remove begin.
                        if "end_time" in colname:
                            if "end" in colname:
                                colname.remove("end")
                df["end-date"] = temp[colname].values.reshape(-1)
            else:  # If no start date provided
                df["end-date"] = np.nan

            # Add stage and date checked columns, if they exist
            if any("pulled" in str for str in temp.columns):
                df["pulled"] = temp["pulled"].values.reshape(-1)
            else:
                df["pulled"] = np.nan

            if any("time_checked" in str for str in temp.columns):
                df["time_checked"] = temp["time_checked"].values.reshape(-1)
            else:
                df["time_checked"] = np.nan

            if any("cleaned" in str for str in temp.columns):
                df["cleaned"] = temp["cleaned"].values.reshape(-1)
            else:
                df["cleaned"] = np.nan

            if any("time_cleaned" in str for str in temp.columns):
                df["time_cleaned"] = temp["time_cleaned"].values.reshape(-1)
            else:
                df["time_cleaned"] = np.nan

            if any("qaqc" in str for str in temp.columns):
                df["qaqc"] = temp["qaqc"].values.reshape(-1)
            else:
                df["qaqc"] = np.nan

            if any("time_qaqc" in str for str in temp.columns):
                df["time_qaqc"] = temp["time_qaqc"].values.reshape(-1)
            else:
                df["time_qaqc"] = np.nan

            if any("merged" in str for str in temp.columns):
                df["merged"] = temp["merged"].values.reshape(-1)
            else:
                df["merged"] = np.nan

            if any("time_merge" in str for str in temp.columns):
                df["time_merge"] = temp["time_merge"].values.reshape(-1)
            else:
                df["time_merge"] = np.nan

            # network name column
            df["network"] = networkname

            # Cleaned variable coverage
            for clean_var in CLEANED_VARS:
                if any(clean_var in str for str in temp.columns):
                    df[clean_var] = temp[clean_var].values.reshape(-1)
                else:
                    df[clean_var] = np.nan

            dffull = pd.concat([dffull, df], sort=False)

        except Exception as e:
            print(e)

    # Organize full dataframe
    # If end date is "active", make this be today's date
    today = datetime.datetime.now()
    dffull["end-date"] = dffull["end-date"].replace("Active", today)

    # # Format dates in datetime format
    dffull["start-date"] = pd.to_datetime(
        dffull["start-date"], utc=True, format="mixed"
    )
    dffull["end-date"] = pd.to_datetime(dffull["end-date"], utc=True, format="mixed")

    # # Remove any duplicates (of network and ID)
    if option == "pull":
        dffull.drop_duplicates(
            subset=["name", "latitude", "longitude", "network"], inplace=True
        )
    else:
        dffull.drop_duplicates(
            subset=["era-id", "latitude", "longitude", "network"], inplace=True
        )

    # Resort by network
    dffull.sort_values(by=["network"], inplace=True)

    # Reset index
    dffull = dffull.reset_index(drop=True)
    print("Stationlists concatenated together...")

    return dffull


def export_stationlist(df_to_save: pd.DataFrame, directory: str, option: str):
    """
    Helper function to export new stationlist.

    Parameters
    ----------
    df_to_save : pd.DataFrame
        concatenated stationlist to export
    directory : str
        s3 bucket corresponding to option
    option : str
        which stage of development to generate a stationlist for

    Returns
    -------
    None
    """

    # set-up correct csv to export
    df_to_save = stationlist_cols(df_to_save, option)

    # Export to s3
    s3_path = f"s3://{BUCKET_NAME}/{directory}all_network_stationlist_{option}.csv"
    df_to_save.to_csv(s3_path, na_rep="NaN")
    print(f"all_network_stationlist_{option}.csv generated and saved to AWS.")

    return None


def stationlist_cols(df_to_save: pd.DataFrame, option: str) -> pd.DataFrame:
    """
    Helper function to that sets which columns to export per stage.

    Parameters
    ----------
    df_to_save : pd.DataFrame
        concatenated stationlist to export
    option : str
        which stage of development to generate a stationlist for

    Returns
    -------
    df_to_export : pd.DataFrame
        concatenated stationlist with appropriate columns
    """

    pull_cols = [
        "name",
        "latitude",
        "longitude",
        "elevation",
        "start-date",
        "end-date",
        "pulled",
        "time_checked",
        "network",
    ]
    clean_cols = pull_cols + ["cleaned", "time_cleaned"] + CLEANED_VARS
    clean_cols = [item.replace("name", "era-id") for item in clean_cols]
    qaqc_cols = clean_cols + ["qaqc", "time_qaqc"]
    merge_cols = qaqc_cols + ["merged", "time_merge"]

    if option == "pull":
        df_to_export = df_to_save[pull_cols]

    if option == "clean":
        df_to_export = df_to_save[clean_cols]

    if option == "qaqc":
        df_to_export = df_to_save[qaqc_cols]

    if option == "merge":
        df_to_export = df_to_save[merge_cols]

    return df_to_export


def generate_stationlist(option: str):
    """
    Generates the all network stationlist, from the input per network stationlists,
    based on user selected stage of development (option).

    Parameters
    ----------
    option : str
        which stage of development to generate a stationlist for

    Returns
    -------
    None
    """
    print(f"Generating {option} all network station list...")

    if option == "pull":
        directory = PULL_DIR
        option_ed = "pulled"
    elif option == "clean":
        directory = CLEAN_DIR
        option_ed = "cleaned"
    elif option == "qaqc":
        directory = QAQC_DIR
        option_ed = option
    elif option == "merge":
        directory = MERGE_DIR
        option_ed = "merged"

    df_to_save = retrieve_and_concat_stnlists(directory, option)
    export_stationlist(df_to_save, directory, option)

    # Print some useful statistics
    if option == "pull":
        # CIMIS and CW3E have complicated data formats and read as NaN in station list
        # 15,730 + 261 CIMIS stations + 13 CW3E stations
        counts = 16004
    else:
        counts = df_to_save[option_ed].value_counts()["Y"]
    print(
        f"Number of {option} stations: {counts} out of {len(df_to_save)} possible stations."
    )
    print(f"Successful {option} station rate: {(counts / len(df_to_save)) * 100}")

    return None


# =======================================================================
if __name__ == "__main__":
    # generate station chart
    generate_stationlist(option="merge")

    # Options: pull, clean, qaqc, merge
