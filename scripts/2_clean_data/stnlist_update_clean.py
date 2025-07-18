"""
stnlist_update_clean.py

This script iterates through a specified network and checks to see what stations have been successfully cleaned,
updating the station list in the 1_raw_wx folder to reflect station availability. Error.csvs in the cleaned bucket are also parsed,
with relevant errors added to the corresponding stations if station files are not cleaned, or if the errors occur during or after the cleaning process.

Functions
---------
- get_station_list: Retrieves network station list.
- get_cleaned_stations: Adds a new column to station list with bool "Y" or "N" status upon cleaning.
- parse_error_csv: Retrieves error csv file for a network.
- clean_qa: Updates station list and saves to AWS, adding cleaned status, time of clean, and relevant errors.
- cwop_stnlist_merge: Merges CWOP stationlists together, overwrites full stationlist file for CWOP.

Intended Use
------------
Run this script after cleaning has been completed for a network to update the network stationlist for cleaning success rate tracking.

Notes
-----
Because errors.csv are parsed, very old errors.csv may want to be removed manually from AWS or thresholded below.
"""

import boto3
import pandas as pd
from io import BytesIO, StringIO
import numpy as np
import xarray as xr
import s3fs
from datetime import datetime
import re

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")
BUCKET_NAME = "wecc-historical-wx"
RAW_WX = "1_raw_wx/"
CLEAN_WX = "2_clean_wx/"


def get_station_list(network: str) -> pd.DataFrame:
    """
    Retrieves network station list.

    Parameters
    ----------
    network : str
        network name

    Returns
    -------
    station_list : pd.DataFrame
        station list for a specified network
    """

    network_prefix = RAW_WX + network + "/"

    # If station list is CIMIS, extension is .xlsx
    if network == "CIMIS":
        station_list = f"stationlist_{network}.xlsx"
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=network_prefix + station_list)
        station_list = pd.read_excel(BytesIO(obj["Body"].read()))
    elif network == "ASOSAWOS":
        # Use merged station list
        network_prefix = CLEAN_WX + network + "/"
        station_list = f"stationlist_{network}_merge.csv"
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=network_prefix + station_list)
        station_list = pd.read_csv(obj["Body"])
    else:
        station_list = f"stationlist_{network}.csv"
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=network_prefix + station_list)
        station_list = pd.read_csv(obj["Body"])
    return station_list


def get_cleaned_stations(network: str) -> pd.DataFrame:
    """
    Adds a new column to station list with bool "Y" or "N" status upon cleaning.

    Parameters
    ----------
    network : str
        network name

    Returns
    -------
    pd.DataFrame
        station list with new flag of whether data was successfully cleaned
    """

    df = {"ID": [], "Time_Cleaned": []}
    network_prefix = CLEAN_WX + network + "/"
    for item in s3.Bucket(BUCKET_NAME).objects.filter(
        Prefix=network_prefix + network + "_"
    ):
        if network == "CW3E":
            # cleaned CW3E data is stored as one file per year per station, handle separately
            stn = item.key.split("_")[3]
            clean_id = network + "_" + stn
            time_mod = item.last_modified
        else:
            # all other networks, one file per station
            clean_id = item.key.split("/")[-1].replace(".nc", "")
            time_mod = item.last_modified
        df["ID"].append(clean_id)
        df["Time_Cleaned"].append(time_mod)

    return pd.DataFrame(df)


def parse_error_csv(network: str) -> pd.DataFrame:
    """
    Retrieves error csv file for a network.

    Parameters
    ----------
    network : str
        network name

    Returns
    -------
    errordf : pd.DataFrame
        errors data per station
    """

    errordf = []
    errors_prefix = f"{CLEAN_WX}{network}/errors"

    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=errors_prefix):
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item.key)
        errors = pd.read_csv(obj["Body"])
        if errors.empty:
            # If file empty
            continue
        else:
            errors = errors[["File", "Time", "Error"]]
            errordf.append(errors)

    if not errordf:
        # If no errors in cleaning
        return pd.DataFrame()

    else:
        errordf = pd.concat(errordf)
        errordf = errordf.drop_duplicates(subset=["File", "Error"])
        # Drop any whole network errors
        errordf = errordf[errordf.File != "Whole network"]

        return errordf


def clean_qa(network: str, clean_var_add: bool = False, cwop_letter: str | None = None):
    """
    Updates station list and saves to AWS, adding cleaned status, time of clean, and relevant errors.

    Parameters
    ----------
    network : str
        network name
    clean_var_add : bool, optional
        adds variable observation count columns to staiton list, default is False
    cwop_letter : str, optional
        CWOP subsetting letter, default is None

    Returns
    -------
    None

    Notes
    -----
    1. clean_var_add will open every cleaned station file, check which variables are present, and flag in station list.
    It is a highly time intensive process, so recommendation is to only run clean_var_add=True after a full clean, or a partial clean update (new data added).
    """

    if "otherisd" in network:
        # Fixing capitalization issues
        network = "OtherISD"
    else:
        network = network.upper()

    # Call functions
    stations = get_station_list(network)
    cleanids = get_cleaned_stations(network)
    errors = parse_error_csv(network)

    if cleanids.empty:
        print(
            "No cleaned files for this network. Please run the relevant cleaning script and try again."
        )
        exit()

    # Retrieve clean station list
    # Drop double index column
    stations = stations.loc[:, ~stations.columns.str.match("Unnamed")]

    # Add standardized ID column - ERA-ID
    if "ASOS" in network or "OtherISD" in network:
        stations["ERA-ID"] = network + "_" + stations["ISD-ID"].str.replace("-", "")
    elif "CIMIS" in network:
        # Drop one last row that is not standardized
        stations = stations.dropna(subset=["Station Number"])
        stations["ERA-ID"] = (
            network + "_" + stations["Station Number"].astype("int").astype("str")
        )
    elif "CW3E" in network:
        stations["ERA-ID"] = network + "_" + stations["STID"].str.replace("C3", "")
    elif network in [
        "CAHYDRO",
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
        "RAWS",
        "SGXWFO",
        "SHASAVAL",
        "VCAPCD",
    ]:
        # MADIS networks
        stations["ERA-ID"] = network + "_" + stations["STID"]
    elif network in ["MARITIME", "NDBC"]:
        stations["ERA-ID"] = network + "_" + stations["STATION_ID"]
    elif network in ["SCAN", "SNOTEL"]:
        stations["ERA-ID"] = (
            network + "_" + stations["stationTriplet"].str.split(":").str[0]
        )

    # Make ERA-ID first column
    eraid = stations.pop("ERA-ID")
    # Standardize names (one outlier station in CWOP)
    eraid = eraid.str.upper()
    stations.insert(0, "ERA-ID", eraid)

    # Join cleaned columns to column list
    stations = stations.merge(cleanids, left_on="ERA-ID", right_on="ID", how="outer")
    if "ID_y" in stations.columns:
        # Make binary cleaned column
        stations["Cleaned"] = np.where(stations.ID_y.isna(), "N", "Y")
        # Drop ID column
        stations = stations.drop(["ID_x", "ID_y"], axis=1)
    else:
        #  Make binary cleaned column
        stations["Cleaned"] = np.where(stations.ID.isna(), "N", "Y")
        # Drop ID column
        stations = stations.drop("ID", axis=1)

    # Move Time_Cleaned to last
    s = stations.pop("Time_Cleaned")
    stations = pd.concat([stations, s], axis=1)

    # Add errors to column by station - only add error if error occurred at or after file clean, if file cleaned.
    stations["Errors"] = np.nan

    # Remove any NAs from ERA-ID
    stations = stations.loc[stations["ERA-ID"].notnull()]

    # If station not in stations list, add it manually
    cleanids_list = [x for x in cleanids.ID if any(y in x for y in stations["ERA-ID"])]
    # Get any ids that aren't in the stations list already
    cleanids_nolist = [x for x in cleanids.ID.tolist() if x not in cleanids_list]
    cleanids_nolist = cleanids.loc[cleanids.ID.isin(cleanids_nolist)]
    if not cleanids_nolist.empty:
        not_in_list = pd.DataFrame(
            {
                "ERA-ID": cleanids_nolist.ID,
                "Cleaned": "Y",
                "Time_Cleaned": cleanids_nolist.Time_Cleaned,
            }
        )
        stations = pd.concat([stations, not_in_list])

    # Get list of station IDs
    ids = [id.split("_")[-1] for id in stations["ERA-ID"].tolist()]

    if errors.empty:
        # If no errors, stop here.
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

        for index, row in stations.iterrows():
            # For each station
            error_sta = errors.loc[errors.ID == row["ERA-ID"]]
            if error_sta.empty:
                # if no errors for station
                continue
            else:
                if not pd.isnull(row["Time_Cleaned"]):
                    # Only keep errors from cleaning at or after time of clean
                    error_sta = error_sta.loc[
                        (error_sta.Time >= row["Time_Cleaned"])
                        | (error_sta.Time.isna()),
                        :,
                    ]

                if len(error_sta) == 1:
                    stations.loc[index, "Errors"] = error_sta["Error"].values[0]
                elif len(error_sta) > 1:
                    values = [
                        f"{x.File}: {x.Error}" for index, x in error_sta.iterrows()
                    ]
                    value = " ".join(values)
                    stations.loc[index, "Errors"] = value

    # Print summary
    clean_y = stations["Cleaned"].value_counts()["Y"]
    clean_n = stations["Cleaned"].value_counts()["N"]

    if network != "CW3E":
        if "Y" in stations["Cleaned"].values:
            # order is important here, if no "N" is present in a cleaned network, it will bark without this
            if "N" not in stations["Cleaned"].values:
                print(
                    f"Station list updated for cleaned {network} stations. All stations cleaned: {clean_y} stations cleaned."
                )
            else:
                print(
                    f"Station list updated for cleaned {network} stations. {clean_y} stations cleaned, {clean_n} stations not cleaned."
                )
        else:
            print(
                f"Station list updated for cleaned {network} stations. No stations cleaned successfully. {clean_n} stations not yet cleaned."
            )
    else:
        # network is CW3E
        if "Y" in stations["Cleaned"].values:
            if "N" not in stations["Cleaned"].values:
                print(
                    f"Station list updated for cleaned {network} stations. All stations cleaned: {clean_y} station-years cleaned."
                )
            else:
                print(
                    f"Station list updated for cleaned {network} stations. {clean_y} station-years cleaned, {clean_n} station-years not cleaned."
                )
        else:
            print(
                f"Station list updated for cleaned {network} stations. No stations cleaned successfully. {clean_n} stations not yet cleaned."
            )

    # clean_var_add, with subsetting for CWOP
    if clean_var_add == True:
        print(
            "Processing all cleaned files to assess variable coverage -- this may take awhile based on size of network!"
        )

        # Set up error handling for identifing cleaned files that can't open
        # Appears that some datetime/index failed to format -- need to reclean
        errors = {"File": [], "Time": [], "Error": []}
        # Set end time to be current time at beginning of download: for error handling csv.
        end_api = datetime.now().strftime("%Y%m%d%H%M")
        # For attributes of netCDF file
        timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

        # add in default columns of "N" to cleaned station list for all core and associated variables
        # also adds column that counts number of valid/non-nan observations
        core_vars = [
            "tas",
            "tdps",
            "tdps_derived",
            "ps",
            "psl",
            "ps_altimeter",
            "ps_derived",
            "pr",
            "pr_5min",
            "pr_1h",
            "pr_24h",
            "pr_localmid",
            "hurs",
            "sfcWind",
            "sfcWind_dir",
            "rsds",
        ]
        for var in core_vars:
            stations[str(var)] = "N"
            stations[str(var + "_nobs")] = 0  # default of 0 to start

        # add column for total length of each record, valid (non-nan) and nans
        stations["total_nobs"] = 0  # default of 0 to start

        # open cleaned datafile
        network_prefix = CLEAN_WX + network + "/"
        files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=network_prefix):
            file = str(item.key)
            files += [file]

        # Get list of cleaned file names
        files = list(filter(lambda f: f.endswith(".nc"), files))

        # get list of all station filenames successfully cleaned, and filter by subsetting (CWOP)
        if network != "CWOP":
            files = list(filter(lambda f: f.endswith(".nc"), files))
        elif network == "CWOP" and cwop_letter == None:
            # in case all CWOP is run at once
            print(
                "Warning: Setting cwop_letter = None is for an entire network update for CWOP, estimated 3+ days to complete."
            )
            files = list(filter(lambda f: f.endswith(".nc"), files))
        elif network == "CWOP" and cwop_letter != None:
            # subsetting in place for CWOP
            # Procedure for grouping of data in CWOP to split up 7k+ stations by first letter
            not_ABCDEFG = (
                "A",
                "B",
                "C",
                "D",
                "E",
                "F",
                "G",
            )
            # catch-all single letter stations (K, L, M, P, S, T, U, W at present)
            if "other" in cwop_letter and len(cwop_letter) == 5:
                # cwop_letter = "other"
                ids = [id for id in files if not id[-8].startswith(not_ABCDEFG)]

            elif "other" in cwop_letter and len(cwop_letter) != 5:
                # additional letters + other category called, ex: cwop_letter = "ABC + other"
                letter_to_clean = cwop_letter.replace(" ", "")
                letter_to_clean = letter_to_clean.replace("other", "")
                # so it doesn't clean "o t h e r"
                letter_to_clean = letter_to_clean.replace("+", "")
                letter_ids = tuple(letter_to_clean)
                other_ids = [id for id in files if not id[-8].startswith(not_ABCDEFG)]
                letter_ids = [id for id in files if id[-8].startswith(letter_ids)]
                ids = other_ids + letter_ids

            if len(cwop_letter) == 1:
                # single letter cleaning, ex: cwop_letter = "A"
                ids = [id for id in files if id[-8].startswith(str(cwop_letter))]

            if "other" not in cwop_letter and len(cwop_letter) != 1:
                # more than one letter provided, but not other category, ex: cwop_letter = "ACD"
                letter_ids = tuple(cwop_letter)
                ids = [id for id in files if id[-8].startswith(letter_ids)]

            print(
                f"CWOP batch variable coverage update for '{cwop_letter}' stations: batch-size of {len(ids)} stations"
            )
            # resetting subset to main file list to be consistent for all networks and subset options
            files = ids

            # generates subsetted list of just the station ids, needed to subset for CWOP
            stnids_to_check = []
            for item in ids:
                stnids = str(item).split("/")
                stnids = stnids[-1].replace(".nc", "")  # drops .nc
                stnids_to_check += [stnids]

        for file in files:
            if file not in files:
                # dont run qa/qc on a station that isn't cleaned
                continue
            else:
                try:
                    print(file)
                    fs = s3fs.S3FileSystem()
                    aws_url = f"s3://wecc-historical-wx/{file}"

                    if network == "CWOP" and cwop_letter != None:
                        # produces stationlist update of only stations within that cwop_letter so it doesnt overwrite
                        # uses subsetted station ids list
                        stations = stations.loc[
                            stations["ERA-ID"].isin(stnids_to_check)
                        ]

                    with fs.open(aws_url) as fileObj:
                        # setting engine=None (default) uses what is best for system, previously engine='h5netcdf'
                        ds = xr.open_dataset(fileObj)
                        # fill total num of obs with length of time var
                        stations.loc[
                            stations["ERA-ID"] == ds.station.values[0], "total_nobs"
                        ] = ds.time.shape[0]

                        # mark each variable as present if in dataset, and count number of valid/non-nan values
                        for var in ds.variables:
                            if var in core_vars:
                                stations.loc[
                                    stations["ERA-ID"] == ds.station.values[0], str(var)
                                ] = "Y"
                                stations.loc[
                                    stations["ERA-ID"] == ds.station.values[0],
                                    str(var + "_nobs"),
                                ] = ds[str(var)].count()

                        # close dataset
                        ds.close()

                except Exception as e:
                    print(f"{file} not opening")
                    errors["File"].append(file)
                    errors["Time"].append(end_api)
                    errors["Error"].append(f"clean_var_add error in opening file: {e}")
                    continue

        # reset index
        stations = stations.reset_index(drop=True)

        # Save errors file to cleaned bucket
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=f"{CLEAN_WX}{network}/add_clean_var_errors_{network}_{end_api}.csv",
        )

    # Save station file to cleaned bucket
    new_buffer = StringIO()
    stations.to_csv(new_buffer, index=False)
    content = new_buffer.getvalue()

    # set different files for CWOP if subsetting
    if clean_var_add == False:
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=f"{CLEAN_WX}{network}/stationlist_{network}_cleaned.csv",
        )

    else:
        # clean_var_add == True
        if network == "CWOP" and cwop_letter != None:
            s3_cl.put_object(
                Bucket=BUCKET_NAME,
                Body=content,
                Key=f"{CLEAN_WX}{network}/stationlist_{network}_cleaned_{cwop_letter}.csv",
            )
        else:
            s3_cl.put_object(
                Bucket=BUCKET_NAME,
                Body=content,
                Key=f"{CLEAN_WX}{network}/stationlist_{network}_cleaned.csv",
            )

    return None


def cwop_stnlist_merge(network: str):
    """
    Merges CWOP stationlists together, overwrites full stationlist file for CWOP.

    Parameters
    ----------
    network : str
        network name

    Returns
    -------
    None
    """

    # Check network is CWOP
    if network != "CWOP":
        print(
            f"Incorrect network ({network}) provided to cwop_stnlist_merge! Please pass 'CWOP' as network."
        )
        exit()

    # should be (A, B, C, D, E, F, G, other) -- assumes that no groupings were applied (i.e., "AB + other")
    station_files = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix="2_clean_wx/CWOP/"):
        file = str(item.key)
        station_files += [file]

    # grabbing only station files
    station_files = [file for file in station_files if "station" in file]
    # grabbing only subset lists
    station_files = [file for file in station_files if len(file) > 45]

    all_df = pd.DataFrame()
    for item in station_files:
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=item)
        body = obj["Body"].read()
        df = pd.read_csv(BytesIO(body), encoding="utf8")
        all_df = pd.concat([all_df, df], ignore_index=True)

    # sorts alphabetically and resets the index
    all_df = all_df.sort_values("ERA-ID", ascending=True).reset_index(drop=True)

    # save to s3 bucket
    csv_buffer = StringIO()
    all_df.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=f"{CLEAN_WX}{network}/stationlist_{network}_cleaned.csv",
    )

    return None


if __name__ == "__main__":
    clean_qa("CWOP", clean_var_add=False, cwop_letter=None)
    # cwop_stnlist_merge("CWOP")  # Use once all CWOP stationlist(s) are updated with variable coverage

    # List of all stations for ease of use here:
    # ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, HOLFUY, HPWREN, LOXWFO
    # MAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, otherisd, RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME
    # NDBC, SCAN, SNOTEL

    # Note: OtherISD only runs as "otherisd"
    # Note: Make sure there is no space in the name CAHYDRO ("CA HYDRO" will not run)
