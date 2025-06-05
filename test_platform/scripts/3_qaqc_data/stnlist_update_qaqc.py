"""
stnlist_update_qaqc.py

This script iterates through a specified network and checks to see what stations have been successfully passed quality control,
updating the station list in the 1_raw_wx folder to reflect station availability. Error.csvs in the cleaned bucket are also parsed,
with relevant errors added to the corresponding stations if station files do not pass QA/QC, or if the errors occur during or after the QA/QC process.

Note that because errors.csv are parsed, very old errors.csv may want to be removed manually from AWS or thresholded below
(removing those produced during code testing). 

Functions
---------
- get_station_list: Retrieves specific network stationlist from clean bucket
- get_zarr_last_mod: Identifies the last modified date from a zarr 
- get_qaqc_stations: Retrieves list of all stations that pass QAQC
- parse_error_csv: Retrieves all processing error files for a network
- qaqc_qa: Processing function that updates the stationlist with QA/QC status

Intended Use
------------
Run this script after QAQC has been completed for a network (via pcluster run) to update the network stationlist for QA/QC success rate tracking.
"""

import pandas as pd
from io import BytesIO, StringIO
import numpy as np
import boto3
import s3fs

# Set environment variables
BUCKET_NAME = "wecc-historical-wx"
RAW_WX = "1_raw_wx/"
CLEAN_WX = "2_clean_wx/"
QAQC_WX = "3_qaqc_wx/"
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")


def get_station_list(network: str) -> pd.DataFrame:
    """
    Given a network name, return a pandas dataframe containing the network's station list from the clean bucket.
    Intentionally grabbing the cleaned version of the stationlist to retain information about whether a station was cleaned going into QA/QC process.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    station_list : pd.DataFrame
        station list of all stations within a network

    Notes
    -----
    Can be updated to read directly from AWS
    """
    network_prefix = CLEAN_WX + network + "/"
    station_list = f"stationlist_{network}_cleaned.csv"
    obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=network_prefix + station_list)
    station_list = pd.read_csv(obj["Body"])
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

    path_no_ext = fn.split(".")[0]  # Grab only path name without extension
    path_no_bucket = path_no_ext.split(BUCKET_NAME)[-1][
        1:
    ]  # Grab only part without bucket name

    # idenitfy last_modified date from metadata date within .zarr
    mod_list = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=path_no_bucket):
        mod_list.append(str(item.last_modified))

    # return most recent datetime value
    last_mod = max(mod_list)
    return last_mod


def get_qaqc_stations(network: str) -> pd.DataFrame:
    """
    Retrieves a list of all stations in a given network that have passed the QA/QC process.

    This function scans the `3_qaqc_wx` folder in the AWS S3 bucket for the specified network,
    identifying which stations have successfully completed QA/QC based on the presence of `.zarr` files.
    It returns a DataFrame containing the station ID, a placeholder for the timestamp when QA/QC occurred,
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
        - 'Time_QAQC': (currently empty; placeholder for QA/QC timestamp)
        - 'QAQC': 'Y' if a .zarr file exists (passed QA/QC), 'N' otherwise
    """
    df = {"ID": [], "Time_QAQC": [], "QAQC": []}  # Initialize results dictionary

    # Construct the S3 path prefix for the network inside the QAQC folder
    parent_s3_path = f"{BUCKET_NAME}/{QAQC_WX}{network}"

    # Use s3fs to list all items under this path
    s3_fs = s3fs.S3FileSystem(anon=False)
    all_paths = s3_fs.ls(parent_s3_path)

    # Filter to only .zarr folders
    zarr_folders = [f"{path}" for path in all_paths if path.endswith(".zarr")]

    for item in zarr_folders:
        # Extract the station ID from the folder name, which is usually the last part of the path
        station_id = item.split("/")[-1].split(".")[-2].replace(".zarr", "")
        df["ID"].append(station_id)
        df["QAQC"].append("Y")  # QA/QC passed

        # Retrieving "last_modified" timestamp on a .zarr requires going into group
        time_mod = get_zarr_last_mod(item)
        df["Time_QAQC"].append(time_mod)

    return pd.DataFrame(df)


def parse_error_csv(network: str) -> pd.DataFrame:
    """
    Given a network name, return a pandas dataframe containing all errors reported for the network in the QAQC stage.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    errordf : pd.Dataframe
        dataframe of all errors produced during QAQC
    """
    errordf = []  # List to store all non-empty error DataFrames found

    # Define the path prefix in the S3 bucket for error CSVs
    errors_prefix = f"{QAQC_WX}{network}/qaqc_errs/errors"

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
        errordf = errordf[
            errordf.File != "Whole network"
        ]  # Drop any whole network errors
        return errordf


def qaqc_qa(network: str):
    """
    Update station list and save to AWS, adding qa/qc status, time of qa/qc pass and any relevant errors.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    None
    """
    if "otherisd" in network:  # Fixing capitalization issues
        network = "OtherISD"
    else:
        network = network.upper()

        # Call functions
        stations = get_station_list(network)  # grabs stationlist_cleaned
        qaqc_ids = get_qaqc_stations(network)  # grabs stations that pass qaqc
        errors = parse_error_csv(network)

        if qaqc_ids.empty:
            print(
                "No QA/QC'd files for this network. Please run the relevant qaqc script and try again."
            )
            exit()

        # Join cleaned columns to column list
        stations = stations.merge(
            qaqc_ids, left_on="ERA-ID", right_on="ID", how="outer"
        )
        if "ID_y" in stations.columns:
            stations["QAQC"] = np.where(
                stations.ID_y.isna(), "N", "Y"
            )  # Make binary qa/qc column
            # Drop ID column
            stations = stations.drop(["ID_x", "ID_y"], axis=1)
        else:
            stations["QAQC"] = np.where(
                stations.ID.isna(), "N", "Y"
            )  # Make binary qa/qc column
            # Drop ID column
            stations = stations.drop("ID", axis=1)

        # Move Time_Cleaned to last
        s = stations.pop("Time_QAQC")
        stations = pd.concat([stations, s], axis=1)

        # Add errors to column by station - only add error if error occurred at or after file clean, if file cleaned.
        stations["Errors_QAQC"] = np.nan

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
            errors["Time"] = pd.to_datetime(
                errors["Time"], format="%Y%m%d%H%M", utc=True
            )

            for index, row in errors.iterrows():
                id = []
                id = [x for x in ids if x in row["File"]]
                if id:
                    errors.loc[index, "ID"] = network + "_" + id[-1]

            for index, row in stations.iterrows():  # For each station
                error_sta = errors.loc[errors.ID == row["ERA-ID"]]
                if error_sta.empty:  # if no errors for station
                    continue
                else:
                    if not pd.isnull(row["Time_QAQC"]):  # If file cleaned
                        error_sta = error_sta.loc[
                            (error_sta.Time >= row["Time_QAQC"])
                            | (error_sta.Time.isna()),
                            :,
                        ]  # Only keep errors from qaqc at or after time of qaqc

                    if len(error_sta) == 1:
                        stations.loc[index, "Errors_QAQC"] = error_sta["Error"].values[
                            0
                        ]
                    elif len(error_sta) > 1:
                        values = [
                            f"{x.File}: {x.Error}" for index, x in error_sta.iterrows()
                        ]
                        value = " ".join(values)
                        stations.loc[index, "Errors_QAQC"] = value

        # Print summary
        if "Y" in stations["QAQC"].values:
            if (
                "N" not in stations["QAQC"].values
            ):  # order is important here, if no "N" is present in a qa/qc'd network, it will bark without this
                print(
                    "Station list updated for {} stations that pass QA/QC. All stations pass: {} stations.".format(
                        network, stations["QAQC"].value_counts()["Y"]
                    )
                )
            else:
                print(
                    "Station list updated for {} stations that pass QA/QC. {} stations passed QA/QC, {} stations were not QA/QC.".format(
                        network,
                        stations["QAQC"].value_counts()["Y"],
                        stations["QAQC"].value_counts()["N"],
                    )
                )
        else:
            print(
                "Station list updated for {} stations. No stations successfully pass QA/QC. {} stations do not yet pass QA/QC.".format(
                    network, stations["QAQC"].value_counts()["N"]
                )
            )

    # Save station file to cleaned bucket
    new_buffer = StringIO()
    stations.to_csv(new_buffer, index=False)
    content = new_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=QAQC_WX + network + "/stationlist_{}_qaqc.csv".format(network),
    )


if __name__ == "__main__":
    qaqc_qa("ASOSAWOS")

# List of all stations for ease of use here:
# ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, HOLFUY, HPWREN, LOXWFO
# MAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME
# NDBC, SCAN, SNOTEL, VALLEYWATER

# Note: OtherISD only runs as "otherisd"
# Note: Make sure there is no space in the name CAHYDRO ("CA HYDRO" will not run)
