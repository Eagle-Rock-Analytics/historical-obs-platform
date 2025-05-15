"""
This script iterates through a specified network and checks to see what stations have been successfully passed quality control,
updating the station list in the 1_raw_wx folder to reflect station availability. Error.csvs in the cleaned bucket are also parsed,
with relevant errors added to the corresponding stations if station files do not pass QA/QC, or if the errors occur during or after the QA/QC process.

Note that because errors.csv are parsed, very old errors.csv may want to be removed manually from AWS or thresholded below
(removing those produced during code testing)
"""

import pandas as pd
from io import BytesIO, StringIO
import numpy as np
import boto3

# Set environment variables
bucket_name = "wecc-historical-wx"
raw_wx = "1_raw_wx/"
clean_wx = "2_clean_wx/"
qaqc_wx = "3_qaqc_wx/"
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")


# ----------------------------------------------------------------------
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
    network_prefix = clean_wx + network + "/"
    station_list = f"stationlist_{network}_cleaned.csv"
    obj = s3_cl.get_object(Bucket=bucket_name, Key=network_prefix + station_list)
    station_list = pd.read_csv(obj["Body"])
    return station_list


# ----------------------------------------------------------------------
def get_qaqc_stations(network: str) -> pd.DataFrame:
    """
    Given a network name, return a pandas dataframe of all stations that pass QA/QC in
    the 3_qaqc_wx AWS bucket, with the date that the file was last modified, and a column that states whether or not it has a .zarr file.

    Parameters
    ----------
    network : str
        name of network

    Returns
    -------
    pd.DataFrame
        pandas dataframe of all stations that pass QA/QC in the 3_qaqc_wx AWS bucket
    """
    df = {"ID": [], "Time_QAQC": [], "QAQC": []}
    network_prefix = f"{qaqc_wx}{network}/"

    zarr_folders = [
        f"{network_prefix}"
        for path in network_prefix
        if path.endswith(".zarr")
        and (
            any(prefix in path.split("/")[-1] for prefix in network_prefix)
            if network_prefix
            else print(len(zarr_folders))
        )
    ]

    print(len(zarr_folders))
    for item in zarr_folders:
        key = item.key
        # Skip any files that are not .nc or .zarr
        # if network == "CW3E":
        #     print("forthcoming")
        #     continue

        # file_path = item.key.split("/")[-2]
        if key.endswith(".nc"):
            station_id = key.split("/")[-1].replace(".nc", "")
            df["ID"].append(station_id)
            df["Time_QAQC"].append(item.last_modified)
            df["QAQC"].append("N")
        elif key.endswith(".zarr"):
            station_id = key.split("/")[-2]  # folder name before trailing slash
            df["ID"].append(station_id)
            df["Time_QAQC"].append(item.last_modified)
            df["QAQC"].append("Y")
        else:
            continue  # skip non-nc/zarr keys

    return pd.DataFrame(df)


# ----------------------------------------------------------------------
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
    errordf = []
    errors_prefix = f"{qaqc_wx}{network}/qaqc_errs/errors"
    for item in s3.Bucket(bucket_name).objects.filter(Prefix=errors_prefix):
        obj = s3_cl.get_object(Bucket=bucket_name, Key=item.key)
        errors = pd.read_csv(obj["Body"])
        if errors.empty:  # If file empty
            continue
        else:
            errors = errors[["File", "Time", "Error"]]
            errordf.append(errors)
    if not errordf:  # If no errors in cleaning
        return pd.DataFrame()
    else:
        errordf = pd.concat(errordf)
        errordf = errordf.drop_duplicates(subset=["File", "Error"])
        errordf = errordf[
            errordf.File != "Whole network"
        ]  # Drop any whole network errors
        return errordf


# ----------------------------------------------------------------------
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
        Bucket=bucket_name,
        Body=content,
        Key=qaqc_wx + network + "/stationlist_{}_qaqc.csv".format(network),
    )


# ---------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    qaqc_qa("NDBC")

    # List of all stations for ease of use here:
    # ASOSAWOS, CAHYDRO, CIMIS, CW3E, CDEC, CNRFC, CRN, CWOP, HADS, HNXWFO, HOLFUY, HPWREN, LOXWFO
    # MAP, MTRWFO, NCAWOS, NOS-NWLON, NOS-PORTS, RAWS, SGXWFO, SHASAVAL, VCAPCD, MARITIME
    # NDBC, SCAN, SNOTEL, VALLEYWATER

    # Note: OtherISD only runs as "otherisd"
    # Note: Make sure there is no space in the name CAHYDRO ("CA HYDRO" will not run)


# ------------------------------------------------------------------------------------------------------------------------------------------------
def _station_has_zarr(network: str, station_id: str) -> str:
    """
    Check if a station has a corresponding .zarr file in the QAQC bucket.

    Parameters
    ----------
    network : str
        Name of network
    station_id : str
        Station ID to check

    Returns
    -------
    str
        "Y" if .zarr file exists, "N" otherwise
    """
    prefix = f"{qaqc_wx}{network}/{station_id}"
    response = s3_cl.list_objects_v2(
        Bucket=bucket_name,
        Prefix=prefix,
        MaxKeys=1,  # We only need to know if at least one exists
    )
    return "Y" if response.get("KeyCount", 0) > 0 else "N"


# -------------------------------------------------------------------------------------------------------------------
# def _list_zarr_files(bucket_name, prefix=""):
#     objects = []
#     paginator = s3_client.get_paginator("list_objects_v2")
#     page_iterator = paginator.paginate(Bucket=bucket_name, Prefix=prefix)

#     for page in page_iterator:
#         for obj in page.get("Contents", []):
#             objects.append(obj["Key"])

#     return objects

# def _get_zarr_files(bucket_name, prefix=""):
#     all_objects = _list_zarr_files(bucket_name, prefix)
#     zarr_files = [obj for obj in all_objects if obj.endswith('.zarr')]
#     return zarr_files
