"""
homr_metadata.py

This script queries the HOMR metadata API and returns information about relevant stations.

Functions
---------
- get_homr_metadata: takes an NCDC ID and returns up to 6 tables with metadata information about the station.
- get_all_homr_ids: returns a table with all NCDC IDs in WECC states and high-level metadata (e.g. location, state, network) for each station, saving this table to AWS
- get_all_homr_metadata: reads in this table and queries all NCDC IDs in these states for detailed information, including sensor maintenance logs, obstructions
    and detailed remarks. These are saved as 6 tables, linked through the shared NCDC ID.

Intended Use
-------------
Generates metadata for relevant stations. Run after a full pull. 

References
----------
https://www.ncei.noaa.gov/access/homr/api
"""

import requests
import pandas as pd
from io import StringIO, BytesIO
import boto3

s3_cl = boto3.client("s3")


def flatten_data(y: dict | list) -> dict:
    """
    Flattens multi-level nested JSONs, splitting columns composed of both lists and dictionaries.

    Parameters
    ----------
    y : (dict | list)
        the nested JSON data to flatten

    Returns
    -------
    out :  dict
        a flattened dictionary with combined keys

    References
    ----------
    https://gist.github.com/davidwarshaw/8d4c74c31371f8f8d5eb00cccb86dd08
    """

    out = {}

    def flatten(x, name=""):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + "_")
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + "_")
                i += 1
        else:
            out[name[:-1]] = x

    flatten(y)
    return out


def get_homr_metadata(
    id: str,
) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    """
    Gets metadata for a given NCDC ID from the HOMR API

    Parameters
    ----------
    id : str
        The NCDC ID

    Returns
    -------
    names : pd.DataFrame
        station names
    identifiers : pd.DataFrame
        station identification codes
    platforms : pd.DataFrame
        platform
    locations : pd.DataFrame
        geographic location details
    remarks : pd.DataFrame
        any remarks listed (may be empty)
    updates : pd.DataFrame
        any updates listed (may be empty)
    """

    # This function takes one NCDC ID (string) as input and returns the station's names, identifiers, platforms, location, remarks, and updates
    # as distinct objects.
    testurl = f"https://www.ncei.noaa.gov/access/homr/services/station/{id}"
    request = requests.get(testurl).json()
    for i in request["stationCollection"]["stations"]:
        # print(i)
        # flat = flatten_data(i) # Ignore dictionaries, keep each row as a station
        pandas = pd.json_normalize(i, max_level=0)

        # Split into 4 tables:
        # Names
        names = pd.json_normalize(pandas["names"][0])
        names.insert(loc=0, column="ncdcid", value=id)

        # Identifiers
        identifiers = pd.json_normalize(pandas["identifiers"][0])
        identifiers.insert(loc=0, column="ncdcid", value=id)

        # Platforms
        if "platforms" in pandas.columns:
            platforms = pd.json_normalize(pandas["platforms"][0])
            platforms.insert(loc=0, column="ncdcid", value=id)
        else:
            platforms = pd.DataFrame()

        # Location
        location = flatten_data(pandas["location"][0])
        location = pd.json_normalize(location)
        location.insert(loc=0, column="ncdcid", value=id)

        # Remarks
        if "remarks" in pandas.columns:
            remarks = pd.json_normalize(pandas["remarks"][0])
            remarks.insert(loc=0, column="ncdcid", value=id)
        else:
            remarks = pd.DataFrame()

        # Updates
        if "updates" in pandas.columns:
            updates = pd.json_normalize(pandas["updates"][0])
            updates.insert(loc=0, column="ncdcid", value=id)
        else:
            updates = pd.DataFrame()

    return names, identifiers, platforms, location, remarks, updates


def get_all_homr_ids(bucket_name: str, savedir: str):
    """
    Iterates through all WECC states and save header HOMR metadata for all stations.

    Parameters
    ----------
    bucket_name : str
        name of the S3 bucket where HOMR metadata will be stored
    savedir : str
        directory path in the bucket where 'homr_ids.csv' file is saved to QAQC folder

    Returns
    -------
    None
    """

    dfs = []
    states = [
        "WA",
        "OR",
        "CA",
        "ID",
        "NV",
        "AZ",
        "UT",
        "NM",
        "CO",
        "TX",
        "WY",
        "MT",
        "SD",
    ]
    for state in states:
        testurl = f"https://www.ncei.noaa.gov/access/homr/services/station/search?headersOnly=true&state={state}"
        try:
            request = requests.get(testurl).json()
            df = pd.json_normalize(request["stationCollection"]["stations"])
            df.columns = [col.replace("header.", "") for col in df.columns]
            df.columns = [col.replace("por.", "") for col in df.columns]
            dfs.append(df)
            # print(df)
        except Exception as e:
            print(f"Unable to process {state}. {e}")
            continue

    homr_df = pd.concat(dfs)
    # Save station list to AWS

    csv_buffer = StringIO()
    homr_df.to_csv(csv_buffer, index=False)
    content = csv_buffer.getvalue()

    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=savedir + "homr_ids.csv")


def get_all_homr_metadata(bucket_name: str, savedir: str):
    """
    Takes all NCDC IDs saved in homr_ids.csv and compiles 5 csvs of names, identifiers, platforms, location, remarks, and updates.

    Parameters
    ----------
    bucket_name : str
        S3 bucket where files will be saved
    savedir : str
        Directory path in the bucket

    Returns
    -------
    None
    """

    # initialize dfs
    namedf = []
    identifierdf = []
    platformdf = []
    locationdf = []
    remarkdf = []
    updatedf = []

    # Read in homr_ids.csv
    obj = s3_cl.get_object(Bucket=bucket_name, Key=savedir + "homr_ids.csv")
    homr_ids = pd.read_csv(BytesIO(obj["Body"].read()))
    count = 0
    for i in homr_ids["ncdcStnId"]:
        count += 1
        names, identifiers, platforms, location, remarks, updates = get_homr_metadata(i)
        namedf.append(names)
        identifierdf.append(identifiers)
        platformdf.append(platforms)
        locationdf.append(location)
        remarkdf.append(remarks)
        updatedf.append(updates)
        print(count)

        # Print progress statement.
        if count % 250 == 0:  # For every two hundred and fifty records
            print(
                f"Metadata for {count} of {len(homr_ids.ncdcStnId)} stations processed."
            )

    namedf = pd.concat(namedf)
    identifierdf = pd.concat(identifierdf)
    platformdf = pd.concat(platformdf)
    locationdf = pd.concat(locationdf)
    remarkdf = pd.concat(remarkdf)
    updatedf = pd.concat(updatedf)

    # Save to AWS
    for dic in [
        {"homr_name.csv": namedf},
        {"homr_identifier.csv": identifierdf},
        {"homr_platform.csv": platformdf},
        {"homr_location.csv": locationdf},
        {"homr_remark.csv": remarkdf},
        {"homr_update.csv": updatedf},
    ]:
        for key, val in dic.items():
            csv_buffer = StringIO()
            val.to_csv(csv_buffer, index=False)
            content = csv_buffer.getvalue()
            s3_cl.put_object(Bucket=bucket_name, Body=content, Key=savedir + key)


if __name__ == "__main__":
    import os, sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from paths import BUCKET_NAME, QAQC_WX

    SAVEDIR = f"{QAQC_WX}/HOMR-Metadata/"

    # get_homr_metadata('20027492') # Run for one station
    get_all_homr_metadata(BUCKET_NAME, SAVEDIR)  # Only run periodically.
