import boto3
import pandas as pd
from io import BytesIO, StringIO
import numpy as np
import xarray as xr
import s3fs
from datetime import datetime

# Set environment variables
BUCKET_NAME = "wecc-historical-wx"
RAW_WX = "1_raw_wx/"
CLEAN_WX = "2_clean_wx/"
QAQC_WX = "3_qaqc_wx/"
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")
# goals from Jira ticket

# PART 1: function needs to be able to identify all qaqc per network stations lists
# piece needs to go through stationlist_NETWORK_qaqc.csv use parent s3 path from get_qaqc_stations work
# instead of searching for zarr file, would search for csv file

# this might be something that could be added to the stnlist_update_clean script, Victoria is on board, what does Neil think?

# PART 2: read all of the stations lists into memory, can be done with pandas

# PART 3: concatenate all of those per network stations lists into a single dataframe, also can be done with pandas,
# be sure to exclude VALLEYWATER, setting up s3 search piece might be the most difficult part, stnlist_update_clean might be helpful with this

# PART 4: sort alphabetically on network

# PART 5: export as a csv to s3 bucket!


def id_all_qaqc_networks(network: str):
    # if statement to search for csv file (can go through stnlist_update_qaqc for bones of this)
    parent_s3_path = f"{bucket_name}/{qaqc_wx}{network}"

    # Use s3fs to list all items under this path
    s3_fs = s3fs.S3FileSystem(anon=False)
    all_paths = s3_fs.ls(parent_s3_path)

    csv_folders = [f"{path}" for path in all_paths if path.endswith(".csv")]

    for item in csv_folders:
        if item.endswith(".nc"):  # Handle .nc files (non-csv)
            station_id = item.split(".")[-2].replace(".nc", "")
            df["ID"].append(
                station_id
            )  # <--- these won't be added, this is just adding code snippets for the moment
            df["Time_QAQC"].append("")  # Placeholder, to be resolved
            df["QAQC"].append("N")  # Not QA/QC passed
        elif item.endswith(".zarr"):
            # Extract the station ID from the folder name, which is usually the last part of the path
            station_id = item.split("/")[-1].split(".")[-2].replace(".zarr", "")
            df["ID"].append(station_id)
            df["Time_QAQC"].append("")  # Placeholder, to be resolved
            df["QAQC"].append("Y")  # QA/QC passed
        elif item.endswith(".csv"):
            station_id = item.split(".")[-2].replace(
                ".csv", ""
            )  # <-- this is also being updated, adding in as placeholder
        else:
            continue


# for loop to concatenate all of the per network stations lists into a single dataframe (done via pandas), excluding Valleywater
# QUESTION: there are older files that have currently been excluded from the cleaning stage, will these be excluded as well?
# QUESTION: is this going to be a one off tool or will it be used regularly? If regular, would be good in the post processing step (after merge step) and then decide where it fits
# otherwise, if this will not be used regularly, can have it as its own separate script
# sort the network alphabetically
# return at end will be exported csv to s3 bucket
