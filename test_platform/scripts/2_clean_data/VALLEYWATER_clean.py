"""
This script performs data cleaning for networks pulled through Valley Water API for historical observational analysis
in AR pattern recognition. These stations may be ingested into the Historical Data Platform in a future iteration, but
processing conforms to the HDP standard of quality.

Approach:
(1) Read through variables and drop unnecessary variables
(2) Infill mis-timed timesteps at the correct 15 min interval with NaN in precipitation record
(3) Add empty elevation variable, to be infilled via DEM
(4) Converts station metadata to standard format, with unique identifier
(5) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(6) Tracks existing qa/qc flag(s) for review
(7) Outputs cleaned variables as a single zarr for each station in an individual network.

Inputs: Raw data for the network's stations, with each csv file representing a station.
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as zarr.

Modification History:
- 11/26/2024: Converted from a python notebook to a python script
- 11/30/2024: Added empty elevation variable with proper attributes
- 12/6/2024: Script now creates csv file that stores station info & cleaning time & upload to s3. This is used in QAQC step 3.
- 12/12/2024: Change datetime conversion from using tz_localize and tz_convert to a simple +8 hr following advice from Valley Water team
- 01/14/2025:
    - Modified to work with new VW data that contains information about gaps.
    - Missing timesteps are filled with NaN instead of 0.
    - Time conversion removed because the new data contains a time column in UTC timezone
    - Data read in using pd.read_csv() instead of tempfile method
- 01/20/2025: Add check for station in API json file csv due to VW sending data that didn't exist in the API and error being raised in script (station ID 6053)
- 04/08/2025: Remove part of script that merge VW station list with existing full station list. This is now done via another script.

"""

## Imports
import xarray as xr
import pandas as pd
import numpy as np
import warnings
import sys
import os
import s3fs
import boto3
from datetime import datetime, timezone
from pathlib import Path  # to get file suffix
import time

# Set name of bucket and folder containing data
bucket = "wecc-historical-wx"
network = "VALLEYWATER"
folder_raw = "1_raw_wx/VALLEYWATER"  # Raw data
folder_clean = (
    "2_clean_wx/VALLEYWATER"  # Where to upload cleaned files. Folder must already exist
)


def main():
    """
    Main function that processes Valley Water precipitation data.

    This function:
    1. Reads raw precipitation data from S3
    2. Cleans and standardizes the data format
    3. Converts units and handles missing values
    4. Adds metadata and quality control flags
    5. Saves cleaned data as zarr files to S3
    6. Generates a station metadata CSV file

    """

    # For attributes of exported file
    timestamp = datetime.now(timezone.utc).strftime("%m-%d-%Y, %H:%M:%S")

    # ERA-formatted variable name
    era_var_name = "pr_15min"

    # Create empty dataframe for storing QAQC and station info
    # This will be saved as a csv file and used in QAQC step 3
    stations_df = pd.DataFrame(
        {
            "era-id": [],
            "longitude": [],
            "latitude": [],
            "elevation": [],
            "start-date": [],
            "end-date": [],
            "cleaned": [],
            "time_cleaned": [],
            "network": [],
            "{0}_nobs".format(era_var_name): [],
            "total_nobs": [],
        }
    )

    filenames = get_filenames_in_s3_folder(bucket=bucket, folder=folder_raw)
    filenames = [file for file in filenames if "Precip_Increm.Final@" in file]

    # Read in csv containing information about each station
    # Lat, lon, name, and watershed :)
    # This info was manually retrieved from the API
    # See scrape_json_file.ipynb for more info (in the valley water repo)
    station_info_csv_filepath = (
        "s3://wecc-historical-wx/1_raw_wx/VALLEYWATER/VALLEYWATER_station_info.csv"
    )
    station_info_df = pd.read_csv(station_info_csv_filepath)

    # Loop through each file, clean/reformat, and upload to s3
    # Print a pretty progress bar to console :)
    for i in progressbar(range(len(filenames))):
        filename = filenames[i]

        # Read in file
        df = pd.read_csv(
            filename, header=14
        )  # Remove header from each file so it can be read in as a dataframe

        # Convert time to datetime so it can be better managed programmatically
        # No time conversion needed since UTC column already exists in data
        # I have to do dt.tz_localize(None) because pandas does a weird thing where it adds the UTC timezone to the dtype of the object
        # xarray didn't know how to convert it to a datetime index when I converted it to a Dataset later in the pipeline
        df["time"] = pd.to_datetime(df["ISO 8601 UTC"]).dt.tz_localize(None)

        # Delete unused columns
        df = df[["time", "Value", "Approval Level"]]

        # Remove rows where Approval Level = Null
        # If theres a big chunk of missing data, only one NaN value is reported in the middle of that chunk
        # This messes up the temporal resolution of the data
        # We will remove that random NaN value in the middle of the chunk, which is identified by Approval Level = Null
        # Then, we will fill the entire missing data chunk with NaNs at the appropriate temporal resolution of the data
        df = df[~df["Approval Level"].isnull()]

        # Infill missing timestamps with -999
        # This helps us trace these values so that we can add the era QAQC flag to show that they were infilled with NaN
        # In a later step, we will replace -999 with NaN (after adding the QAQC flag)
        expected_tdelta = pd.Timedelta("0 days 00:15:00")  # 15 minutes
        df = df.set_index("time").resample(expected_tdelta).asfreq().fillna(value=-999)

        # Flag with ERA QAQC variable appropriate for NaN to correct timestamp
        # Where the Value is -999, fill the QAQC with the qaqc number
        # Where the value is NOT -999, fill the QAQC with NaN
        df[era_var_name + "_eraqc"] = (
            df["Value"]
            .where(df["Value"] == -999, other=np.NaN)
            .where(df["Value"] != -999, other=30)
        )

        # Also replace -999 in Values with NaN
        # Double checking on original NaNs not being accidentally flagged
        df = df.replace({"Value": {-999: np.NaN}})

        # Set Approval Level to appropriate values
        df = df.replace({"Approval Level": {-999: np.NaN}}).rename(
            columns={"Approval Level": "raw_qc"}
        )

        # Convert precip units in --> mm
        precip_col_name = "Value"  # Variable name in raw data
        df[era_var_name] = df[precip_col_name] * 25.4  # Convert in --> mm

        # Retrieve station ID from the filename
        station_id = int(filename.split("@")[1].split(".EntireRecord.csv")[0])

        # Get station info from csv file
        station_info_i = station_info_df[
            station_info_df["station_id"] == int(station_id)
        ]

        # Check if that station is in the csv! If not, print a warning, skip the station, and move to next loop iteration
        # Implementing this check due to station ID #6053 not found in API
        if len(station_info_i) < 1:
            print(
                "WARNING: no information found in csv file {0} for station ID {1}. Skipping and moving to next station".format(
                    station_info_csv_filepath, station_id
                )
            )
            continue  # End this loop iteration completely, continue to next station

        # Remove date strings from filename
        # Convert to uppercase characters to match existing filename formatting conventions
        # Current format: VALLEYWATER_[station_id].csv
        relative_filepath = "{}_{}.zarr".format(network, station_id)

        # Convert data type to xarray object
        df_pr = df[
            [era_var_name, era_var_name + "_eraqc", "raw_qc"]
        ]  # Just preserve the precip data
        ds = xr.Dataset.from_dataframe(df_pr)
        ds = ds.expand_dims(
            {"station": ["{}_{}".format(network, station_id)]}
        )  # Add singleton dimension for station network + id

        # Add lat and lon as coordinates
        ds = ds.assign_coords(
            {
                "lat": (
                    ("station", "time"),
                    np.full((1, len(df)), float(station_info_i["lat"].item())),
                ),
                "lon": (
                    ("station", "time"),
                    np.full((1, len(df)), float(station_info_i["lon"].item())),
                ),
            }
        )

        # Add elevation data filled with NaNs
        nan_data = np.full(
            ds["pr_15min"].shape, np.nan
        )  # Array of same shape as existing data variable
        ds = ds.assign(
            {"elevation": (ds.dims, nan_data)}
        )  # Add new data variable "elevation" to Dataset

        # Assign appropriate variable & coordinate attributes
        ds[era_var_name].attrs = {
            "long_name": "15_minute_precipitation_amount",
            "units": "mm/15min",
            "comment": "Precipitation accumulated in previous 15 minutes.",
        }
        ds["elevation"].attrs = {
            "standard_name": "height_above_mean_sea_level",
            "long_name": "station_elevation",
            "units": "meter",
            "positive": "up",
        }
        ds["station"].attrs = {"long_name": "station_id"}
        ds["time"].attrs = {
            "long_name": "time",
            "standard_name": "time",
            "comment": "in UTC",
        }
        ds["lat"].attrs = {
            "long_name": "latitude",
            "standard_name": "latitude",
            "units": "degrees_north",
        }
        ds["lon"].attrs = {
            "long_name": "longitude",
            "standard_name": "longitude",
            "units": "degrees_east",
        }

        # Assign appropriate global attributes
        ds.attrs = {
            "title": "{} cleaned".format(network),
            "institution": "Eagle Rock Analytics",
            "license": "",
            "source": "",
            "station_name": station_info_i["station_name"].item(),
            "watershed": station_info_i["watershed"].item(),
            "Networks": "VALLEYWATER",
            "sensor_height_m": np.nan,
            "disclaimer": "This document was prepared as a result of work funded by the Santa Clara Valley Water District. It does not necessarily represent the views of the Santa Clara Valley Water District or its employees. Neither the Santa Clara Valley Water District, nor it's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by Santa Clara Valley Water District, nor has the Santa Clara Valley Water District passed upon the accuracy of the information in this document.",
            "history": "VALLEYWATER_clean.py script run on {} UTC".format(timestamp),
            "comment": "Intermediate data product: may not have been subject to any cleaning or QA/QC processing",
            "raw_files_merged": 1,
        }

        # Set the complete filepath to the s3 bucket
        filepath_s3 = "s3://{}/{}/{}".format(bucket, folder_clean, relative_filepath)

        # Upload as zarr to bucket
        ds.to_zarr(
            filepath_s3,
            consolidated=True,  # https://docs.xarray.dev/en/stable/internals/zarr-encoding-spec.html
            mode="w",  # Write & overwrite if file with same name exists already
        )

        # Now, save info for this station to stations dataframe
        # Info for each station is saved as a single row, which is then appended as a row to the master dataframe
        time = pd.to_datetime(ds.time.values)
        nobs = int(len(df[era_var_name]))  # Number of observations
        stations_df_i = pd.DataFrame(
            {
                "era-id": ["{0}_{1}".format(network, station_id)],
                "longitude": [station_info_i["lon"].item()],
                "latitude": [station_info_i["lat"].item()],
                "elevation": [np.nan],
                "start-date": [time[0].strftime("%m-%d-%Y, %H:%M:%S")],
                "end-date": [time[-1].strftime("%m-%d-%Y, %H:%M:%S")],
                "cleaned": ["Y"],
                "time_cleaned": [timestamp],
                "network": [network],
                "{0}_nobs".format(era_var_name): [nobs],
                "total_nobs": [
                    nobs
                ],  # Total nobs is the same as single variable nobs because we just have one variable
            }
        )
        stations_df = pd.concat([stations_df, stations_df_i], ignore_index=True)

        # Re-sort into alphabetical/numerical order
        stations_df = stations_df.sort_values("era-id", ignore_index=True)

    print(
        "Completed uploading all zarrs for VALLEYWATER QAQC step 2. Saving csv station list..."
    )

    # Save csv file with station info to AWS
    csv_s3_filepath = "s3://{0}/{1}/stationlist_VALLEYWATER_cleaned.csv".format(
        bucket, folder_clean
    )
    stations_df.to_csv(csv_s3_filepath, index=False)
    print("Station list csv saved to s3 path: {0}".format(csv_s3_filepath))

    print("SCRIPT COMPLETE")


def get_filenames_in_s3_folder(bucket: str, folder: str) -> list[str]:
    """
    Get a list of files in s3 bucket.
    Make sure you follow the naming rules exactly for the two function arguments.
    See example in the function docstrings for more details.

    Parameters
    ---------
    bucket: str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc.
    folder: str
        Folder within the bucket that you want the filenames from

    Returns
    -------
    files_in_s3: list of str
        List of filenames in the bucket

    Example
    -------
    You want to get all the filenames in a s3 bucket with the following path:
    s3 URI: "s3://wecc-historical-wx/1_raw_wx/VALLEYWATER/"
    >>> get_filenames_in_s3_folder(
    >>>    bucket = "wecc-historical-wx",
    >>>    folder = "1_raw_wx/VALLEYWATER"
    >>> )
    ['ValleyWater_6001_1900-01-01_2024-11-11.csv','ValleyWater_6004_1900-01-01_2024-11-11.csv']

    References
    ----------
    https://stackoverflow.com/questions/59225939/get-only-file-names-from-s3-bucket-folder

    """

    if folder[:-1] != "/":  # add slash
        folder = folder + "/"

    s3 = boto3.resource("s3")
    s3_bucket = s3.Bucket(bucket)

    # Get all the filenames
    files_in_s3 = [f.key for f in s3_bucket.objects.filter(Prefix=folder).all()]

    # Get filenames with URI
    # This allows you to load the file directly using the path
    filenames_with_uri = ["s3://{}/{}".format(bucket, file) for file in files_in_s3]

    return filenames_with_uri


if __name__ == "__main__":
    with warnings.catch_warnings():
        # Annoying RuntimeWarning is raised... ignore here so we can better evaluate the print output of the function, which includes a progress bar
        warnings.simplefilter("ignore")
        main()
