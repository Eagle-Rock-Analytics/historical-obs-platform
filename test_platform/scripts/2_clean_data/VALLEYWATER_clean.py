""" Valley Water Data: QAQC 2
1. Download raw data from AWS bucket 
2. Infill missing timesteps with 0
3. Convert time PST --> UTC 
4. Add empty elevation variable
5. Reformat data, add attributes, etc. 
6. Upload cleaned file to s3 as a zarr store 

Author: Nicole Keeney
Creation Date: 11/18/2024
Modification History: 
- 11/26/2024: Converted from a python notebook to a python script
- 11/30/2024: Added empty elevation variable with proper attributes
- 12/6/2024: Script now creates csv file that stores station info & cleaning time & upload to s3. This is used in QAQC step 3.  
- 12/12/2024: Change datetime conversion from using tz_localize and tz_convert to a simple +8 hr following advice from Valley Water team
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
import tempfile

# Set name of bucket and folder containing data
bucket = "wecc-historical-wx"
folder_raw = "1_raw_wx/VALLEYWATER"  # Raw data
folder_clean = (
    "2_clean_wx/VALLEYWATER"  # Where to upload cleaned files. Folder must already exist
)


def main():
    # Name of variable
    var_name = "pr_15min"

    # For attributes of netCDF file.
    timestamp = datetime.now(timezone.utc).strftime("%m-%d-%Y, %H:%M:%S")

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
            "{0}_nobs".format(var_name): [],
            "total_nobs": [],
        }
    )

    # Define temporary directory in local drive for downloading data from S3 bucket
    # If the directory doesn't exist, it will be created
    temp_dir = "../../data/tmp"
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # Read csv file containing station data from s3 bucket
    station_info_filename = "VALLEYWATER_station_info.csv"
    fs = s3fs.S3FileSystem()
    with fs.open(
        "s3://{0}/{1}/{2}".format(bucket, folder_raw, station_info_filename)
    ) as f:
        station_info_df = pd.read_csv(f)

    # Get all the filenames
    filenames = get_filenames_in_s3_folder(bucket=bucket, folder=folder_raw)
    filenames = [
        file for file in filenames if "station_info" not in file
    ]  # Ignore station info csv file

    # Loop through each file, clean/reformat, and upload to s3
    # Print a pretty progress bar to console :)
    for i in progressbar(range(len(filenames))):
        filename = filenames[i]

        # Read in and clean the data
        df = read_and_clean_data(
            bucket=bucket,
            folder=folder_raw,
            filename=filename,
            temp_dir=temp_dir,
            delete=True,
        )

        # Get station info from filename
        network, station_id, time_start, time_end = filename.split("_")
        network = (
            network.upper()
        )  # Capitalize the network following hist-obs conventions

        # Get station info from csv file
        station_info_i = station_info_df[
            station_info_df["station_id"] == int(station_id)
        ]

        # Convert precip units in --> mm
        precip_col_name = "Precipitation (in.)"
        if precip_col_name not in list(df.columns):  # Ensure that the column exists!
            raise ValueError(
                "'{}' expected as a column name in station data but not found".format(
                    precip_col_name
                )
            )
        df = df.rename(
            columns={precip_col_name: var_name}
        )  # Rename column following naming conventions
        df[var_name] = df[var_name] * 25.4  # Convert in --> mm

        # Convert data type to xarray object
        ds = xr.Dataset.from_dataframe(df)
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
        ds[var_name].attrs = {
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
        # Remove date strings from filename
        # Convert to uppercase characters to match existing filename formatting conventions
        # Current format: VALLEYWATER_[station_id].csv
        relative_filepath = "{}_{}.zarr".format(network, station_id)

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
        nobs = int(len(df[var_name]))  # Number of observations
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
                "{0}_nobs".format(var_name): [nobs],
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

    # Merge VW station list with existing full station list
    full_station_list_s3_path = (
        "s3://wecc-historical-wx/2_clean_wx/temp_clean_all_station_list.csv"
    )
    print(
        "Adding VW station data to full station list located at {0}".format(
            full_station_list_s3_path
        )
    )
    stations_all_df = pd.read_csv(full_station_list_s3_path, index_col=0)

    # Concatenate VW data to main dataframe
    stations_df_new = pd.concat([stations_df, stations_all_df], ignore_index=True)

    # Re-sort into alphabetical order, just in case!
    stations_df_new = stations_df_new.sort_values("era-id", ignore_index=True)

    # Move nobs back to final column position
    # When merging, it gets relocated
    stations_df_new["total_nobs"] = stations_df_new.pop("total_nobs")

    # Write to s3
    stations_df_new.to_csv(full_station_list_s3_path)
    print("SUCCESS: csv file updated with VW info in s3")

    print("SCRIPT COMPLETE")


def read_and_clean_data(
    bucket, folder, filename, temp_dir="", tdelta="0h15min0s", delete=True
):
    """
    1. Read in file from s3 bucket
    2. Infill missing timestamps with 0
    3. Rename/reformat columns

    Parameters
    ----------
    bucket: str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc...
    folder: str
        Folder within the bucket containing the data you want to read
    filename: str
        Name of the file in the bucket (i.e. "ValleyWater_6004_1900-01-01_2024-11-11.csv")
    temp_dir: str, optional
        Temporary directory for storing downloaded data
        Default to current directory
    tdelta: str, optional
        Time delta between station observations
        String must be formatted as such: "[n]h[x]min[y]s" where n,x,y represent the hours, minutes, and seconds
        Default to "0h15min0s" (15 minutes)
        Function will check that the station minimum time delta == tdelta
    delete: boolean, optional
        Delete the local version of the data after downloading?
        Default to True

    Returns
    -------
    df: pd.DataFrame

    """

    # Read raw data from bucket
    df = read_data_from_s3(
        bucket=bucket, folder=folder, filename=filename, temp_dir=temp_dir
    )

    # Convert "Timestamp" column to column "Time"
    # Convert values to datetime object
    # Convert PST to UTC (add 8 hr)
    df["time"] = pd.to_datetime(df["Timestamp"].values) + pd.DateOffset(hours=8)
    df.drop("Timestamp", axis=1, inplace=True)

    # Check that minimum time delta between observations is the same as the function argument for tdelta
    timedelta_row = df["time"].diff()
    min_tdelta = strftdelta(timedelta_row.min())
    if (min_tdelta) != tdelta:
        print(
            "WARNING: Your input timedelta is not the same as the minimum timedelta found in the data. \nYour input timedelta: tdelta = {0}\nMinimum timedelta found in the data: {1}".format(
                tdelta, min_tdelta
            )
        )

    # Infill missing timestamps with 0
    df = df.set_index("time").resample(min_tdelta).asfreq().fillna(value=0)

    return df


def strftdelta(tdelta):
    """Get timedelta string for pd.Timedelta object

    Arguments
    ---------
    tdelta: pd.Timedelta

    Returns
    -------
    str: timedelta in string format
        Format is "[n]h[x]min[y]s" where n,x,y represent the hours, minutes, and seconds
        This format can be used by the pd.resample function

    """
    d = {"days": tdelta.days}
    d["hours"], rem = divmod(tdelta.seconds, 3600)
    d["minutes"], d["seconds"] = divmod(rem, 60)
    return "{}h{}min{}s".format(d["hours"], d["minutes"], d["seconds"])


def get_filenames_in_s3_folder(bucket, folder):
    """Get a list of files in s3 bucket.
    Make sure you follow the naming rules exactly for the two function arguments.
    See example in the function docstrings for more details.

    Parameters
    ---------
    bucket: str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc...
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

    s3 = boto3.resource("s3")
    s3_bucket = s3.Bucket(bucket)

    # Get all the filenames
    # Just get relative path (f.key.split(folder + "/")[1])
    files_in_s3 = [
        f.key.split(folder + "/")[1]
        for f in s3_bucket.objects.filter(Prefix=folder).all()
    ]

    # Delete empty filenames
    # I think the "empty" filename/s is just the bucket path, which isn't a file but is read as an object by the objects.filter function
    files_in_s3 = [f for f in files_in_s3 if f != ""]

    return files_in_s3


def read_data_from_s3(bucket, folder, filename, temp_dir=""):
    """Read data from AWS s3 bucket

    Parameters
    ----------
    bucket: str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc...
    folder: str
        Folder within the bucket containing the data you want to read
    filename: str
        Name of the file in the bucket (i.e. "ASOSAWOS_72012200114.nc")
    temp_dir: str, optional
        Temporary directory for storing downloaded data
        Default to current directory

    Returns
    -------
    station_data: xr.Dataset (data_type = "netcdf") or pd.DataFrame (data_type = "csv")

    Example
    -------
    Proper arguments for data with the s3 URI: s3://wecc-historical-wx/1_raw_wx/VALLEYWATER/ValleyWater_5007_1980-01-01_2024-11-08.csv
    >>> read_data_from_s3(
    >>>     bucket = "wecc-historical-wx",
    >>>     folder = "1_raw_wx/VALLEYWATER",
    >>>     filename = "ValleyWater_6001_1900-01-01_2024-11-11.csv"
    >>> )

    """
    # Get suffix of filename (i.e. .csv, .nc)
    suffix = Path(filename).suffix

    # Temp file for downloading from s3
    temp_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, prefix="", suffix=suffix, delete=True
    )

    # Create s3 file system
    s3 = s3fs.S3FileSystem(anon=False)

    # Get URL to netcdf in S3
    s3_url = "s3://{}/{}/{}".format(bucket, folder, filename)

    # Read data
    s3_file_obj = s3.get(s3_url, temp_file.name)

    if suffix in [".nc", ".h5"]:
        station_data = xr.open_dataset(temp_file.name, engine="h5netcdf").load()
    elif suffix == ".csv":
        station_data = pd.read_csv(temp_file.name)

    # Close temporary file
    temp_file.close()

    return station_data


def progressbar(it, prefix="", size=60, out=sys.stdout):
    """
    Print a progress bar to console

    Parameters
    ----------
    it: int
        iternation of list
    size: int, optional
        size (length) of progress bar

    Returns
    -------
    progress bar printed to console

    Example
    -------
    >>> for i in progressbar(10): # Progress bar of length 10 is printed after each iteration i
    >>> # Loop does something

    References
    ----------
    https://stackoverflow.com/questions/3160699/python-progress-bar

    """
    count = len(it)
    start = time.time()  # time estimate start

    def show(j):
        x = int(size * j / count)
        # time estimate calculation and string
        remaining = ((time.time() - start) / j) * (count - j)
        mins, sec = divmod(remaining, 60)  # limited to minutes
        time_str = f"{int(mins):02}:{sec:03.1f}"
        print(
            f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} Est wait {time_str}",
            end="\r",
            file=out,
            flush=True,
        )

    show(0.1)  # avoid div/0
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    print("\n", flush=True, file=out)


if __name__ == "__main__":
    with warnings.catch_warnings():
        # Annoying RuntimeWarning is raised... ignore here so we can better evaluate the print output of the function, which includes a progress bar
        warnings.simplefilter("ignore")
        main()
