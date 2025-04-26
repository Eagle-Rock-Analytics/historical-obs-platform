"""
CW3E_combine.py

This script merges individual NetCDF files for the same weather station in the CW3E network.
Originally, CW3E station data were split into multiple NetCDF files due to size limitations during the initial cleaning step.
This script combines them into single zarr files for easier ingestion and analysis.

Overview:
---------
1. Read the list of cleaned stations from a CSV stored in S3.
2. Identify all NetCDF files corresponding to each station.
3. Open all NetCDF files, concatenate them along the "time" dimension, and sort chronologically.
4. Save the merged dataset as a zarr file back to S3.

Outputs:
--------
- One consolidated zarr file per station, stored in the same S3 folder.

Example usage:
--------------
python CW3E_combine.py
"""

import xarray as xr
import s3fs
import pandas as pd
import boto3
from time import time


def get_filenames_in_s3_folder(bucket, folder):
    """Get a list of files in s3 bucket.
    Make sure you follow the naming rules exactly for the two function arguments.
    See example in the function docstrings for more details.

    Parameters
    ---------
    bucket : str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc...
    folder : str
        Folder within the bucket that you want the filenames from
        MAKE SURE folder doesn't have a trailing "/"
        i.e. it should be "[folder]", not "[folder]/"

    Returns
    -------
    files_in_s3 : list of str
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
    [1] https://stackoverflow.com/questions/59225939/get-only-file-names-from-s3-bucket-folder
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


def open_multiple_netcdf_files(filepaths, fs):
    """Open a list of NetCDF files from S3 using xarray.

    Parameters
    ----------
    filepaths : list of str
        List of S3 filepaths to NetCDF files.
    fs : s3fs.S3FileSystem
        Open S3 filesystem object.

    Returns
    -------
    ds_list : list of xarray.Dataset
        List of loaded xarray datasets.
    """
    ds_list = []
    for filepath in filepaths:
        try:
            with fs.open(filepath) as fileObj:
                ds = xr.open_dataset(fileObj).load()
                ds_list.append(ds)
        except Exception:
            print(f"File {filepath} could not be read in")
            continue
    return ds_list


def concat_and_sort_datasets(ds_list):
    """Concatenate a list of datasets along 'time' and sort chronologically.

    Parameters
    ----------
    ds_list : list of xarray.Dataset
        List of datasets to concatenate.

    Returns
    -------
    station_ds : xarray.Dataset
        Concatenated and sorted dataset.
    """
    if len(ds_list) == 0:
        return None
    station_ds = xr.concat(ds_list, dim="time").sortby("time")
    return station_ds


def main():
    """Main processing workflow."""

    start_time = time()

    # Initialize an S3 file system client to interact with S3 storage
    fs = s3fs.S3FileSystem()

    # Define s3 paths and such
    bucket = "wecc-historical-wx"
    CW3E_cleaned_folder = "2_clean_wx/CW3E"
    csv_filepath_s3 = (
        "s3://wecc-historical-wx/2_clean_wx/temp_clean_all_station_list.csv"
    )

    # Read list of cleaned stations from CSV
    stations_df = pd.read_csv(csv_filepath_s3)

    # Filter for CW3E network and cleaned stations
    network_df = stations_df[
        (stations_df["network"] == "CW3E") & (stations_df["cleaned"] == "Y")
    ]

    # Get list of all CW3E NetCDF files in S3
    cw3e_files_all = get_filenames_in_s3_folder(bucket, CW3E_cleaned_folder)
    cw3e_nc_files = [file for file in cw3e_files_all if file.split(".")[1] == "nc"]

    # Get the IDs of all the stations in CW3E
    station_ids = network_df["era-id"].values

    print(f"Processing {len(station_ids)} stations for network: CW3E")
    for i in range(len(station_ids)):
        station_id = station_ids[i]
        print(f"Processing station {station_id}: station {i+1}/{len(station_ids)}")

        # Find all NetCDF files for the current station
        filenames_in_s3 = [file for file in cw3e_nc_files if station_id in file]
        if len(filenames_in_s3) == 0:
            print(f"No NetCDFs found for station: {station_id}")
            continue
        else:
            print(f"Found {len(filenames_in_s3)} files: {filenames_in_s3}")

        filepaths_in_s3 = [
            f"s3://{bucket}/{CW3E_cleaned_folder}/{filename}"
            for filename in filenames_in_s3
        ]

        # Open and load all datasets
        print("Opening all files... ")
        ds_list = open_multiple_netcdf_files(filepaths_in_s3, fs)
        if len(ds_list) == 0:
            continue

        # Concatenate and sort datasets
        print(
            "Concatenating files along time dimension, ensuring chronological sorting... "
        )
        station_ds = concat_and_sort_datasets(ds_list)
        if station_ds is None:
            continue

        # Write to zarr format
        print("Writing to zarr... ")
        zarr_s3_path = f"s3://{bucket}/{CW3E_cleaned_folder}/{station_id}.zarr"
        try:
            station_ds.to_zarr(
                zarr_s3_path,
                mode="w",
                consolidated=True,
            )
            print(f"Zarr written to: {zarr_s3_path}")
        except Exception:
            print(
                f"Zarr for station {station_id} could not be successfully written to S3 bucket"
            )
            continue

    print("Script complete!")
    end_time = time()
    elapsed_time = (end_time - start_time) / 60  # divide by 60 to get minutes
    print(f"Total elapsed time: {elapsed_time:.2f} minutes")


if __name__ == "__main__":
    main()
