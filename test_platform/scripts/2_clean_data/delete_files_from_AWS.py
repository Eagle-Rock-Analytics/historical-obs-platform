"""
delete_files_from_AWS.py

This script deletes files in batch from AWS. Options for which files to delete are provided.

Functions
---------
- delete_files_from_AWS: Batch delete files from AWS.

Intended Use
------------
Script is intended to be used CAREFULLY and on a rare as-needed basis. 
User is forced to input specific network and Y/N answer in order to proceed. 
Options:
    - which_to_delete = "empty"
        When grabbing data for an updated pull from Synoptic/MADIS networks for end of 2022 data,
        our update_pull script produced an update raw data file for every station, regardless if
        that station had new data to add. Therefore, there are many "empty" raw data files in our
        /1_raw_wx bucket that will cause issues in cleaning and updating the cleaned station list
        for cleaned variable coverage.
        This option *specifically* identifies these empty files and deletes them.

    - which_to_delete = "update_pull"
        When grabbing data for an updated pull from Synoptic/MADIS networks for end of 2022 data,
        our update_pull script produced an update raw data file for every station. This file has a
        suffix with the date of the update pull. When the next full pull is performed, the original
        raw data files will be overwritten, but the old update pull files with the date suffix will
        not, even though that data should be encapsulated in the latest new full pull.
        This option *specifically* identifies the files with a date suffix and deletes them.
        Only perform after a new full pull has occured.

    - which_to_delete = "monthly"
        NDBC, MARITIME networks produce a monthly file if the year of data coverage is within the
        current year, rather than a single annual file. Once the annual file is available for the year,
        these monthly files should be deleted.

    - which_to_delete = "no_wx_data"
        Some networks (e.g., NDBC, MARITIME, ASOSAWOS) do not have any valid meteorological variables
        and should not proceed through cleaning. 
"""

import pandas as pd
import s3fs
import boto3
from clean_utils import get_file_paths

s3 = boto3.resource("s3")
BUCKET_NAME = "wecc-historical-wx"


def delete_files_from_AWS(network: str, which_to_delete: str):
    """
    Batch delete files form AWS.
    This function should be used CAREFULLY as it will PERMANENTLY batch delete files from AWS!!!

    Parameters
    ----------
    network : str
        name of network
    which_to_delete : str
        batch delete files from AWS options: "empty", "update_pull", "monthly", "no_wx_data" (see notes)

    Returns
    -------
    None

    Notes
    -----
    Options of which files to delete:
    which_to_delete = "empty"
        When grabbing data for an updated pull from Synoptic/MADIS networks for end of 2022 data,
        our update_pull script produced an update raw data file for every station, regardless if
        that station had new data to add. Therefore, there are many "empty" raw data files in our
        /1_raw_wx bucket that will cause issues in cleaning and updating the cleaned station list
        for cleaned variable coverage.
        This option *specifically* identifies these empty files and deletes them.
        Note: we are no longer pulling data from Synoptic in the future.

    which_to_delete = "update_pull"
        When grabbing data for an updated pull from Synoptic/MADIS networks for end of 2022 data,
        our update_pull script produced an update raw data file for every station. This file has a
        suffix with the date of the update pull. When the next full pull is performed, the original
        raw data files will be overwritten, but the old update pull files with the date suffix will
        not, even though that data should be encapsulated in the latest new full pull.
        This option *specifically* identifies the files with a date suffix and deletes them.
        Only perform after a new full pull has occured.

    which_to_delete = "monthly"
        NDBC, MARITIME networks produce a monthly file if the year of data coverage is within the
        current year, rather than a single annual file. Once the annual file is available for the year,
        these monthly files should be deleted.

    which_to_delete = "no_wx_data"
        Some networks (e.g., NDBC, MARITIME, ASOSAWOS) do not have any valid meteorological variables
        and should not proceed through cleaning.

    If which_to_delete is not provided an option, nothing will delete.

    References
    ----------
    [1] Modified from: https://stackoverflow.com/questions/32501995/boto3-s3-renaming-an-object-using-copy-object
    """

    # deletes only update pull files with date extensions that are empty
    if which_to_delete == "empty":
        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        # get files
        all_files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=rawdir):
            file = str(item.key)
            all_files += [file]

        # identify files with date in filename
        files_to_delete = []
        for file in all_files:
            file_parts = file.split("/")
            stn_file = file_parts[-1]  # grabs just file name, no path
            stn_file = stn_file.split(".")[0]  # removes file extenstion

            if (
                ("_" in stn_file)
                and "errors" not in stn_file
                and "stationlist" not in stn_file
            ):
                fn_date = file.split("_")[-1]
                # grabs part of filename after _ character
                if "-" in fn_date:
                    # must follow YYYY-MM-DD format
                    files_to_delete += [file]

        # open each file, and assess if empty
        empty = []
        for file in files_to_delete:
            print(f"Checking {file} for empty size")
            fs = s3fs.S3FileSystem()
            aws_url = f"s3://{BUCKET_NAME}/{file}"

            with fs.open(aws_url) as fileObj:
                df = pd.read_csv(fileObj, header=6, low_memory=False)

                # assess if empty
                if len(df.index) == 1 and df["Station_ID"].isnull().values.any():
                    # any station that is missing the station_id is empty
                    empty += [file]

        print(f"Number of empty files to delete in {rawdir}: {len(empty)}")
        print(empty)

        resp = input(
            f"\nWARNING: this will permanently delete {len(empty)} files from {rawdir}. Are you ready to delete these files? (Y/N)"
        )
        if resp == "N":
            # do not delete files!!
            print("Files not deleted. Exiting function.")

        elif resp == "Y":
            # proceed to delete files
            for file in empty:
                # delete file -- THIS IS THE LINE OF CODE THAT MATTERS
                s3.Object(BUCKET_NAME, file).delete()
                print(f"File {file} deleted from AWS bucket")

        else:
            # any other character thrown in, do not delete
            print("Invalid response. Exiting function.")

    # deletes ALL update pull files with date extensions
    elif which_to_delete == "update_pull":
        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        # get files
        all_files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=rawdir):
            file = str(item.key)
            all_files += [file]

        # identify files with date in filename
        files_to_delete = []
        for file in all_files:
            file_parts = file.split("/")
            stn_file = file_parts[-1]  # grabs just file name, no path
            stn_file = stn_file.split(".")[0]  # removes file extenstion

            if (
                ("_" in stn_file)
                and "errors" not in stn_file
                and "stationlist" not in stn_file
            ):
                # grabs part of filename after _ character
                fn_date = file.split("_")[-1]
                if "-" in fn_date:
                    # must follow YYYY-MM-DD format
                    files_to_delete += [file]

        print(f"Number of files to delete in {rawdir}: {len(files_to_delete)}")
        print(files_to_delete)  # should have date in filename

        resp = input(
            f"\nWARNING: this will permanently delete {len(files_to_delete)} files from {rawdir}. Are you ready to delete these files? (Y/N)"
        )
        if resp == "N":
            # do not delete files!!
            print("Files not deleted. Exiting function.")

        elif resp == "Y":
            # proceed to delete files
            for file in files_to_delete:
                # delete file -- THIS IS THE LINE OF CODE THAT MATTERS
                s3.Object(BUCKET_NAME, file).delete()
                print(f"File {file} deleted from AWS bucket")

        else:
            # any other character thrown in, do not delete
            print("Invalid response. Exiting function.")

    # deletes the monthly files from NDBC, MARITIME after year file is in
    # note: CIMIS has monthly files, but should not delete, as these files retain data for all stations
    elif which_to_delete == "monthly":
        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        if network == "NDBC" or network == "MARITIME":
            # get files
            all_files = []
            for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=rawdir):
                file = str(item.key)
                all_files += [file]

            all_files = [file for file in all_files if "stationlist" not in file]
            all_files = [file for file in all_files if "error" not in file]

            # identify monthly files
            files_to_delete = []
            for file in all_files:
                file_parts = file.split("/")
                stn_file = file_parts[-1]  # grabs just file name, no path
                stn_file = stn_file.split(".")[0]
                # removes file extenstion (.txt.gz, or .zip)

                if "csv" in stn_file:
                    # skip canadian files, they're annual
                    continue

                else:
                    stn_year = stn_file[-4]
                    # use to identify if annual file is already downloaded
                    stn_name = stn_file[:5]
                    # NDBC and MARITIME station names are always 5 characters long
                    yr_file_to_check = f"{stn_name}h{stn_year}"
                    if stn_file[-5] != "h":
                        if yr_file_to_check in all_files:
                            files_to_delete += [file]

                        else:
                            # annual file for that station is not downloaded yet
                            print(
                                f"Warning! The annual file for this {stn_name} is not yet downloaded. Please check before proceeding."
                            )

            print(f"Number of files to delete in {rawdir}: {len(files_to_delete)}")
            print(files_to_delete)  # should have date in filename

            resp = input(
                f"\nWARNING: this will permanently delete {len(files_to_delete)} files from {rawdir}. Are you ready to delete these files? (Y/N)"
            )
            if resp == "N":
                # do not delete files!!
                print("Files not deleted. Exiting function.")

            elif resp == "Y":
                # proceed to delete files
                for file in files_to_delete:
                    # delete file -- THIS IS THE LINE OF CODE THAT MATTERS
                    s3.Object(BUCKET_NAME, file).delete()
                    print(f"File {file} deleted from AWS bucket")

            else:
                # any other character thrown in, do not delete
                print("Invalid response. Exiting function.")

    return None


if __name__ == "__main__":
    network = "MTRWFO"
    delete_files_from_AWS(network, which_to_delete="empty")
    # intentionally going to have user do this network by network, as we should be careful as to what we delete for now
