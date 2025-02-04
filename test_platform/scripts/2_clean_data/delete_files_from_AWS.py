"""This script deletes files in batch from AWS. Options for which files to delete are provided. """

import pandas as pd
import s3fs
import boto3
from cleaning_helpers import get_file_paths

s3 = boto3.resource("s3")


def delete_files_from_AWS(network, which_to_delete):
    """
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
        this function does not return a value

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

    which_to_delete = "no_wx_data" ?
        Could be set-up to delete raw files that do not contain any meteorological info, i.e., buoys
        Not continuing for now.

    If which_to_delete is not provided an option, nothing will delete.

    References
    ----------
    [1] Modified from: https://stackoverflow.com/questions/32501995/boto3-s3-renaming-an-object-using-copy-object
    """

    bucket_name = "wecc-historical-wx"

    # deletes only update pull files with date extensions that are empty
    if which_to_delete == "empty":
        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        # get files
        all_files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix=rawdir):
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
                fn_date = file.split("_")[
                    -1
                ]  # grabs part of filename after _ character
                if "-" in fn_date:  # must follow YYYY-MM-DD format
                    # print(file)
                    files_to_delete += [file]

        # open each file, and assess if empty
        empty = []
        for file in files_to_delete:
            print("Checking {} for empty size".format(file))
            fs = s3fs.S3FileSystem()
            aws_url = "s3://" + bucket_name + "/" + file

            with fs.open(aws_url) as fileObj:
                df = pd.read_csv(fileObj, header=6, low_memory=False)

                # assess if empty
                if (
                    len(df.index) == 1 and df["Station_ID"].isnull().values.any()
                ):  # any station that is missing the station_id is empty
                    empty += [file]

        print("Number of empty files to delete in {0}: {1}".format(rawdir, len(empty)))
        print(empty)

        resp = input(
            "\nWARNING: this will permanently delete {0} files from {1}. Are you ready to delete these files? (Y/N)".format(
                len(empty), rawdir
            )
        )
        if resp == "N":  # do not delete files!!
            print("Files not deleted. Exiting function.")

        elif resp == "Y":  # proceed to delete files
            for file in empty:
                s3.Object(
                    bucket_name, file
                ).delete()  # delete file -- THIS IS THE LINE OF CODE THAT MATTERS
                print("File {0} deleted from AWS bucket".format(file))

        else:  # any other character thrown in, do not delete
            print("Invalid response. Exiting function.")

    # deletes ALL update pull files with date extensions
    elif which_to_delete == "update_pull":
        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        # get files
        all_files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix=rawdir):
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
                fn_date = file.split("_")[
                    -1
                ]  # grabs part of filename after _ character
                if "-" in fn_date:  # must follow YYYY-MM-DD format
                    files_to_delete += [file]

        print(
            "Number of files to delete in {0}: {1}".format(rawdir, len(files_to_delete))
        )
        print(files_to_delete)  # should have date in filename

        resp = input(
            "\nWARNING: this will permanently delete {0} files from {1}. Are you ready to delete these files? (Y/N)".format(
                len(files_to_delete), rawdir
            )
        )
        if resp == "N":  # do not delete files!!
            print("Files not deleted. Exiting function.")

        elif resp == "Y":  # proceed to delete files
            for file in files_to_delete:
                s3.Object(
                    bucket_name, file
                ).delete()  # delete file -- THIS IS THE LINE OF CODE THAT MATTERS
                print("File {0} deleted from AWS bucket".format(file))

        else:  # any other character thrown in, do not delete
            print("Invalid response. Exiting function.")

    # deletes the monthly files from NDBC, MARITIME after year file is in
    # note: CIMIS has monthly files, but should not delete, as these files retain data for all stations
    elif which_to_delete == "monthly":
        # grab raw data directory for specific network
        rawdir, cleandir, qaqcdir = get_file_paths(network)

        if network == "NDBC" or network == "MARITIME":
            # get files
            all_files = []
            for item in s3.Bucket(bucket_name).objects.filter(Prefix=rawdir):
                file = str(item.key)
                all_files += [file]

            all_files = [file for file in all_files if "stationlist" not in file]
            all_files = [file for file in all_files if "error" not in file]

            # identify monthly files
            files_to_delete = []
            for file in all_files:
                # print(file)
                file_parts = file.split("/")
                stn_file = file_parts[-1]  # grabs just file name, no path
                stn_file = stn_file.split(".")[
                    0
                ]  # removes file extenstion (.txt.gz, or .zip)

                if "csv" in stn_file:
                    continue  # skip canadian files, they're annual

                else:
                    stn_year = stn_file[
                        -4:
                    ]  # use to identify if annual file is already downloaded
                    stn_name = stn_file[
                        :5
                    ]  # NDBC and MARITIME station names are always 5 characters long

                    yr_file_to_check = stn_name + "h" + stn_year
                    if stn_file[-5] != "h":
                        if yr_file_to_check in all_files:
                            files_to_delete += [file]

                        else:  # annual file for that station is not downloaded yet
                            print(
                                "Warning! The annual file for this {0} is not yet downloaded. Please check before proceeding.".format(
                                    stn_name
                                )
                            )

            print(
                "Number of files to delete in {0}: {1}".format(
                    rawdir, len(files_to_delete)
                )
            )
            print(files_to_delete)  # should have date in filename

            resp = input(
                "\nWARNING: this will permanently delete {0} files from {1}. Are you ready to delete these files? (Y/N)".format(
                    len(files_to_delete), rawdir
                )
            )
            if resp == "N":  # do not delete files!!
                print("Files not deleted. Exiting function.")

            elif resp == "Y":  # proceed to delete files
                for file in files_to_delete:
                    s3.Object(
                        bucket_name, file
                    ).delete()  # delete file -- THIS IS THE LINE OF CODE THAT MATTERS
                    print("File {0} deleted from AWS bucket".format(file))

            else:  # any other character thrown in, do not delete
                print("Invalid response. Exiting function.")

    return None


if __name__ == "__main__":
    network = "MTRWFO"
    delete_files_from_AWS(network, which_to_delete="empty")
    # intentionally going to have user do this network by network, as we should be careful as to what we delete for now
