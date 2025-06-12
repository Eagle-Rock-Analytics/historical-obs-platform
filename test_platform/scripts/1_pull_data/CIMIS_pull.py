"""
CIMIS_pull.py

This script downloads CIMIS data from the CADWR using ftp.
Approach:
(1) Get station list (does not need to be re-run constantly)
(2) Download data using station list.

Functions
---------
- get_cimis_stations: Get up to date station list of CIMIS stations and downloads to AWS S3 the Stations List file from the FTP server
- get_cimis_data_ftp: Query FTP server for CIMIS data and download raw csv files, use for the full retrieval
- get_cimis_update_ftp: Query ftp server for CIMIS data and download csv files. Use to update data.

Intended Use
------------ 
Retrieves raw data for an individual network, all variables, all times. Organized by station, with 1 file per year.

Notes
-----
1. The file for each station-year is updated daily for the current year.
To pull real-time data, we may want to write just an API call with date ranges and stations and update the most recent year folder only.
This is a separate function/branch.
2. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

from ftplib import FTP
from datetime import datetime, timezone, date
import pandas as pd
import boto3
from io import BytesIO, StringIO

try:
    from calc_pull import ftp_to_aws
except RuntimeError as e:
    print(f"Error importing calc_pull: {e}")


s3 = boto3.client("s3")
BUCKET_NAME = "wecc-historical-wx"
DIRECTORY = "1_raw_wx/CIMIS/"


def get_cimis_stations(directory: str):
    """
    Get up to date station list of CIMIS stations and downloads to AWS S3 the Stations List file from the FTP server.
    All stations in CA, so no need to filter out stations here. File just downloaded for reference.

    Parameters
    ----------
    directory : str
        path to AWS bucket

    Returns
    -------
    None
    """

    # Login using ftplib, get list of stations as csv
    ftp = FTP("ftpcimis.water.ca.gov")
    ftp.login()  # user anonymous, password anonymous
    ftp.cwd("pub2")

    # Get station list
    ftp_to_aws(
        ftp,
        "CIMIS Stations List (January20).xlsx",
        directory,
        rename="stationlist_CIMIS.xlsx",
    )

    # Get recent units
    ftp_to_aws(
        ftp,
        "readme-ftp-Revised5units.txt",
        "1_raw_wx/CIMIS/",
        rename="ftpunits_post2014.txt",
    )

    # Get older units
    ftp_to_aws(
        ftp,
        "readme (prior to June 2014)units.txt",
        "1_raw_wx/CIMIS/",
        rename="ftpunits_pre2014.txt",
    )

    return None


def get_cimis_data_ftp(directory: str, years: list[str] | None, get_all: bool = True):
    """
    Query FTP server for CIMIS data and download raw csv files, use for the full retrieval.

    Parameters
    ----------
    directory : str
        name of folder to save raw data into
    years : list[str], optional
        specify years to download data from as a range of years mapped to strings (e.g. ['1998', '1999'])
    get_all : bool
        If False, only download files whose last edit date is newer than the most recent files downloaded in the save folder.
        Only use to update a complete set of files.

    Returns
    -------
    None
    """

    # Set up error handling
    errors = {"File": [], "Time": [], "Error": []}
    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    # Login using ftplib
    ftp = FTP("ftpcimis.water.ca.gov")
    ftp.login()  # user anonymous, password anonymous
    ftp.cwd("pub2/")
    pwd = ftp.pwd()

    try:
        objects = s3.list_objects(Bucket=BUCKET_NAME, Prefix=directory)
        all = objects["Contents"]

        # Get date of last edited file
        latest = max(all, key=lambda x: x["LastModified"])
        last_edit_time = latest["LastModified"]

        # Get list of all file names
        alreadysaved = []
        for item in all:
            files = item["Key"]
            alreadysaved.append(files)
        alreadysaved = [ele.replace(directory, "") for ele in alreadysaved]

    except:
        # If folder empty or there's an error with the "last downloaded" metadata, re-download all data
        get_all = True

    try:
        # First all data not from current year through annual folder
        ftp.cwd("annualMetric/")

        # Get list of all file names in folder
        filenames = ftp.nlst()
        filenames = [i for i in filenames if i.endswith(".zip")]
        # Only keep hourly data, drop daily files
        filenames = [i for i in filenames if i.startswith("hourly")]

        if years is not None:
            # filter files by year
            filenames = [str for str in filenames if any(sub in str for sub in years)]

        for filename in filenames:
            # Returns time modified (in UTC)
            modifiedTime = ftp.sendcmd("MDTM " + filename)[4:].strip()
            # Convert to datetime
            modifiedTime = datetime.strptime(modifiedTime, "%Y%m%d%H%M%S").replace(
                tzinfo=timezone.utc
            )

            # If get_all is False, only download files whose last edit date has changed since the last download or whose filename is not in the folder
            if get_all is False:
                # If filename already in bucket
                if filename in alreadysaved:
                    # If file new since last run-through, write to folder
                    if modifiedTime > last_edit_time:
                        ftp_to_aws(ftp, filename, directory)
                    else:
                        print(f"{filename} already saved")
                else:
                    # Else, if filename not saved already, save
                    ftp_to_aws(ftp, filename, directory)

            # If get_all is true, download all files in folder
            elif get_all is True:
                ftp_to_aws(ftp, filename, directory)

        # Now, repeat to download present year's data (housed in the 'hourly' folder)
        pres_year = date.today().strftime("%Y")

        # If years includes present year or years not specified.
        if (years is None) or (pres_year in years):
            ftp.cwd(pwd)
            ftp.cwd("monthlyMetric/")

            # Get list of all file names in folder.
            filenames = ftp.nlst()
            filenames = [i for i in filenames if i.endswith(".zip")]
            # Only keep hourly data, drop daily files.
            filenames = [i for i in filenames if i.startswith("hourly")]

            for filename in filenames:
                if get_all is True:
                    # We don't save hourlyallstns.zip because it'll be less compatible with the get_all = False option or future updating runs. Instead, download all files.
                    ftp_to_aws(ftp, filename, directory)
                else:
                    # Returns time modified (in UTC)
                    modifiedTime = ftp.sendcmd("MDTM " + filename)[4:].strip()
                    # Convert to datetime
                    modifiedTime = datetime.strptime(
                        modifiedTime, "%Y%m%d%H%M%S"
                    ).replace(tzinfo=timezone.utc)

                    # If get_all is False, only download files whose last edit date has changed since the last download or whose filename is not in the folder.
                    if get_all is False:
                        # If filename already in saved bucket
                        if filename in alreadysaved:
                            # If file new since last run-through, write to folder
                            if modifiedTime > last_edit_time:
                                ftp_to_aws(ftp, filename, directory)
                            else:
                                print(f"{filename} already saved")
                        else:
                            # Else, if filename not saved already, save.
                            ftp_to_aws(ftp, filename, directory)

    except Exception as e:
        print(f"Error in downloading file {filename}: {e}")
        errors["File"].append(filename)
        errors["Time"].append(end_api)
        errors["Error"].append(e)

    # close connection
    ftp.quit()

    # Write errors to csv
    csv_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"errors_cimis_{end_api}.csv",
    )

    return None


def get_cimis_update_ftp(directory: str, start_date: str | None, end_date: str | None):
    """
    Query ftp server for CIMIS data and download csv files. Use to update data.

    Parameters
    ----------
    directory : str
        folder within bucket
    start_date : str
        start year to search, optional
    end_date : str
        end year to search, optional

    Returns
    -------
    None
    """
    # Set up error handling
    errors = {"File": [], "Time": [], "Error": []}
    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    # Login using ftplib
    ftp = FTP("ftpcimis.water.ca.gov")
    ftp.login()  # user anonymous, password anonymous
    ftp.cwd("pub2/")
    pwd = ftp.pwd()

    try:
        # First all data not from current year through annual folder
        ftp.cwd("annualMetric/")

        # Get list of all file names in folder
        filenames = ftp.nlst()
        filenames = [i for i in filenames if i.endswith(".zip")]
        # Only keep hourly data, drop daily files
        filenames = [i for i in filenames if i.startswith("hourly")]

        # Get present year
        pres_year = date.today().strftime("%Y")

        # If years includes previous years
        if pres_year != (start_date[0:4]):
            if start_date is not None:
                # Only keep filenames after start year
                filenames = [
                    str for str in filenames if int(str[-8:-4]) >= int(start_date[0:4])
                ]

            if end_date is not None:
                # Only keep filenames after start year
                filenames = [
                    str for str in filenames if int(str[-8:-4]) <= int(end_date[0:4])
                ]

            for filename in filenames:
                ftp_to_aws(ftp, filename, directory)

        if (end_date is None) or (pres_year == end_date[0:4]):
            # If years includes present year or years not specified
            ftp.cwd(pwd)
            ftp.cwd("monthlyMetric/")

            # Get list of all file names in folder
            filenames = ftp.nlst()
            filenames = [i for i in filenames if i.endswith(".zip")]
            # Only keep hourly data, drop daily files
            filenames = [i for i in filenames if i.startswith("hourly")]

            # Filter months by start and end dates
            if start_date is not None and start_date[0:4] == pres_year:
                if end_date is None:
                    # Get all months from start date to present
                    dates = (
                        pd.date_range(
                            start_date[0:7], date.today(), freq="MS", inclusive="both"
                        )
                        .strftime("%b")
                        .tolist()
                    )
                else:
                    # Get all months from start to end date
                    dates = (
                        pd.date_range(
                            start_date[0:7], end_date[0:7], freq="MS", inclusive="both"
                        )
                        .strftime("%b")
                        .tolist()
                    )
            elif end_date is not None:
                # Get all months from January up through end date
                dates = (
                    pd.date_range(
                        datetime(year=date.today().year, month=1, day=1),
                        end_date[0:7],
                        freq="MS",
                        inclusive="both",
                    )
                    .strftime("%b")
                    .tolist()
                )
            elif end_date is None:  # Get all months available
                dates = (
                    pd.date_range(
                        datetime(year=date.today().year, month=1, day=1),
                        date.today(),
                        freq="MS",
                        inclusive="both",
                    )
                    .strftime("%b")
                    .tolist()
                )

            filenames = [
                file
                for file in filenames
                if any(date.lower() in file for date in dates)
            ]

            for filename in filenames:
                ftp_to_aws(ftp, filename, directory)

    except Exception as e:
        print(f"Error in downloading file {filename}: {e}")
        errors["File"].append(filename)
        errors["Time"].append(end_api)
        errors["Error"].append(e)

    # close connection
    ftp.quit()

    # Write errors to csv
    csv_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"errors_cimis_{end_api}.csv",
    )

    return None


if __name__ == "__main__":
    # Run functions
    get_cimis_stations(DIRECTORY)
    get_cimis_data_ftp(DIRECTORY, years=None, get_all=True)

# Note, for first full data pull, set get_all = True
# For all subsequent data pulls/update with newer data, set get_all = False
