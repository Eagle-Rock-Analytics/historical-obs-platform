"""
OtherISD_pull.py

This script downloads all non-ASOS/AWOS data from ISD using ftp.
Approach:
(1) Download ISD station list and get ASOSAWOS station list from AWS.
(2) Download data using station list.

Functions
---------
- get_wecc_stations: Retrieves to get up to date station list of ISD stations in WECC, and remove all asos-awos stations
- get_otherisd_data_ftp: Query ftp server for non-ASOS/AWOS ISD data and download zipped files.

Intended Use
-------------
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
from datetime import datetime, timezone
import pandas as pd
from shapely.geometry import Point
import pandas as pd
import geopandas as gp
from geopandas.tools import sjoin
import boto3  # For AWS integration.
from io import BytesIO, StringIO

from calc_pull import ftp_to_aws, get_wecc_poly

s3 = boto3.client("s3")
BUCKET_NAME = "wecc-historical-wx"
DIRECTORY = "1_raw_wx/OtherISD/"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


def get_wecc_stations(terrpath: str, marpath: str, directory: str) -> pd.DataFrame:
    """
    Retrieves to get up to date station list of ISD stations in WECC, and remove all asos-awos stations

    Parameters
    ----------
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    directory : str
        AWS path name

    Returns
    -------
    weccstations : pd.DataFrame
        stations within WECC
    """

    # Login using ftplib, get list of stations as csv
    filename = "isd-history.csv"
    ftp = FTP("ftp.ncdc.noaa.gov")
    ftp.login()  # user anonymous, password anonymous
    ftp.cwd("pub/data/noaa/")

    # Read in ISD stations
    r = BytesIO()
    ftp.retrbinary("RETR " + filename, r.write)
    r.seek(0)

    # Read in csv and only filter to include US stations
    stations = pd.read_csv(r)
    weccstations = stations[(stations["CTRY"] == "US")]

    # Use spatial geometry to only keep points in wecc marine / terrestrial areas
    geometry = [Point(xy) for xy in zip(weccstations["LON"], weccstations["LAT"])]
    weccgeo = gp.GeoDataFrame(weccstations, crs="EPSG:4326", geometry=geometry)
    # get bbox of WECC to use to filter stations against
    t, m, bbox = get_wecc_poly(terrpath, marpath)

    # Get terrestrial stations
    weccgeo = weccgeo.to_crs(t.crs)
    # Only keep stations in terrestrial WECC region
    terwecc = sjoin(weccgeo.dropna(), t, how="left")
    terwecc = terwecc.dropna()

    # Get marine stations
    marwecc = sjoin(weccgeo.dropna(), m, how="left")
    marwecc = marwecc.dropna()

    # Join and remove duplicates using USAF and WBAN as combined unique identifier
    weccstations = pd.concat(
        [terwecc.iloc[:, :11], marwecc.iloc[:, :11]], ignore_index=True, sort=False
    ).drop_duplicates(["USAF", "WBAN"], keep="first")

    # Generate ID from USAF/WBAN combo for API call. This follows the naming convention used by FTP/AWS for file names
    # Add leading zeros where they are missing from WBAN stations
    weccstations["ISD-ID"] = (
        weccstations["USAF"]
        + "-"
        + weccstations["WBAN"].astype("str").str.pad(5, side="left", fillchar="0")
    )

    # Reformat time strings for FTP/API call
    weccstations["start_time"] = [
        datetime.strptime(str(i), "%Y%m%d").strftime("%Y-%m-%d")
        for i in weccstations["BEGIN"]
    ]
    weccstations["end_time"] = [
        datetime.strptime(str(i), "%Y%m%d").strftime("%Y-%m-%d")
        for i in weccstations["END"]
    ]

    # Read in ASOS and AWOS station files and use to filter to remove ASOS/AWOS stations
    # Note, this relies on having run the ASOSAWOS pull script prior
    asosawos = pd.read_csv(f"s3://{BUCKET_NAME}/1_raw_wx/ASOSAWOS/stationlist_ASOSAWOS.csv")

    # create mask and filter
    m1 = weccstations.WBAN.isin(asosawos.WBAN)
    weccstations = weccstations[~m1]
    weccstations.reset_index(inplace=True, drop=True)

    # Write non-ASOS AWOS station list to CSV
    csv_buffer = StringIO()
    weccstations.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3.put_object(
        Bucket=BUCKET_NAME, Body=content, Key=directory + "stationlist_OtherISD.csv"
    )

    return weccstations


def get_otherisd_data_ftp(
    station_list: pd.DataFrame,
    directory: str,
    start_date: str | None = None,
    get_all: bool = True,
):
    """
    Query ftp server for non-ASOS/AWOS ISD data and download zipped files.

    Parameters
    ----------
    station_list : pd.DataFrame
        stations to retrieve
    directory : str
        AWS folder to save to
    start_date : str, optional
        subset date start, format "YYYY-MM-DD"
    get_all : bool, optional
        If False, only download files whose last edit date is newer than
        the most recent files downloaded in the save folder. Only use to update a complete set of files.

    Returns
    -------
    None
    """
    # Set up error handling
    errors = {"Date": [], "Time": [], "Error": []}
    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    # Login using ftplib
    ftp = FTP("ftp.ncdc.noaa.gov")
    ftp.login()  # user anonymous, password anonymous
    ftp.cwd("pub/data/noaa")
    pwd = ftp.pwd()

    # Get list of folders (by year) in main FTP folder
    years = ftp.nlst()

    # If no start date specified, manually set to be Jan 01 1980
    if start_date is None:
        start_date = "1980-01-01"

    # Remove depracated stations if filtering by time
    if start_date is not None:
        try:
            # Filter to ensure station is not depracated before time period of interest
            station_list = station_list[station_list["end_time"] >= start_date]
        except Exception as e:
            print(f"Error:", {e})
            # function will use years to filter station files
            years = [i for i in years if (len(i) < 5 and int(i) > 1979)]

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
        # If folder empty or there's an error with the "last downloaded" metadata, redownload all data
        get_all = True

    for i in years:
        # If folder is the name of a year (and not metadata file)
        if len(i) < 5:
            if (
                start_date is not None and int(i) >= int(start_date[0:4])
            ) or start_date is None:
                # If no start date specified or year of folder is within start date range, download folder
                try:
                    ftp.cwd(pwd)
                    ftp.cwd(i)
                    # Get list of all file names in folder
                    filenames = ftp.nlst()
                    # Reformat station IDs to match file names
                    filefiltlist = station_list["ISD-ID"] + "-" + i + ".gz"
                    filefiltlist = filefiltlist.tolist()

                    # Only pull all file names that are contained in station_list ID column
                    fileswecc = [x for x in filenames if x in filefiltlist]

                    for filename in fileswecc:
                        # Returns time modified (in UTC)
                        modifiedTime = ftp.sendcmd("MDTM " + filename)[4:].strip()
                        # Convert to datetime
                        modifiedTime = datetime.strptime(
                            modifiedTime, "%Y%m%d%H%M%S"
                        ).replace(tzinfo=timezone.utc)

                        # If get_all is False, only download files whose last edit date has changed since the last download or whose filename is not in the folder
                        if get_all is False:
                            if filename in alreadysaved:
                                # If filename already in saved bucket
                                if modifiedTime > last_edit_time:
                                    # If file new since last run-through, write to folder
                                    ftp_to_aws(ftp, filename, directory)
                                else:
                                    print(f"{filename} already saved")
                            else:
                                # Else, if filename not saved already, save
                                ftp_to_aws(ftp, filename, directory)

                        elif get_all is True:
                            # If get_all is true, download all files in folder
                            ftp_to_aws(ftp, filename, directory)

                except Exception as e:
                    print(f"Error in downloading date {i}: {e}")
                    errors["Date"].append(i)
                    errors["Time"].append(end_api)
                    errors["Error"].append(e)
                    continue

            else:
                # If year of folder not in start date range, skip folder
                continue

        else:
            # Skip if file or folder isn't a year. Can change to print file/folder name, or to save other metadata files as desired
            continue

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
        Key=directory + f"errors_otherisd_{end_api}.csv",
    )

    return None


if __name__ == "__main__":
    # Run functions
    stations = get_wecc_stations(WECC_TERR, WECC_MAR, DIRECTORY)
    get_otherisd_data_ftp(stations, DIRECTORY, start_date=None, get_all=True)
