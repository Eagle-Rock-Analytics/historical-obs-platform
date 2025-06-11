"""
stnlist_update_pull.py

This function iterates through all networks and checks for missing files or stations, attempting to redownload them
and updating the station list to reflect station availability.

For station-based file systems (SCAN/SNOTEL, MADIS networks), it does the following:
1) Compare station lists to all files and identify any missing station files, attempting to redownload them.
2) MADIS only: opens the last line of each file to locate timeout errors, and redownloads from the last timestamp if error found.

For time-based station files (e.g. ASOSAWOS/OtherISD, MARITIME/NDBC), it does the following:
1) Read in station lists and file names from AWS
2) Compare station lists to all files, and identify any stations that are completely missing from downloaded data.
3) ISD-only: using start and end dates, identify any months of missing data for redownload for all stations,
and update station list to remove stations whose end date precedes the time period of analysis.

For time-based files (CIMIS), no method has been developed.
For mixed-system files (CW3E), no method has been developed.

Functions
---------
- before: Reads byte by byte backwards to return previous line from current cursor position.
- get_timeouts: Given a list of files, generate URL, read last lines and add to timeout list for re-download if error found.
- get_madis_station_timeout_csv: Timeout function search for MADIS specifically.
- madis_retry_downloads: Check for missing files or stations in MADIS networks and attempts to redownload them. 
- scan_retry_downloads: Identifies if there are any missing station files in SCAN or SNOTEL due to timeout errors and attempts to re-download them.
- maritime_retry_downloads: Identifies if there are any missing station files in MARITIME or NDBC and attempts to redownload them. 
- isd_retry_downloads: Identifies if there are any missing station files in ASOSAOWS and OtherISD due to timeout errors and attempts to re-download them
- ids_get_missing_files: Identifies missing ISD files and attempts to redownload them from ISD server into AWS bucket.
- download_comparison: Comparison of which stations downloaded, updating the station_list csv with y/n to download column.
- update_station_list: Reads in station list and updates with "Download" pass/fail flag. 
- retry_downloads: Primary function to attempt to retry timeout or missing station files per network. 

Intended Use
------------
"""

from MADIS_pull import get_madis_station_csv
from SCANSNOTEL_pull import get_scan_station_data
from ASOSAWOS_pullftp import get_asosawos_data_ftp, ftp_to_aws
from OtherISD_pull import get_otherisd_data_ftp
from MARITIME_pull import get_maritime
import boto3
import pandas as pd
import config
from smart_open import open
import os
from datetime import datetime
from ftplib import FTP
from io import StringIO
import numpy as np

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")

BUCKET_NAME  = "wecc-historical-wx"
WECC_TERR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Define lists of networks, according to type
MADIS = [
    "CAHYDRO",
    "CDEC",
    "CNRFC",
    "CRN",
    "CWOP",
    "HADS",
    "HNXWFO",
    "HOLFUY",
    "HPWREN",
    "LOXWFO",
    "MAP",
    "MTRWFO",
    "NCAWOS",
    "NOS-NWLON",
    "NOS-PORTS",
    "RAWS",
    "SGXWFO",
    "SHASAVAL",
    "VCAPCD",
]
SNTL = ["SNOTEL", "SCAN"]
ISD = ["ASOSAWOS", "OtherISD"]
MARITIME = ["MARITIME", "NDBC"]


def before(self) -> str:
    """
    Reads byte by byte backwards to return previous line from current cursor position.

    Parameters
    ----------

    Returns
    -------
    str 
        decoded string from previous line
    
    References
    ----------
    https://stackoverflow.com/questions/8721040/python-read-previous-line-and-compare-against-current-line
    """

    # get the current cursor position, then store it in memory
    current_position = self.tell()
    _tmp = b""  # Specify data in byte form
    _rea = 0  # Set parameter

    while True:
        _rea += 1
        # Move back 1 byte
        self.seek(current_position - _rea)

        # Add the one byte to the string
        _tmp += self.read(1)

        # Move back one byte again
        self.seek(current_position - _rea)

        # If the string in memory contains the "\n" more than once, we will return the string
        # alternatively break if the current position is zero (start of the string)
        if _tmp.count(b"\n") == 2 or self.tell() == 0:
            break

    # Because we're reading backwards, reverse the order and decode into a string
    return _tmp[::-1].decode()


def get_timeouts(files: list[str]) -> pd.DataFrame:
    """
    Given a list of files, generate URL, read last lines and add to timeout list for re-download if error found.

    Parameters
    ----------
    files : list[str]
        list of files to check for timeouts
    
    Returns
    -------
    ids_split : pd.DataFrame

    References
    ----------
    https://stackoverflow.com/questions/46258499/how-to-read-the-last-line-of-a-file-in-python
    """

    ids_split = []
    for file in files:
        url = f"s3://{BUCKET_NAME}/{file}"

        # Use the open method from smart_open
        with open(url, "rb") as f:  
            try:  
                # Use seek here, much faster than reading entire file into memory
                f.seek(-2, os.SEEK_END)  
                
                while f.read(1) != b"\n":
                    f.seek(-2, os.SEEK_CUR)

            except OSError:
                # catch OSError in case of a one line file
                f.seek(0)

            last_line = f.readline().decode()

            if "Timeout" in last_line:  
                # If last line is a timeout error, use f.seek to avoid reading entire file into memory
                print(f"Timeout error in {file}: processing for secondary download.")

                # Go backwards one line
                bytelen = len(last_line)
                f.seek(-bytelen, 1) 
                # Use function to seek backwards byte by byte, saving file when next '\n' reached
                prev_line = before(f)  

                # Get last completed timestamp
                time = prev_line.split(",")[1]
                # Get station ID, and add to dataframe
                station = file.split("/")[-1].replace(".csv", "")
                ids_split.append([station, time]) 

            elif ("status" in last_line):  
                # Otherwise, if other errors found in last line (not seen yet)
                print(f"Timeout error in {file}: processing for secondary download.")

                # Go backwards one line
                bytelen = len(last_line)  
                f.seek(-bytelen, 1)
                # Use function to seek backwards byte by byte, saving file when next '\n' reached
                prev_line = before(f)  

                # If station ID in previous line, treat as functional end point
                # If there are any prefixes in file name, remove
                stid = file.split("/")[-1]
                stid = stid.split("_")[-1]  
                stid = stid.replace(".csv", "")

                if stid in prev_line:
                    # Get last completed timestamp
                    time = prev_line.split(",")[1]  
                      # Get station ID and add to dataframe
                    station = file.split("/")[-1].replace(".csv", "")
                    ids_split.append([station, time])

                else:  
                    # Otherwise, just print an error to console. The last rows will be cleaned in the next stage.
                    print(f"Error: last lines of {file} not parseable.")
                    continue

    ids_split = pd.DataFrame(ids_split, columns=["STID", "start"])
    ids_split["start"] = pd.to_datetime(ids_split["start"], format="%Y-%m-%dT%H:%M:%SZ")
    ids_split["start"] = ids_split["start"].dt.strftime("%Y%m%d%H%M")

    return ids_split


def get_madis_station_timeout_csv(token: str, directory: str):
    """
    Timeout function search for MADIS specifically. Read AWS directory, 
    generate list of files, and iterate through the get_timeouts() script 
    until no more timeout errors found
    
    Parameters
    ----------
    token : str
        Synoptic API token (private)
    directory : str
        MADIS directory name

    Returns
    -------
    None
    """

    round = 1
    files = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=directory):
        file = str(item.key)
        files += [file]
    # Get list of file names
    files = list(filter(lambda f: f.endswith(".csv"), files))  

    files = [file for file in files if "errors" not in file]
    # Remove error and station list files
    files = [file for file in files if "station" not in file]
    # Only run on primary files first 
    files = [file for file in files if "2_" not in file]  

    # Get list of timeout IDs from primary downloads
    ids_split = get_timeouts(files)

    # Keep running this function as long as needed
    while ids_split.empty is False:
        round += 1
        print(f"Round {round} of downloading timeout errors")

        # Rerun pull script on timeout files from date of timeout.
        get_madis_station_csv(
            token, ids_split, directory, timeout=True, round=str(round)
        )  

        # Check to see if any of the split files needs to be split again
        # Filter by the round, checking the files just downloaded for other status errors
        files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(
            Prefix=directory + f"{str(round)}_"
        ):
            file = str(item.key)
            files += [file]

        # Get list of file names
        timeout_files = list(
            filter(lambda f: f.endswith(".csv"), files)
        ) 
        # Run on timeout files
        ids_split = get_timeouts(timeout_files)  

    print(f"No more timeout errors found in {directory} network")
    return None


def madis_retry_downloads(token: str, network: str):
    """
    Check for missing files or stations in MADIS networks and attempts to redownload them. 
    
    Parameters
    ----------
    token : str
        Synoptic API token (private)
    network : str
        name of network

    Returns
    -------
    None
    """

    # Get list of files in folder
    prefix = "1_raw_wx/" + network
    files = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=prefix):
        file = str(item.key)
        files += [file]
    files = list(filter(lambda f: f.endswith(".csv"), files))

    # ID standard station file
    station_file = [file for file in files if "stationlist" in file]  
    station_file = str(station_file[0])
    files = [file for file in files if "errors" not in file]
    # Remove error and station list files
    files = [file for file in files if "station" not in file]

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file.replace(".csv", "") for file in stations]

    # Read in station list
    station_list = s3_cl.get_object(Bucket=BUCKET_NAME, Key=station_file)
    station_list = pd.read_csv(station_list["Body"])

    # Get list of IDs not in download folder
    missed_stations = [id for id in station_list["STID"] if id not in stations]
    # Format list in way that MADIS_pull script wants it
    missed_ids = station_list[["STID", "start"]]  
    missed_ids = missed_ids[missed_ids.STID.isin(missed_stations)]

    # Check for duplicate stations in station list
    dup_stations = station_list.loc[station_list.duplicated(subset="STID", keep=False)]
    dup_stations = dup_stations[["STID", "start"]]
    # Keep first start date
    dup_todownload = dup_stations.sort_values(by="start", ascending=True).drop_duplicates(keep="first", subset=["STID"])

    # If any stations duplicated, take the earlier start data and add to the redownload queue
    download_ids = pd.merge(missed_ids, dup_todownload, how="outer")
    # If download_ids has duplicate STIDs, take first record
    download_ids = download_ids.sort_values(by="start", ascending=True).drop_duplicates( keep="first", subset=["STID"])  

    # Reorganize start date to meet API specs
    download_ids["start"] = pd.to_datetime(
        download_ids["start"], format="%Y-%m-%dT%H:%M:%SZ"
    )
    download_ids["start"] = download_ids["start"].dt.strftime("%Y%m%d%H%M")

    # Print list of stations to download
    if download_ids.empty is False:
        print("Downloading IDs:")
        print(download_ids)
        # Note here we ignore the end date of files, since we will be trimming the last 2 months of data anyways.
        # This could be changed down the road as these dates diverge.
        errors = get_madis_station_csv(
            token=token,
            directory=prefix + "/",
            ids=download_ids,
            timeout=False,
        )

        # Manually print out errors for immediate verification of success. Will also save to AWS.
        print(errors)
    else:
        print("No missing station files. Checking for timeout errors.")
    
    return None


def scan_retry_downloads(network: str, terrpath: str, marpath: str):
    """
    Identifies if there are any missing station files in SCAN or SNOTEL due to timeout errors and attempts to re-download them
    
    Parameters
    ----------
    network : str
        network name
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries

    Returns
    -------
    None
    """

    prefix = "1_raw_wx/" + network[0]
    files = []

    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=prefix):
        file = str(item.key)
        files += [file]
    # Get list of file names
    files = list(filter(lambda f: f.endswith(".csv"), files))  

    # ID station file
    station_file = [file for file in files if "station" in file]  
    station_file = str(station_file[0])
    files = [file for file in files if "errors" not in file]
     # Remove error and station list files
    files = [file for file in files if "station" not in file] 

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file.replace(".csv", "") for file in stations]

    # Read in station list
    station_list = s3_cl.get_object(Bucket=BUCKET_NAME, Key=station_file)
    station_list = pd.read_csv(station_list["Body"])

    # Get list of IDs not in download folder
    missed_ids = [id for id in station_list["stationTriplet"] if id not in stations]
    missed_stations = station_list.loc[station_list["stationTriplet"].isin(missed_ids)]

    # Print list of stations to download
    print(missed_stations)

    # Note here we ignore the end date of files, since we will be trimming the last 2 months of data anyways
    get_scan_station_data(
        terrpath,
        marpath,
        start_date=None,
        stations=missed_stations,
        primary=True,
    )

    return None


def maritime_retry_downloads(network: str) -> pd.Series:
    """
    Identifies if there are any missing station files in MARITIME or NDBC and attempts to redownload them. 

    Parameters
    ----------
    network : str
        network name
        
    Returns
    -------
    missed_stations : pd.Series
        list of staiton ids not within AWS pull folder for network

    Notes
    -----
    1. There are 2 different file formats:
        Canadian stations: STID#_csv.zip
        NDBC stations: STID#hYYYY.txt.gz
    """

    # Get list of files in folder
    prefix = "1_raw_wx/" + network
    files = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=prefix):
        file = str(item.key)
        files += [file]

    station_file = [file for file in files if "station" in file]

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file[0:5] for file in stations]  # Drop trailing metadata

    # Read in station list
    station_list = s3_cl.get_object(Bucket=BUCKET_NAME, Key=str(station_file[0]))
    station_list = pd.read_csv(station_list["Body"])

    # Get list of IDs not in download folder
    missing_ids = [id for id in station_list["STATION_ID"] if id not in stations]
    missed_stations = station_list[station_list["STATION_ID"].isin(missing_ids)]

    return missed_stations


def isd_retry_downloads(network: str) -> tuple(pd.DataFrame, pd.DataFrame):
    """
    Identifies if there are any missing station files in ASOSAOWS and OtherISD due to timeout errors and attempts to re-download them
    
    Parameters
    ----------
    network : str
        network name

    Returns
    -------
    missed_stations: pd.DataFrame
        missing stations
    missing_files: pd.DataFrame
        filepaths to missing stations
    """

    # Get list of files in folder
    prefix = "1_raw_wx/" + network
    files = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=prefix):
        file = str(item.key)
        files += [file]

    station_file = [file for file in files if "ISD" in file]
    station_file = [file for file in station_file if "station" in file]
    # Get list of file names
    files = list(filter(lambda f: f.endswith(".gz"), files))  

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file.replace(".gz", "") for file in stations]
    stations = [file[0:-5] for file in stations]

    # Read in station list
    station_list = s3_cl.get_object(Bucket=BUCKET_NAME, Key=str(station_file[0]))
    station_list = pd.read_csv(station_list["Body"])

    # Get list of IDs not in download folder
    missing_ids = [id for id in station_list["ISD-ID"] if id not in stations]
    missed_stations = station_list[station_list["ISD-ID"].isin(missing_ids)]

    # Fix formatting
    # Drop first column of dataframe
    missed_stations = missed_stations.iloc[:, 1:]
    missed_stations["WBAN"] = (
        missed_stations["WBAN"].astype(str).str.pad(5, fillchar="0")
    )

    # Get list of filenames where IDs have partially downloaded (years missing)
    downloaded_ids = station_list[~station_list["ISD-ID"].isin(missing_ids)]
    missing_files = pd.DataFrame()
    for index, id in downloaded_ids.iterrows():
        if int(id["start_time"][0:4]) < 1980:
            years = range(1980, int(id["end_time"][0:4]) + 1)

        else:
            years = range(int(id["start_time"][0:4]), int(id["end_time"][0:4]) + 1)

        id_files = [file for file in files if id["ISD-ID"] in file]
        id_success = [
            file for file in id_files if any(str(year) in file for year in years)
        ]
        missing_years = [
            year for year in years if not any(str(year) in id for id in id_success)
        ]
        if missing_years:
            filenames = [
                id["ISD-ID"] + "-" + str(missing_year) + ".gz"
                for missing_year in missing_years
            ]
            missing_files = pd.concat(
                [missing_files, pd.DataFrame(zip(missing_years, filenames))]
            )

    # Add column name
    missing_files.columns = ["year", "file_name"]

    return missed_stations, missing_files



def isd_get_missing_files(missing_files: pd.DataFrame, network: str):
    """
    Identifies missing ISD files and attempts to redownload them from ISD server into AWS bucket. 
    
    Parameters
    ----------
    missing files : pd.DataFrame
        missing ISD files (['year', 'file_name'])
    network : str
        name of network to retrieve

    Returns
    -------
    None
    """

    if missing_files.empty:
        print("No missing files to download.")
        return None
    else:
        # Set up error handling
        errors = {"Date": [], "Time": [], "Error": []}
        # Set end time to be current time at beginning of download
        end_api = datetime.now().strftime("%Y%m%d%H%M")  

        # Set up directory
        directory = "1_raw_wx/" + network + "/"

        # Login using ftplib
        ftp = FTP("ftp.ncdc.noaa.gov")
        ftp.login()  # user anonymous, password anonymous
        ftp.cwd("pub/data/noaa")
        pwd = ftp.pwd() 

        # Get list of folders (by year) in main FTP folder
        years = ftp.nlst()

        # Iterate through years
        for i in years: 
            # If folder is the name of a year (and not metadata file)
            if len(i) < 5:  
                # If folder is in our missing file list
                if int(i) in list(missing_files["year"]):  
                    try:
                        ftp.cwd(pwd)  
                        ftp.cwd(i)  
                        filenames = ftp.nlst() 

                        # Get relevant missing files
                        sub_missing = missing_files[
                            missing_files["year"] == int(i)]  
                        
                        # Get list to download
                        to_download = [
                            file
                            for file in filenames
                            if file in list(sub_missing["file_name"])
                        ] 

                        # If there are files in the list
                        if to_download:  
                            # Iterate through each one and download to AWS directory
                            for file in to_download:
                                ftp_to_aws(ftp, file, directory)

                        else:
                            print(f"No missing files available for download in folder {i}")
                            continue

                    except Exception as e:
                        print(f"Error in downloading date {i}: {e}")
                        errors["Date"].append(i)
                        errors["Time"].append(end_api)
                        errors["Error"].append(e)
                    pass
                else:
                    continue  
            else:
                continue  

        # Write errors to AWS
        csv_buffer = StringIO()
        errors = pd.DataFrame(errors)
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=directory + f"errors_{network.lower()}_{end_api}.csv",
        )

    return None


def download_comparison(network: str):
    """
    Comparison of which stations downloaded, updating the station_list csv with y/n to download column.
    In general, if a station cannot be downloaded (has a N for download) it is an ocean-observing buoy ONLY, 
    or no data is provided (optimization/testing buoy)

    Parameters
    ----------
    network : str
        network name to check

    Returns
    -------
    None
    """

    # Get list of files in folder
    directory = "1_raw_wx/" + network
    files = []

    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=directory):
        file = str(item.key)
        files += [file]

    # Get station list
    station_file = [file for file in files if "stationlist" in file]

    # Get list of file names
    files = [file for file in files if "stationlist" not in file]  

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    # Drop hYYYYtxt.gz or _csv.zip
    downloaded_stations = [file[0:5] for file in stations]  

    # Read in station list
    station_csv = s3_cl.get_object(Bucket=BUCKET_NAME, Key=str(station_file[0]))
    station_csv = pd.read_csv(station_csv["Body"])

    # Adds download column so we can compare post full data pull, will get filled after full pull
    # Mainly important for the oceanographic buoys that do not contain wx obs but are flagged as a part of WECC
    station_csv["Pulled"] = np.where(
        station_csv["STATION_ID"].isin(downloaded_stations), "Y", "N"
    )
    station_csv["Time_Checked"] = pd.to_datetime("now", utc=True).replace(microsecond=0)

    # Moves Note column to last column because of weird lat-lon column overwriting -- metadata issue for cleaning
    station_csv = station_csv.reindex(
        columns=[col for col in station_csv.columns if col != "NOTE"] + ["NOTE"]
    )

    # Reorders the indices in both stationlists
    # Previously was the full index from station_table, so there was a mismatch in index and actual number of provided stations
    station_csv.reset_index(inplace=True, drop=True)

    # Write stations to respective AWS bucket
    new_buffer = StringIO()
    station_csv.to_csv(new_buffer)
    content = new_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"/stationlist_{network}.csv",
    )

    pull_y = station_csv["Pulled"].value_counts()["Y"]
    pull_n = station_csv["Pulled"].value_counts()["N"]

    print(f"{network} station_list updated. {pull_y} successful downloads, {pull_n} failed downloads.")  
    print(station_csv.head(5))

    return None


def update_station_list(network: str) -> pd.DataFrame:
    """
    Reads in station list and updates with "Download" pass/fail flag. 
    
    Parameters
    ----------
    network : str
        network name to generate an updated station list for

    Returns
    -------
    station_csv : pd.DataFrame
        udpated stationlist with "Download" pass/fail flag
    """

    directory = "1_raw_wx/" + network + "/"
    
    # Read in station list to compare against
    files = []
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=directory):
        file = str(item.key)
        files += [file]

    # Get station file and list of files
    station_file = [file for file in files if "stationlist" in file]

    # Get list of file names
    if network in MADIS or network in SNTL:
        files = list(filter(lambda f: f.endswith(".csv"), files))  

    elif network in ISD:
        files = list(filter(lambda f: f.endswith(".gz"), files))

    if network == "ASOSAWOS":
        station_file = [
            file for file in station_file if "ISD" in file
        ]

    station_file = str(station_file[0])

    # Read in station file
    test = s3.Bucket(BUCKET_NAME).Object(station_file).get()
    station_csv = pd.read_csv(test["Body"])

    # Filter files
    files = [file for file in files if "errors" not in file]
    # Remove error and station list files
    files = [file for file in files if "station" not in file]  

    # Get list of downloaded IDs from file names
    downloaded_stns = set()
    for file in files:
        # Grab station_id from each filename
        dn_file = file.split("/")[-1].replace(".csv", "")  
        # Grabs station_id from each filename
        dn_file = dn_file.replace(".gz", "")  
        
        if network in MADIS:
            # Remove leading extension, if file is one of multiple
            dn_file = dn_file.split("_")[-1]

        elif network in ISD:
            # Remove year suffix (e.g. '-2015')
            dn_file = dn_file.rsplit("-", 1)[0]  

        downloaded_stns.add(dn_file)
    downloaded_stns = list(downloaded_stns)

    # Adds download column so we can compare post full data pull, will get filled after full pull
    # For each network, specify ID column to check against.
    if network in MADIS:
        station_csv["Pulled"] = np.where(
            station_csv["STID"].isin(downloaded_stns), "Y", "N"
        )

    elif network in SNTL:
        # Split triplet ID into 3 subcomponents
        station_csv[["id", "state", "subnetwork"]] = station_csv[
            "stationTriplet"
        ].str.split(":", expand=True)
        station_csv["Pulled"] = np.where(
            station_csv["stationTriplet"].isin(downloaded_stns), "Y", "N"
        )
    elif network in ISD:
        station_csv["Pulled"] = np.where(
            station_csv["ISD-ID"].isin(downloaded_stns), "Y", "N"
        )

    station_csv["Time_Checked"] = pd.to_datetime("now", utc=True).replace(microsecond=0)

    # Drop first column from station_csv
    station_csv = station_csv.drop(columns=station_csv.columns[0], axis=1)

    # For ISD stations, also double check that stations in list meet time criteria and filter any that aren't, updating station list
    # (This gets done in practice in the code, but isn't yet reflected in our station lists)
    if network in ISD:
        # Remove rows where data ends before 1980
        station_csv = station_csv[
            station_csv["end_time"].str[0:4].astype(int) >= 1980
        ]  

    # Write stations to respective AWS bucket
    new_buffer = StringIO()
    station_csv.to_csv(new_buffer)
    content = new_buffer.getvalue()

    s3_cl.put_object(Bucket=BUCKET_NAME, Body=content, Key=station_file)
    pull_y = station_csv["Pulled"].value_counts()["Y"]
    pull_n = station_csv["Pulled"].value_counts()["N"]

    if "N" in station_csv["Pulled"].values:
        print(f"{network} station_list updated. {pull_y} successful downloads, {pull_n} failed downloads.") 
    else:
        print(f"{network} station_list updated. {pull_y} successful downloads.")

    return station_csv


def retry_downloads(token: str, networks: list[str] | None=None):
    """
    Primary function to attempt to retry timeout or missing station files per network. 

    Parameters
    ----------
    token : str
        Synoptic API token (private)
    networks : list[str]
        name of network to retry download

    Returns
    -------
    None
    """

    # If network not provided, get list of all networks from AWS bucket and iterate
    if networks is None:
        print("All networks specified.")
        response = s3_cl.list_objects_v2(
            Bucket=BUCKET_NAME, Prefix="1_raw_wx/", Delimiter="/"
        )
        # generate list of all network prefixes from AWS bucket response
        networks = [
            prefix["Prefix"][:-1].replace("1_raw_wx/", "")
            for prefix in response["CommonPrefixes"]
        ]

    # List of networks to attempt redownload is provided
    else:
        # Remove any spaces or slashes in inputted network names
        networks = [i.replace(" ", "") for i in networks]
        networks = [i.replace("/", "") for i in networks]

    for network in networks:
        print(f"Attempting to download missing files for {network} network")
        directory = "1_raw_wx/" + network + "/"

        # MADIS
        if network in MADIS:
            # Retry any missing downloads
            madis_retry_downloads(token, network)  

            # Get timeout CSVs
            print(f"Identifying file timeouts for {network} network")
            get_madis_station_timeout_csv(token=token,  directory=directory)
            update_station_list(network)

        # SNOTEL/SCAN
        elif network in SNTL:
            # Get IDs
            scan_retry_downloads(
                network=[network],
                terrpath=WECC_TERR,
                marpath=WECC_MAR,
            )  
            update_station_list( network)

        # ASOSAWOS / OtherISD
        elif network in ISD:
            missed_stations, file_list = isd_retry_downloads(network=network)

            if len(missed_stations) == 0:  # If no missing stations
                print("All selected stations have data downloaded.")
            else:
                print(f"Attempting to download data for {len(missed_stations)} missing stations.")
                if network == "ASOSAWOS":
                    # Download all missing stations
                    get_asosawos_data_ftp(
                        missed_stations, directory, get_all=True
                    )
    
                elif network == "OtherISD":
                    # Download all missing stations
                    get_otherisd_data_ftp(
                        missed_stations,  directory, get_all=True
                    )

            # Regenerate file_list after full station download
            missed_stations, file_list = isd_retry_downloads(network=network)

            # Download missing files
            print(f"Attempting to download {len(file_list)} missing files.")
            isd_get_missing_files(file_list,  network)

            # update station list
            update_station_list(network)

        # MARITIME / NDBC
        elif network in MARITIME:
            # Get list of missing stations, missing files
            missed_stations = maritime_retry_downloads(network)
            get_maritime(missed_stations,  network, years=None)
            download_comparison( network)

        else:
            print(f"{network} network not currently configured for download retry.")
            continue

    return None


if __name__ == "__main__":
    retry_downloads(token=config.token, networks=["CWOP"])
    # If networks not specified, will attempt all networks (generating list from folders in raw bucket.)
