"""
update_pull.py

This script contains the functions needed to retrieve updated, or "newer" data per network, following a "full pull" of data. 
Data will be downloaded from the last date updated in AWS until 45 days prior to the present day (the window of final data).

Functions
---------
- get_last_date: Get last download date of files in AWS folder.
- update_asosawos: Retrieves newer data for ASOSAWOS since the date of the last data retrieval
- update_cimis: Retrieves newer data for CIMIS since the date of the last data retrieval
- update_cw3e: Retrieves newer data for CW3E since the date of the last data retrieval
- update_hads: Retrieves newer data for HADS since the date of the last data retrieval
- update_madis: Retrieves newer data for MADIS networks since the date of the last data retrieval
- update_maritime: Retrieves newer data for MARITIME / NDBC since the date of the last data retrieval
- update_otherisd: Retrieves newer data for OtherISD since the date of the last data retrieval
- update_SCAN: Retrieves newer data for SCAN / SNOTEL since the date of the last data retrieval

Intended Use
------------
This script will be run to retrieve "newer" data since the date of the last data retrieval. 
Example: 
>> A "Full Pull" for HADS occurred Jan 15, 2024. 
>> A "Update Pull" using this script would retrieve all data from Jan 15, 2024 up to 45 days prior to present day. 

Notes
-----
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

from datetime import datetime, timezone, timedelta
import ASOSAWOS_pullftp
from SCANSNOTEL_pull import get_scan_station_data
from OtherISD_pull import get_wecc_stations, get_otherisd_data_ftp
from MADIS_pull import get_madis_metadata, get_madis_station_csv
from CW3E_pull import get_cw3e_metadata, get_cw3e_update
from CIMIS_pull import get_cimis_update_ftp
from HADS_pull import get_hads_update
from MARITIME_pull import get_maritime_update, get_maritime_station_ids
from MADIS_pull import madis_update
from stnlist_update_pull import retry_downloads
import boto3
import config

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")
BUCKET_NAME = "wecc-historical-wx"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

today = datetime.now(timezone.utc).date()
download_window = 45  # Preliminary data lag
download_date = today - timedelta(days=download_window)  # Set to midnight UTC time


def get_last_date(folder: str, file_ext: str | None = None) -> str:
    """Get last download date of files in AWS folder.

    Parameters
    ----------
    folder : str
        name of network directory
    file_ext: str, optional
        file extension to look for (e.g. .gz, only useful if not .csv or multiple file types)

    Returns
    -------
    datetime.datetime.date
        last modified date
    """

    files = s3.Bucket(BUCKET_NAME).objects.filter(Prefix=folder)
    if file_ext:
        files = [
            [obj.key, obj.last_modified]
            for obj in sorted(files, key=lambda x: x.last_modified, reverse=True)
            if obj.key.endswith(file_ext)
            and not any(x in obj.key for x in ("errors", "station"))
        ]
    else:
        files = [
            [obj.key, obj.last_modified]
            for obj in sorted(files, key=lambda x: x.last_modified, reverse=True)
            if not any(x in obj.key for x in ("errors", "station"))
        ]

    # Get most recent and nth most recent time
    last_time_modified_list = [x[1] for x in files]
    # most recent date of data pull from stations that require an updated data

    last_time_modified = max(last_time_modified_list)
    return last_time_modified.date()


def update_asosawos(last_time_mod: str | None = None):
    """Retrieves newer data for ASOSAWOS since the date of the last data retrieval.
    As currently written, this will overwrite all files for the current year.

    Parameters
    ----------
    last_time_mod : str, optional
        time of last modification

    Returns
    -------
    None
    """

    network = "ASOSAWOS"
    directory = f"1_raw_wx/{network}/"

    if last_time_mod is None:
        last_time_mod = get_last_date(folder=directory, file_ext=".gz")

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = ASOSAWOS_pullftp.get_wecc_stations(WECC_TERR, WECC_MAR)
        ASOSAWOS_pullftp.get_asosawos_data_ftp(
            stations, directory, start_date=str(last_time_mod)
        )
        retry_downloads(token=config.token, networks=[network])

    else:
        print(f"{network} station files up to date.")

    return None


def update_cimis(last_time_mod: str | None = None):
    """
    Retrieves newer data for CIMIS since the date of the last data retrieval.
    No retry download method avaialble. May overwrite most recent if pull is
    repeated more frequently than monthly.

    Parameters
    ----------
    last_time_mod : str, optional
        string of last modified time

    Returns
    -------
    None
    """

    network = "CIMIS"
    directory = f"1_raw_wx/{network}/"

    if last_time_mod is None:
        last_time_mod = get_last_date(folder=directory, file_ext=".zip")

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_cimis_update_ftp(
            directory,
            start_date=str(last_time_mod),
            end_date=str(download_date),
        )
    else:
        print(f"{network} station files up to date.")

    return None


def update_cw3e(last_time_mod: str | None = None):
    """
    Retrieves newer data for CW3E since the date of the last data retrieval.
    Will download multiple byte files for each station for each day selected, resulting
    in a long run-time.

    Parameters
    ----------
    last_time_mod : str, optional
        string of last modified time

    Returns
    -------
    None

    Notes
    -----
    1. The LBH station will alwyas be completely redownloaded. At present LBH
    gets dropped during the cleaning phase, so should not affect continued processing.
    """

    network = "CW3E"
    directory = f"1_raw_wx/{network}/"

    if last_time_mod is None:
        last_time_mod = get_last_date(folder=directory, file_ext="m")

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_cw3e_metadata(
            token=config.token,
            terrpath=WECC_TERR,
            marpath=WECC_MAR,
            directory=directory,
        )
        get_cw3e_update(
            directory,
            start_date=str(last_time_mod),
            end_date=str(download_date),
        )
    else:
        print(f"{network} station files up to date.")

    return None


def update_hads(last_time_mod: str | None = None):
    """
    Retrieves newer data for HADS since the date of the last data retrieval.

    Parameters
    ----------
    last_time_mod : str, optional
        string of last modified time

    Returns
    -------
    None
    """

    network = "HADS"
    directory = f"1_raw_wx/{network}/"

    if last_time_mod is None:
        last_time_mod = get_last_date(folder=directory, file_ext=".gz")

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        get_hads_update(
            directory,
            start_date=str(last_time_mod),
            end_date=str(download_date),
        )

    else:
        print(f"{network} station files up to date.")

    return None


def update_madis(network: str, last_time_mod: str | None = None):
    """
    Retrieves newer data for MADIS networks since the date of the last data retrieval.

    Parameters
    ----------
    network : str
        name of MADIS network to retrieve
    last_time_mod : str, optional
        time of last modification

    Returns
    -------
    None
    """

    directory = f"1_raw_wx/{network}/"

    # CAHYDRO has special handling due to space character
    if network == "CAHYDRO":
        if last_time_mod is None:
            last_time_mod = get_last_date(folder=directory, file_ext=".csv")

        if last_time_mod < download_date:
            print(
                f"Downloading {network} data from {last_time_mod} to {download_date}."
            )
            madis_update(
                token=config.token,
                networks=["CA HYDRO"],
                pause=None,
                start_date=str(last_time_mod),
                end_date=str(download_date),
            )
        else:
            print(f"{network} station files up to date.")

    # all other MADIS networks
    else:
        if last_time_mod is None:
            last_time_mod = get_last_date(
                folder=directory, n=int(n[network]), file_ext=".csv"
            )

        if last_time_mod < download_date:
            print(
                f"Downloading {network} data from {last_time_mod} to {download_date}."
            )
            madis_update(
                token=config.token,
                networks=[network],
                pause=None,
                start_date=str(last_time_mod),
                end_date=str(download_date),
            )
        else:
            print(f"{network} station files up to date.")

    return None


def update_maritime(network: str, last_time_mod: str | None = None):
    """
    Retrieves newer data for MARITIME / NDBC since the date of the last data retrieval.

    Parameters
    ----------
    network : str
        name of network to download, MARITIME or NDBC
    last_time_mod : str, optional
        time of last modification

    Returns
    -------
    None

    Notes
    -----
    1. Update delay by (minimum of) 45 days. This means 2022 data isn't available in yearly format until ~ Feb 15 2022, e.g.
    2. NDBC archive has a different file extenstion for preliminary data that is not updated until the full month has passed the preliminary period.
    Update/automation will therefore need to restrict to only pulling complete months, which may be beyond 45 days
    """

    directory = f"1_raw_wx/{network}/"

    if last_time_mod is None:
        last_time_mod = get_last_date(folder=directory)
        # Note: File extension is either .gz or .zip, so do not filter by file extension here as errors & station files will be dropped automatically.

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = get_maritime_station_ids(
            WECC_TERR, WECC_MAR, "1_raw_wx/MARITIME/", "1_raw_wx/NDBC/"
        )
        get_maritime_update(
            stations,
            network,
            start_date=str(last_time_mod),
            end_date=str(download_date),
        )
    else:
        print(f"{network} station files up to date.")

    return None


def update_otherisd(last_time_mod: str | None = None):
    """
    Retrieves newer data for OtherISD since the date of the last data retrieval.
    As currently written, this will overwrite all files for the current year.

    Parameters
    ----------
    last_time_mod : str, optional
        time of last modification

    Returns
    -------
    None

    Notes
    -----
    Timeout may occasionaly occur on API end, but resolves upon another attempt after waiting a short amount of time (~hours)
    """

    network = "OtherISD"
    directory = f"1_raw_wx/{network}/"

    if last_time_mod is None:
        last_time_mod = get_last_date(folder=directory, file_ext=".gz")

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")
        stations = get_wecc_stations(WECC_TERR, WECC_MAR, directory)
        get_otherisd_data_ftp(
            stations,
            directory,
            start_date=str(last_time_mod),
            get_all=True,
        )
        retry_downloads(token=config.token, networks=[network])

    else:
        print(f"{network} station files up to date.")

    return None


# Update script: SCAN
def update_SCAN(network: str, last_time_mod: str | None = None):
    """
    Retrieves newer data for SCAN / SNOTEL since the date of the last data retrieval.

    Parameters
    ----------
    network : str
        name of network to return, SCAN or SNOTEL
    last_time_mod : str, optional
        time of last modification

    Returns
    -------
    None
    """

    today = datetime.now(timezone.utc).date()

    if last_time_mod is None:
        last_time_mod = get_last_date(f"1_raw_wx/{network}/")

    if last_time_mod < download_date:
        print(f"Downloading {network} data from {last_time_mod} to {download_date}.")

        # handling for SNOTEL naming
        if network == "SCAN":
            get_scan_station_data(
                WECC_TERR,
                WECC_MAR,
                start_date=str(last_time_mod),
                end_date=str(download_date),
                networks=[network],
                fileext=str(today),
            )
        elif network == "SNOTEL":
            usda_network = "SNTL"
            get_scan_station_data(
                WECC_TERR,
                WECC_MAR,
                start_date=str(last_time_mod),
                end_date=str(download_date),
                networks=[usda_network],
                fileext=str(today),
            )
    else:
        print(f"{network} station files up to date.")

    return None


if __name__ == "__main__":
    # option to custom select start date of download (format: datetime(year, month, day).date())
    last_time_mod = None

    update_asosawos(last_time_mod)
    # update_cimis(last_time_mod)
    # update_cw3e(last_time_mod)
    # update_hads(last_time_mod)
    # update_madis(network = 'CAHYDRO', last_time_mod) # modify network name
    # update_maritime(network = "MARITIME", last_time_mod) # modify network name: MARITIME or NDBC
    # update_otherisd(last_time_mod)
    # update_SCAN(network = "SCAN", last_time_mod) # modify network name: SCAN or SNOTEL
