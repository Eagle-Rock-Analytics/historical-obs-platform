"""
MADIS_pull.py

This script scrapes MADIS network data for ingestion into the Historical Observations Platform via API.

Functions
---------
- get_network_metadata: Calls Synoptic API to get network metadata and save to AWS
- get_madis_metadata: produces a list of station IDs filtered by bounding box and network
- get_madis_station_csv: Download network data from Synoptic API. 
- get_madis_station_csv_update: Download updated network data from Synoptic API. 
- madis_pull: Pulls the raw data from a MADIS network.
- madis_update: Pulls the updated raw data from a MADIS network.

Intended Use
------------
Retrieves raw data for an individual network, all variables, all times. Organized by station, with 1 file per year.

Notes
-----
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

import requests
import pandas as pd
from datetime import datetime
import re
import boto3
from io import StringIO
import config  # Import API keys.

from calc_pull import get_wecc_poly

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes
BUCKET_NAME = "wecc-historical-wx"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"
RAW_PATH = "1_raw_wx/"


def get_network_metadata(token: str) -> pd.DataFrame:
    """Calls Synoptic API to get network metadata and save to AWS

    Parameters
    ----------
    token : str
        Synoptic API token (private)

    Returns
    -------
    networks : pd.DataFrame
        table of station metadata
    """

    # Produce table of station metadata from Synoptic API
    r = requests.get(f"https://api.synopticdata.com/v2/networks?token={token}").json()
    networks = pd.DataFrame(r["MNET"])

    # Sort by number of stations
    networks = networks.sort_values("REPORTING_STATIONS", ascending=False)

    # Save to AWS bucket
    csv_buffer_err = StringIO()
    networks.to_csv(csv_buffer_err)
    content = csv_buffer_err.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME, Body=content, Key=RAW_PATH + "madis_network_metadata.csv"
    )

    return networks


def get_madis_metadata(
    token: str, terrpath: str, marpath: str, networkid: str, directory
) -> pd.DataFrame | None:
    """Returns metadata for stations from Synoptic API.

    Parameters
    ----------
    token : str
        Synoptic API token (private)
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    networkid : str
        network name identifier
    directory : str
        AWS directory folder to save

    Returns
    -------
    ids : pd.DataFrame
        stations returned per network
    """

    try:
        t, m, bbox = get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0, :].tolist()  # [lonmin,latmin,lonmax,latmax]
        bbox_api = ",".join([str(elem) for elem in bbox_api])
        # Access station metadata to get list of IDs in bbox and network
        # Using: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
        url = "https://api.synopticdata.com/v2/stations/metadata?token={token}&network={networkid}&bbox={bbox_api}&recent=20&output=json"

        request = requests.get(url).json()
        ids = []
        station_list = pd.DataFrame(request["STATION"])
        # Split Period of Record column
        station_list = pd.concat(
            [station_list, station_list["PERIOD_OF_RECORD"].apply(pd.Series)], axis=1
        )
        station_list = station_list.drop("PERIOD_OF_RECORD", axis=1)

        # Sort by start date (note some stations return 'None' here)
        ids = station_list[["STID", "start"]].sort_values("start")
        # Reformat date to match API format
        ids["start"] = pd.to_datetime(ids["start"], format="%Y-%m-%dT%H:%M:%SZ")
        ids["start"] = ids["start"].dt.strftime("%Y%m%d%H%M")

        # Save station list to AWS
        csv_buffer_err = StringIO()
        station_list.to_csv(csv_buffer_err)
        content = csv_buffer_err.getvalue()
        # Get network name from directory name
        networkname = directory.replace("1_raw_wx/", "")
        networkname = networkname.replace("/", "")
        networkname = networkname.replace(" ", "")

        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=directory + f"stationlist_{networkname}.csv",
        )

        return ids

    except Exception as e:
        print(f"Error: {e}")
        return None


def get_madis_station_csv(
    token: str,
    ids: pd.DataFrame,
    directory: str,
    start_date: str | None = None,
    **options: dict[str, bool],
):
    """Download network data from Synoptic API.

    Parameters
    ----------
    token : str
        Synoptic API token (private)
    ids : pd.DataFrame
        network ids
    directory : str
        AWS directory to save
    start_date : str, optional
        specific start date to subset for; format "YYYYMMDDHHMM"
    options : dict[bool, str]
        timeout : bool
            will identify and download any station data that timed out the API request
        extension : str
            looking for csv files
        round : str
            prefix (directory) name within API request

    Returns
    -------
    None
    """

    # Set up error handling df.
    errors = {"Station ID": [], "Time": [], "Error": []}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    for index, id in ids.iterrows():
        # Set start date
        if start_date is None:
            # Adding error handling for NaN dates
            if id["start"] == "NaN" or pd.isna(id["start"]):
                # If any of the stations in the chunk have a start time of NA, pull data from 01-01-1980 OH:OM and on
                start_api = "198001010000"
                # To check: all stations seen so far with NaN return empty data. If true universally, can skip these instead of saving.

            else:
                start_api = id["start"]
                # If start_api starts prior to 1980, set manually to be 1980-01-01.
                if start_api[0:4] < "1980":
                    start_api = "198001010000"

        else:
            start_api = start_date

        # If station file is 2_/3_ etc., get station name
        if options.get("timeout") == True:
            id["STID"] = id["STID"].split("_")[-1]

        # Generate URL
        # Note: decision here to use full flag suite of MesoWest and Synoptic data.
        # See Data Checks section here for more information: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
        id_stn = id["STID"]
        url = f"https://api.synopticdata.com/v2/stations/timeseries?token={token}&stid={id_stn}&start={start_api}&end={end_api}&output=csv&qc=on&qc_remove_data=off&qc_flags=on&qc_checks=synopticlabs,mesowest"

        try:
            # request = requests.get(url)
            s3_obj = s3.Object(BUCKET_NAME, directory + f"{id_stn}.csv")

            # If **options timeout = True, save file as 2_STID.csv
            if options.get("timeout") == True:
                prefix = options.get("round")
                s3_obj = s3.Object(BUCKET_NAME, directory + f"{prefix}_{id_stn}.csv")

            with requests.get(url, stream=True) as r:
                # If API call returns a response
                if r.status_code == 200:
                    # If error response returned. Note that this is formatted differently depending on error type
                    if "RESPONSE_MESSAGE" in r.text:
                        # Get error message and clean
                        error_text = str(
                            re.search("(RESPONSE_MESSAGE.*)", r.text).group(0)
                        )
                        error_text = re.sub("RESPONSE_MESSAGE.: ", "", error_text)
                        error_text = re.sub(",.*", "", error_text)
                        error_text = re.sub('"', "", error_text)

                        # Append rows to dictionary
                        errors["Station ID"].append(id["STID"])
                        errors["Time"].append(end_api)
                        errors["Error"].append(error_text)
                        continue

                    else:
                        s3_obj.put(Body=r.content)
                        print(f"Saving data for station {id_stn}")

                else:
                    errors["Station ID"].append(id["STID"])
                    errors["Time"].append(end_api)
                    errors["Error"].append(r.status_code)
                    print(f"Error: {r.status_code}")

        except Exception as e:
            print(f"Error: {e}")

    # Write errors to csv for AWS
    csv_buffer_err = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer_err)
    content = csv_buffer_err.getvalue()
    # Get network name from directory name
    networkname = directory.replace("1_raw_wx/", "")
    networkname = networkname.replace("/", "")

    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"errors_{networkname}_{end_api}.csv",
    )

    return None


def get_madis_station_csv_update(
    token: str,
    ids: pd.DataFrame,
    directory: str,
    start_date: str | None = None,
    end_date: str | None = None,
    **options: dict[bool, str],
):
    """Download updated network data from Synoptic API.

    Parameters
    ----------
    token : str
        Synoptic API token (private)
    ids : pd.DataFrame
        station ids
    directory : str
        AWS directory to save
    start_date : str, optional
        date to subset start
    end_date : str, optional
        date to subset end
    options : dict[bool, str]
        timeout : bool
            will identify and download any station data that timed out the API request
        extension : str
            looking for csv files
        round : str
            prefix (directory) name within API request

    Returns
    -------
    None
    """
    # Set up error handling df
    errors = {"Station ID": [], "Time": [], "Error": []}

    # Set end time to be current time at beginning of download
    if end_date is None:
        end_api = datetime.now().strftime("%Y%m%d%H%M")
    else:
        end_api = datetime.strptime(end_date, "%Y-%m-%d").strftime("%Y%m%d%H%M")

    for index, id in ids.iterrows():  # For each station
        # Set start date
        if start_date is None:
            # Adding error handling for NaN dates
            if id["start"] == "NaN" or pd.isna(id["start"]):
                # If any of the stations in the chunk have a start time of NA, pull data from 01-01-1980 OH:OM and on.
                start_api = "198001010000"
                # To check: all stations seen so far with NaN return empty data. If true universally, can skip these instead of saving.

            else:
                start_api = id["start"]
                # If start_api starts prior to 1980, set manually to be 1980-01-01
                if start_api[0:4] < "1980":
                    start_api = "198001010000"
        else:
            start_api = datetime.strptime(start_date, "%Y-%m-%d").strftime("%Y%m%d%H%M")

        # If station file is 2_/3_ etc., get station name
        if options.get("timeout") == True:
            id["STID"] = id["STID"].split("_")[-1]

        # Generate URL
        # Note: decision here to use full flag suite of MesoWest and Synoptic data.
        # See Data Checks section here for more information: https://developers.synopticdata.com/mesonet/v2/stations/timeseries/
        id_stn = id["STID"]
        url = f"https://api.synopticdata.com/v2/stations/timeseries?token={token}&stid={id_stn}&start={start_api}&end={end_api}&output=csv&qc=on&qc_remove_data=off&qc_flags=on&qc_checks=synopticlabs,mesowest"

        try:
            # Set up filename, given parameters
            filename = directory + f"{id_stn}.csv"
            if options.get("extension") is not None:
                if options.get("timeout") == True:
                    prefix = options.get("round")
                    ext = options.get("extension")
                    filename = directory + f"{prefix}_{id_stn}_{ext}.csv"
                else:
                    filename = directory + f"{id_stn}_{ext}.csv"

            elif options.get("timeout") == True:
                prefix = options.get("round")
                filename = directory + f"{prefix}_{id_stn}.csv"

            s3_obj = s3.Object(BUCKET_NAME, filename)

            with requests.get(url, stream=True) as r:
                # If API call returns a response
                if r.status_code == 200:
                    # If error response returned. Note that this is formatted differently depending on error type
                    if "RESPONSE_MESSAGE" in r.text:
                        # Get error message and clean
                        error_text = str(
                            re.search("(RESPONSE_MESSAGE.*)", r.text).group(0)
                        )
                        error_text = re.sub("RESPONSE_MESSAGE.: ", "", error_text)
                        error_text = re.sub(",.*", "", error_text)
                        error_text = re.sub('"', "", error_text)

                        # Append rows to dictionary
                        errors["Station ID"].append(id["STID"])
                        errors["Time"].append(end_api)
                        errors["Error"].append(error_text)
                        continue
                    else:
                        s3_obj.put(Body=r.content)
                        print(f"Saving data for station {id_stn}")

                else:
                    errors["Station ID"].append(id["STID"])
                    errors["Time"].append(end_api)
                    errors["Error"].append(r.status_code)
                    print(f"Error: {r.status_code}")

        except Exception as e:
            print(f"Error: {e}")

    # Write errors to csv for AWS
    csv_buffer_err = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(csv_buffer_err)
    content = csv_buffer_err.getvalue()
    # Get network name from directory name
    networkname = directory.replace("1_raw_wx/", "")
    networkname = networkname.replace("/", "")
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"errors_{networkname}_{end_api}.csv",
    )

    return None


def madis_pull(
    token: str, networks: list[str] | None = None, pause: bool | None = None
):
    """Pulls the raw data from a MADIS network.

    Parameters
    ----------
    token : str
        Synoptic API token (private)
    networks : list[str]
        list of network names
    pause : bool, optional
        If true, prompts user to indicate "yes" to continue before downloading large networks.
        Automatically set to yes when networks is not specified

    Returns
    -------
    None
    """

    # If no networks provided, pull full list
    if networks is None:
        networkdf = get_network_metadata(token)
        # Remove restricted networks
        networkdf = networkdf[networkdf["TOTAL_RESTRICTED"] == 0]
        # Only keep ID, shortname, expected station #s
        networkdf = networkdf[["ID", "SHORTNAME", "REPORTING_STATIONS"]]
        # Set pause to be true (this will be a large runtime.)
        pause = True

    else:
        networkdf = get_network_metadata(token)
        mask = [
            any([re.search(r"\b" + kw + r"\b", r) for kw in networks])
            for r in networkdf["SHORTNAME"]
        ]
        networkdf = networkdf[mask]
        shortname = list(networkdf["SHORTNAME"])
        resp = input(
            f"Networks to be downloaded:{shortname}. Would you like to continue? (Y/N)"
        )
        if resp == "N":
            return
        if resp == "Y":
            pass
        else:
            print("Invalid response, exiting function.")
            return

    # By network, download data
    for index, row in networkdf.iterrows():
        dirname = row["SHORTNAME"]
        print(f"Downloading data for {dirname} network")
        if pause:
            # Set up pause function for large networks
            if row["REPORTING_STATIONS"] >= 1000:
                resp = input(
                    f"Warning: This network contains {row["REPORTING_STATIONS"]} stations. Are you ready to download it? (Y/N)"
                )
                if resp == "N":
                    print("Skipping to next network.")
                    continue
                if resp == "Y":
                    pass
                else:
                    print("Invalid response. Skipping network.")

        # clean up directory name
        dirname = dirname.replace(" ", "")
        dirname = dirname.replace("/", "-")
        directory = RAW_PATH + dirname + "/"

        # Except if the network is CWOP, then manually set to be CWOP
        if row["SHORTNAME"] == "APRSWXNET/CWOP":
            dirname = "CWOP"
            directory = RAW_PATH + "CWOP/"

        # Get list of station IDs and start date
        ids = get_madis_metadata(
            token=config.token,
            terrpath=WECC_TERR,
            marpath=WECC_MAR,
            networkid=row["ID"],
            directory=directory,
        )

        # Get station CSVs
        get_madis_station_csv(token=config.token, directory=directory, ids=ids)

    return None


def madis_update(
    token: str,
    networks: list[str] | None = None,
    pause: bool | None = None,
    start_date=None,
    end_date=None,
):
    """Pulls the updated raw data from a MADIS network.

    Parameters
    ----------
    token : str
        Synoptic API token (private)
    networks : list[str]
        list of network names
    pause : bool, optional
        If true, prompts user to indicate "yes" to continue before downloading large networks.
        Automatically set to yes when networks is not specified
    start_date : str, optional
        specific start date to subset for; format "YYYYMMDDHHMM"
    end_date : str, optional
        specific end date to subset for; format "YYYYMMDDHHMM"

    Returns
    -------
    None
    """

    # If no networks provided, pull full list
    if networks is None:
        networkdf = get_network_metadata(token)
        # Remove restricted networks
        networkdf = networkdf[networkdf["TOTAL_RESTRICTED"] == 0]
        # Only keep ID, shortname, expected station #s
        networkdf = networkdf[["ID", "SHORTNAME", "REPORTING_STATIONS"]]
        # Set pause to be true (this will be a large runtime)
        pause = True

    else:
        networkdf = get_network_metadata(token)
        mask = [
            any([re.search(r"\b" + kw + r"\b", r) for kw in networks])
            for r in networkdf["SHORTNAME"]
        ]
        networkdf = networkdf[mask]
        shortname_list = list(networkdf["SHORTNAME"])
        resp = input(
            f"Networks to be downloaded:{shortname_list}. Would you like to continue? (Y/N)"
        )
        if resp == "N":
            return
        if resp == "Y":
            pass
        else:
            print("Invalid response, exiting function.")
            return

    # By network, download data
    for index, row in networkdf.iterrows():
        dirname = row["SHORTNAME"]
        print(f"Downloading data for network: {dirname}")
        # set up pause for large networks
        if pause:
            if row["REPORTING_STATIONS"] >= 1000:
                resp = input(
                    f"Warning: This network contains {row["REPORTING_STATIONS"]} stations. Are you ready to download it? (Y/N)"
                )
                if resp == "N":
                    print("Skipping to next network.")
                    continue
                if resp == "Y":
                    pass
                else:
                    print("Invalid response. Skipping network.")

        # clean up directory name
        dirname = dirname.replace(" ", "")
        dirname = dirname.replace("/", "-")
        directory = RAW_PATH + dirname + "/"

        # Except if the network is CWOP, then manually set to be CWOP
        if dirname == "APRSWXNET/CWOP":
            dirname = "CWOP"
            directory = RAW_PATH + "CWOP/"

        # Get list of station IDs and start date
        ids = get_madis_metadata(
            token=config.token,
            terrpath=WECC_TERR,
            marpath=WECC_MAR,
            networkid=row["ID"],
            directory=directory,
        )

        # Get station CSVs
        extension = datetime.now().strftime("%Y-%m-%d")
        get_madis_station_csv_update(
            token=config.token,
            directory=directory,
            ids=ids,
            start_date=start_date,
            end_date=end_date,
            extension=extension,
        )

    return None


if __name__ == "__main__":
    madis_pull(config.token)
    # Note: this will download all non-restricted networks in MADIS. Specify to subset; example: madis_pull(config.token, networks=["CWOP"])
    # Warning: running this without specification may overwrite metadata for networks that are in MADIS, but downloaded directly.
