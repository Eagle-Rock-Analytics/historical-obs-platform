"""
SCANSNOTEL_pull.py

This script scrapes SNOTEL and SCAN network data for ingestion into the Historical Observations Platform via SOAP API.
It may in the future be extended to pull other hourly data from networks available through the platform (e.g. BOR, USGS)
Approach:
(1) get_wecc_poly generates a bounding box from WECC shapefiles.
(3) get_SCAN_stations identifies station ids and provides metadata for SCAN stations in WECC.
(3) get SCAN_station_data saves a csv for each station ID, with the option to select a start date for downloads
    (defaults to station start date) and networks of interest, and primary or secondary sensors.

Functions
---------
- get_SCAN_stations: Get SCAN/SNOTEL station list in WECC region from SOAP API and save to AWS.
- get_scan_station_data: Download USDA station data using SOAP API

Notes
-----
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

import pandas as pd
from datetime import datetime
import numpy as np
import boto3
from io import BytesIO, StringIO
from zeep import Client  # For calling SOAP APIs
from zeep.helpers import serialize_object

# import requests # if necessary

try:
    from calc_pull import get_wecc_poly
except RuntimeError as e:
    print(f"Error importing calc_pull: {e}")

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes
BUCKET_NAME = "wecc-historical-wx"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Connect to SCAN API
url = "https://wcc.sc.egov.usda.gov/awdbWebService/services?WSDL"
client = Client(url)


def get_SCAN_stations(
    terrpath: str, marpath: str, networks: list[str] | None = None
) -> pd.DataFrame | None:
    """
    Get SCAN station list in WECC region from SOAP API and save to AWS.

    Parameters
    ----------
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    networks : list[str]
        name of network(s)

    Returns
    -------
    station_metadata : pd.DataFrame
        station metadata
    """
    try:
        t, m, bbox = get_wecc_poly(terrpath, marpath)
        bbox_api = bbox.loc[0, :].tolist()  # [lonmin,latmin,lonmax,latmax]

        # Format bounding box for SOAP request
        request_data = {
            "minLatitude": str(bbox_api[1]),
            "maxLatitude": str(bbox_api[3]),
            "minLongitude": str(bbox_api[0]),
            "maxLongitude": str(bbox_api[2]),
            "logicalAnd": "true",  # and add other required parameter
        }

        # get station ids from bounding box
        stations = client.service.getStations(**request_data)

        # use station ids to return metadata
        request_data = {"stationTriplets": stations}

        # Call SOAP API
        metadata = client.service.getStationMetadataMultiple(**request_data)
        metadata = serialize_object(metadata)  # Reformat object
        station_metadata = pd.DataFrame(metadata)  # Convert to pandas df

        # If networks is none, set default networks
        if networks is None:
            networks = ["SNTL", "SCAN"]

        for i in networks:
            subset = station_metadata[
                station_metadata["stationTriplet"].str.contains(i)
            ]
            # Fix name differences between SNOTEL/SNTL
            if i == "SNTL":
                i = "SNOTEL"
            subdir = "1_raw_wx/{}/".format(i)

            # Save to AWS
            csv_buffer = StringIO()
            subset.to_csv(csv_buffer, index=False)
            content = csv_buffer.getvalue()
            s3_cl.put_object(
                Bucket=BUCKET_NAME,
                Body=content,
                Key=subdir + "stationlist_{}.csv".format(i),
            )

        return station_metadata

    except Exception as e:
        print(f"Error: {e}")
        return None


def get_scan_station_data(
    terrpath: str,
    marpath: str,
    start_date: str | None = None,
    end_date: str | None = None,
    stations: pd.Series | None = None,
    primary: bool = False,
    networks: list[str] | None = ["SNTL", "SCAN"],
    fileext: str | None = None,
):
    """Download USDA station data using SOAP API. Data is organized by station, by sensor.

    Parameters
    ----------
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    start_ddate : str, optional
        If none, download data starting on 1980-01-01. Otherwise, download data after start date (format "YYYY-MM-DD")
    end_date : str, optional
        If none, download data through present. Otherwise, download data after start date (format "YYYY-MM-DD") through end date
    stations : pd.Series, optional
        station metadata to download, required in format with "stationTriplet" and "beginDate"
    primary : bool, optional
        download primary data only, if primary is False, download secondary sensor data.
    networks : list[str], optional
        specify which network to download. If None, get all networks: ["SNTL", "SCAN"]
    fileext :
        file extension

    Returns
    -------
    None
    """

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    # Specify default networks
    if networks is None:
        networks = ["SNTL", "SCAN"]

    # Get stations metadata if not provided
    if stations is None:
        stations = get_SCAN_stations(terrpath, marpath, networks)
    stationTriplets = list(stations["stationTriplet"])

    # Get start and end dates
    if start_date is None:
        start_date = "1980-01-01"

    if end_date is None:
        end_date = datetime.now().strftime("%Y-%m-%d")  # yyyy-mm-dd

    # Select primary or secondary sensors
    if primary:  # If primary True
        ordinal = "1"
    else:
        ordinal = "2"

    # Select codes for desired variables
    sensorcodes = [
        "TAVG",  # Air temp avg
        "TOBS",  # Air temp observed
        "PREC",  # Precipitation accumulation
        "PRCP",  # Precipitation increment
        "PRCPSA",  # Precipitation increment snow-adjusted
        "PRES",  # Barometric pressure
        "DPTP",  # Dew point temperature
        "RHUM",  # Relative humidity
        "RHUMV",  # Relative humidity (avg),
        "SRAD",  # Solar radiation
        "SRADV",  # Solar radiation (avg)
        "SRADT",  # Solar radiation (total)
        "PVPV",  # Partial vapor pressure
        "SVPV",  # Saturated vapor pressure
        "WDIRV",  # Wind direction (avg)
        "WDIR",  # Wind direction
        "WSPD",  # Wind speed
        "WSPDV",  # Wind speed (avg)
    ]

    for i in networks:
        # Set up error handling df
        errors = {"Station ID": [], "Time": [], "Error": []}

        networkTriplets = [k for k in stationTriplets if i in k]
        if i == "SNTL":
            # Maintain consistency across methods
            directory = "1_raw_wx/SNOTEL/"
        else:
            # Define directory path for AWS
            directory = "1_raw_wx/" + i + "/"

        for j in networkTriplets:
            try:
                # Set up empty pd dataframe
                df = pd.DataFrame()

                for sensor in sensorcodes:
                    # Instantaneous:
                    # Required parameters:
                    if i in ["SNTL", "SCAN"]:
                        request_data = {
                            "stationTriplets": j,  # Station ID(s) # Subset for testing.
                            "elementCd": sensor,  # Sensor ID
                            "ordinal": ordinal,  # Primary or secondary sensors
                            "beginDate": start_date,  # Set earliest date
                            "endDate": end_date,  # Set end date as today
                            "filter": "ALL",  # Select time of day to be all
                            "unitSystem": "ENGLISH",  # Convert to consistent unit system
                        }
                        inst_data = client.service.getInstantaneousData(**request_data)
                    else:
                        request_data = {
                            "stationTriplets": j,  # Station ID(s) # Subset for testing.
                            "elementCd": sensor,  # Sensor ID
                            "ordinal": ordinal,  # Primary or secondary sensors
                            "beginDate": start_date,  # Set earliest date
                            "endDate": end_date,  # Set end date as today
                        }
                        inst_data = client.service.getHourlyData(**request_data)
                        # Note: hourly data response not yet tested (only empty dfs returned),
                        # but should be used for non snotel/scan networks.

                    # Reformat object and convert to pandas df
                    data = serialize_object(inst_data[0])
                    station_data = pd.DataFrame(data["values"])

                    # If no data returned for sensor, skip to next sensor
                    if station_data.empty:
                        continue
                    else:
                        # Add sensor code to column titles
                        station_data = station_data.add_prefix(sensor + "_")

                        if df.empty:
                            df = pd.concat([df, station_data])
                            timecol = sensor + "_time"
                            df["time"] = df[timecol]
                        else:
                            df = df.merge(
                                station_data,
                                left_on="time",
                                right_on=sensor + "_time",
                                how="outer",
                            )
                            timecol = sensor + "_time"

                            # If any times in timecol are na, update with times from newest data
                            df["time"] = np.where(
                                df["time"].isnull(), df[timecol], df["time"]
                            )

                # If df empty, skip to next station
                if df.empty:
                    print(f"No data found for station {j}. Skipping to next station.")
                    continue

                # Make time column first in order
                time = df.pop("time")
                df.insert(0, "time", time)

                # Sort by time
                df = df.sort_values(by="time")

                # Write df to csv
                csv_buffer = StringIO()
                df.to_csv(csv_buffer, index=False)
                content = csv_buffer.getvalue()
                if fileext is None:
                    s3_cl.put_object(
                        Bucket=BUCKET_NAME,
                        Body=content,
                        Key=directory + "{}.csv".format(j),
                    )
                else:
                    s3_cl.put_object(
                        Bucket=BUCKET_NAME,
                        Body=content,
                        Key=directory + "{}_{}.csv".format(j, fileext),
                    )
                print(f"Saved data for station {j} in network {i}")
            except Exception as e:
                print(e)
                errors["Station ID"].append(j)
                errors["Time"].append(end_api)
                errors["Error"].append(e)

        # Save errors to network folder
        csv_buffer = StringIO()
        errors = pd.DataFrame(errors)
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        if i == "SNTL":
            i = "SNOTEL"  # Fix name for errors file
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=directory + "errors_{}_{}.csv".format(i, end_api),
        )

    return None


if __name__ == "__main__":
    get_scan_station_data(WECC_TERR, WECC_MAR, BUCKET_NAME, networks=["SNTL"])

# Notes
# 1. Neither BOR nor USGS have hourly data for any of our variables of interest.
# 2. Code below downloads additional metadata lists (with different station names) for SCAN/SNOTEL stations.
# Refer to if names aren't matching correctly down the line.

## SCAN
# url = 'https://wcc.sc.egov.usda.gov/nwcc/yearcount?network=scan&counttype=listwithdiscontinued&state='
# html = requests.get(url).content
# df_list = pd.read_html(html)
# scandf = df_list[-1]
# scandf['huccode'] = scandf['huc'].str.split("(").str[-1]
# scandf['huccode'] = scandf['huccode'].str.replace(")", "")
# # Remove all scan stations without huc code
# scandf['huccode'].replace('', np.nan, inplace=True)
# scandf.dropna(subset = ['huccode'], inplace = True) # Dropped one station.

## SNOTEL
# url = 'https://wcc.sc.egov.usda.gov/nwcc/yearcount?network=sntl&counttype=listwithdiscontinued&state='
# html = requests.get(url).content
# df_list = pd.read_html(html)
# sntldf = df_list[-1]
# sntldf['huccode'] = sntldf['huc'].str.split("(").str[-1]
# sntldf['huccode'] = sntldf['huccode'].str.replace(")", "")
# # Remove all scan stations without huc code
# sntldf['huccode'].replace('', np.nan, inplace=True)
# sntldf.dropna(subset = ['huccode'], inplace = True) # Dropped two stations.
