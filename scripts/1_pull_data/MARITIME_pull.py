"""
MARITIME_pull.py

This script downloads MARITIME & NDBC data from NDBC using http.
Approach:
(1) Download data using station list.

Functions
---------
- get_maritime_station_ids: Returns list of station ids as CSV
- get_maritime: Read in buoys / C-MAN stations via HTTP access.
- get_maritime_update: Read in updated buoys using HTTP access.

Intended Use
------------
Retrieves raw data for an individual network, all variables, all times. Organized by station, with 1 file per year.

Notes
-----
1. This function assumes users have configured the AWS CLI such that their access key / secret key pair are stored in ~/.aws/credentials.
See https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html for guidance.
"""

import requests
from datetime import datetime, date
import pandas as pd
import boto3
from io import StringIO
from shapely.geometry import Point
import geopandas as gp
from geopandas.tools import sjoin
import numpy as np

from calc_pull import get_wecc_poly

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes
BUCKET_NAME = "wecc-historical-wx"
DIRECTORY_MAR = "1_raw_wx/MARITIME/"
DIRECTORY_NDBC = "1_raw_wx/NDBC/"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


def get_maritime_station_ids(
    terrpath: str, marpath: str, directory_mar: str, directory_ndbc: str
) -> pd.DataFrame:
    """Returns list of station ids as CSV

    Parameters
    ----------
    terrpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    marpath : str
        shapefiles for maritime and terrestrial WECC boundaries
    directory_mar : str
        directory path for MARITIME stations
    directory_ndbc : str
        directory path for NDBC stations

    Returns
    -------
    weccstations : pd.DataFrame
        buoy / CMAN stations within WECC

    Notes
    -----
    Gliders (non-stationary) location reported here as 30 N -90 W
    """

    # Identify NDBC buoys or C-MAN stations
    url = "https://www.ndbc.noaa.gov/data/stations/station_table.txt"
    r = requests.get(url)
    lines = r.content.split(b"\n")
    df = []

    for line in lines:
        row = line.split(b"|")
        row = [x.decode("utf-8") for x in row]
        row = list(map(str.strip, row))
        df.append(row)

    stations = pd.DataFrame(
        columns=[
            "STATION_ID",
            "OWNER",
            "TTYPE",
            "HULL",
            "NAME",
            "PAYLOAD",
            "LOCATION",
            "TIMEZONE",
            "FORECAST",
            "NOTE",
        ],
        data=df[1:],
    )
    stations = stations.dropna()
    stations = stations.drop(
        columns=["TTYPE", "HULL", "PAYLOAD", "FORECAST", "TIMEZONE"]
    )

    # Use spatial geometry to only keep points within WECC marine/terrestrial region (should only be marine)
    # Locations listed as one value: "30.000 N 90.000 W" as an example, split on spaces
    mar_lon = []
    mar_lat = []
    location = stations["LOCATION"]
    for line in location:
        line = line.split(" ")
        if "S" in line:
            mar_lat.append(-1 * float(line[0]))
        else:
            mar_lat.append(float(line[0]))

        if "E" in line:
            mar_lon.append(float(line[2]))
        else:
            mar_lon.append(-1 * float(line[2]))

    stations["LATITUDE"] = mar_lat
    stations["LONGITUDE"] = mar_lon

    # Zip lat lon coords and convert to geodataframe
    geometry = [Point(xy) for xy in zip(stations["LONGITUDE"], stations["LATITUDE"])]
    weccgeo = gp.GeoDataFrame(stations, crs="EPSG:4326", geometry=geometry)

    # get bbox of WECC to use to filter stations against
    t, m, bbox = get_wecc_poly(terrpath, marpath)

    # Get terrestrial stations
    weccgeo = weccgeo.to_crs(t.crs)
    # Only keep stations in terrestrial WECC region
    terwecc = sjoin(weccgeo.dropna(), t, how="left")
    terwecc = terwecc.dropna()

    # Get marine stations, keeping stations in maritime region
    marwecc = sjoin(weccgeo.dropna(), m, how="left")
    marwecc = marwecc.dropna()

    weccstations = pd.concat(
        [terwecc, marwecc], ignore_index=True, sort=False
    ).drop_duplicates(["STATION_ID"], keep="first")

    # drop columns and rename in_wecc to in_terr
    weccstations.drop(
        [
            "OBJECTID_1",
            "OBJECTID",
            "Shape_Leng",
            "geometry",
            "FID_WECC_B",
            "BUFF_DIST",
            "index_right",
        ],
        axis=1,
        inplace=True,
    )
    weccstations.rename(
        columns={"in_WECC": "in_terr_wecc", "in_marine": "in_mar_wecc"}, inplace=True
    )

    # Identify which buoys are moored buoys (46xxx) and which are C-MAN/water level obs network/other
    network = []
    for item in weccstations["STATION_ID"]:
        if item[:2] == "46":
            network.append("NDBC")
        else:
            network.append("MARITIME")
    weccstations["NETWORK"] = network

    # Moves Note column to last column because of weird lat-lon column overwriting -- metadata issue for cleaning
    weccstations = weccstations.reindex(
        columns=[col for col in weccstations.columns if col != "NOTE"] + ["NOTE"]
    )

    # Splits dataframe into respective networks - this is somewhat unnecessary
    maritime_network = weccstations[weccstations["NETWORK"] == "MARITIME"]
    ndbc_network = weccstations[weccstations["NETWORK"] == "NDBC"]

    # Write stations to respective AWS bucket - requires separate buffers
    mar_buffer = StringIO()
    maritime_network.to_csv(mar_buffer)
    content = mar_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME, Body=content, Key=directory_mar + "stationlist_MARITIME.csv"
    )

    ndbc_buffer = StringIO()
    ndbc_network.to_csv(ndbc_buffer)
    content = ndbc_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME, Body=content, Key=directory_ndbc + "stationlist_NDBC.csv"
    )

    # Purposely returning the full weccstations, instead of two separte dfs for get_maritime function
    return weccstations


def get_maritime(stations: pd.DataFrame, network: str, years: list[str] | None = None):
    """Read in buoys / C-MAN stations via HTTP access.

    Parameters
    ----------
    stations : pd.DataFrame
        stations to retrieve
    network : str
        name of network to return
    years : list(int)
        years to subset for

    Returns
    -------
    None
    """

    # Set up error handling df
    errors = {"Station ID": [], "Time": [], "Error": []}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    # HTTP access instead of FTP: https://www.ndbc.noaa.gov/data/historical/stdmet/
    directory = "1_raw_wx/" + network + "/"
    dir_stations = stations.loc[stations["NETWORK"] == network]

    # Default for years
    if years is None:
        # Get list of years from 1980 to current year
        years = list(map(str, range(1980, datetime.now().year + 1)))
    else:
        years = years

    # Identifies which stations are owned by Canadian Dept of Environment and Climate Change (not qaqc'd by NDBC)
    canadian_owners = list(stations.loc[stations["OWNER"] == "C"]["STATION_ID"])
    # Targeted pull fix for "CM" owner, if needed in future
    # cm_owner = list(stations.loc[stations['OWNER'] == 'CM']["STATION_ID"])

    for year in years:
        if year < str(datetime.now().year):
            for filename in dir_stations["STATION_ID"]:
                try:
                    if filename in canadian_owners:
                        # Environment and Climate Change Canadian Moored buoy
                        url = f"https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/csvData/c{filename}_csv.zip"
                        # all years in one file, including current year
                        s3_obj = s3.Object(
                            BUCKET_NAME, directory + f"{filename}_csv.zip"
                        )

                    else:
                        # NOAA NDBC Buoy archive
                        url = f"https://www.ndbc.noaa.gov/data/historical/stdmet/{filename}h{year}.txt.gz"
                        s3_obj = s3.Object(
                            BUCKET_NAME,
                            directory + f"{filename}h{year}.txt.gz",
                        )

                    with requests.get(url, stream=True) as r:
                        if r.status_code == 404:
                            # Catches any stations that don't have specific years, could be cleaner potentially
                            continue
                        elif r.status_code == 200:
                            s3_obj.put(Body=r.content)
                            # Note: The Canadian buoy all-years-file still gets downloaded/overwritten for len(years) times, where it could just be downloaded once
                            print(f"Saving data for station {filename} for {year}")
                        else:
                            errors["Station ID"].append(filename)
                            errors["Time"].append(end_api)
                            errors["Error"].append(r.status_code)
                            print(f"Error: {r.status_code}")

                except Exception as e:
                    print(f"Error: {e}")

        else:
            # only applicable for the NDBC stored data (non-Canadian source)
            for filename in dir_stations["STATION_ID"]:
                try:
                    year = int(year)
                    start_month = date(year, 1, 1)  # Jan
                    current_month = date.today()
                    i = current_month.month

                    while i >= start_month.month:
                        current_date = date(year, i, 1)
                        x = current_date.strftime("%b")
                        # NDBC folder names is 3-letter month code
                        url = f"https://www.ndbc.noaa.gov/data/stdmet/{x}/{filename}{i}{year}.txt.gz"
                        s3_obj = s3.Object(
                            BUCKET_NAME,
                            directory + f"{filename}{i}{year}.txt.gz",
                        )

                        with requests.get(url, stream=True) as r:
                            if r.status_code == 404:
                                pass
                            elif r.status_code == 200:
                                s3_obj.put(Body=r.content)
                                print(
                                    f"Saving data for station {filename} for {x} {year}"
                                )
                            else:
                                errors["Station ID"].append(filename)
                                errors["Time"].append(end_api)
                                errors["Error"].append(r.status_code)
                                print(f"Error: {r.status_code}")

                        i -= 1

                except Exception as e:
                    print(f"Error: {e}")

    # Write errors to csv with respective networks
    error_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(error_buffer)
    content = error_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"errors_{network}_{end_api}.csv",
    )

    return None


def get_maritime_update(
    stations: pd.DataFrame,
    network: str,
    start_date: str | None = None,
    end_date: str | None = None,
):
    """Read in updated buoys using HTTP access.

    Can only filter download by start and end month, so we download all
    data for month even if only a few days are requested.

    Parameters
    ----------
    station : pd.DataFrame
        stations to retrieve
    network : str
        name of network to retrieve
    start_date : str, optional
        date to start filtering by
    end_date : str, optional
        date to end filtering by

    Returns
    -------
    None
    """

    # Set up error handling df
    errors = {"Station ID": [], "Time": [], "Error": []}

    # Set end time to be current time at beginning of download
    end_api = datetime.now().strftime("%Y%m%d%H%M")

    # HTTP access instead of FTP: https://www.ndbc.noaa.gov/data/historical/stdmet/
    directory = "1_raw_wx/" + network + "/"
    dir_stations = stations.loc[stations["NETWORK"] == network]

    # Get range of years to download
    if start_date is None:
        start_month = date(1980, 1, 1).month  # Jan
        if end_date is None:
            years = list(map(str, range(1980, int(datetime.now().year) + 1)))
        else:
            years = list(map(str, range(1980, int(end_date[0:4]) + 1)))
            end_month = datetime.strptime(end_date, "%Y-%m-%d").month
    else:
        start_year = int(start_date[0:4])
        start_month = datetime.strptime(start_date, "%Y-%m-%d").month
        if end_date is None:
            years = list(map(str, range(start_year, int(datetime.now().year) + 1)))
        else:
            years = list(map(str, range(start_year, int(end_date[0:4]) + 1)))
            end_month = datetime.strptime(end_date, "%Y-%m-%d").month

    # Set up cross year flag
    if int(datetime.now().strftime("%j")) <= 45:
        # if download occurring in first 45 days of year, and includes previous year's data
        if (str(datetime.now().year - 1)) in years:
            # get list of months from start to end
            date_list = pd.period_range(start=start_date, end=end_date, freq="M")
            month_list = [month.strftime("%b") for month in date_list]
            # Set flag to be true
            cross_year = 1

    # Identifies which stations are owned by Canadian Dept of Environment and Climate Change (not qaqc'd by NDBC)
    canadian_owners = list(stations.loc[stations["OWNER"] == "C"]["STATION_ID"])

    for year in years:
        if year < str(datetime.now().year):
            # If not in current year
            for filename in dir_stations["STATION_ID"]:
                try:
                    # Try to get station txt.gz.
                    if filename in canadian_owners:
                        # Environment and Climate Change Canadian Moored buoy
                        url = f"https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/FBYEARS/{filename}/{filename}_{year}.zip"
                        # all years file format
                        s3_obj = s3.Object(
                            BUCKET_NAME, directory + f"{filename}_{year}_csv.zip"
                        )

                    else:
                        # NOAA NDBC Buoy archive
                        url = f"https://www.ndbc.noaa.gov/data/historical/stdmet/{filename}h{year}.txt.gz"
                        s3_obj = s3.Object(
                            BUCKET_NAME, directory + f"{filename}h{year}.txt.gz"
                        )

                    with requests.get(url, stream=True) as r:
                        if r.status_code == 404:
                            # Catches any stations that don't have specific years, could be cleaner potentially
                            continue
                        elif r.status_code == 200:
                            s3_obj.put(Body=r.content)
                            # Note: The Canadian buoy all-years-file still gets downloaded/overwritten for len(years) times, where it could just be downloaded once
                            print("Saving data for station {filename} for {year}")
                        else:
                            errors["Station ID"].append(filename)
                            errors["Time"].append(end_api)
                            errors["Error"].append(r.status_code)
                            print(f"Error: {r.status_code}")

                except Exception as e:
                    print(f"Error: {e}")

        else:
            # if year in present year, download by month. only applicable for the NDBC stored data (non-Canadian source)
            for filename in dir_stations["STATION_ID"]:
                try:
                    # Environment and Climate Change Canadian Moored buoy
                    if filename in canadian_owners:
                        # Data for current year
                        url = f"https://www.meds-sdmm.dfo-mpo.gc.ca/alphapro/wave/waveshare/FBYEAR2DATE/{filename}_Y2D.zip"
                        # Overwrite current year data file
                        s3_obj = s3.Object(
                            BUCKET_NAME, directory + f"{filename}_{year}_csv.zip"
                        )
                    else:
                        year = int(year)
                        start_month = date(year, 1, 1)  # Jan
                        current_month = date.today()
                        i = current_month.month

                        while i >= start_month.month:
                            current_date = date(year, i, 1)
                            # NDBC folder names is 3-letter month code
                            x = current_date.strftime("%b")
                            url = f"https://www.ndbc.noaa.gov/data/stdmet/{x}/{filename}{i}{year}.txt.gz"
                            s3_obj = s3.Object(
                                BUCKET_NAME, directory + f"{filename}{i}{year}.txt.gz"
                            )

                            with requests.get(url, stream=True) as r:
                                if r.status_code == 404:
                                    pass
                                elif r.status_code == 200:
                                    s3_obj.put(Body=r.content)
                                    print(
                                        f"Saving data for station {filename} for {x} {year}"
                                    )
                                else:
                                    errors["Station ID"].append(filename)
                                    errors["Time"].append(end_api)
                                    errors["Error"].append(r.status_code)
                                    print(f"Error: {r.status_code}")

                        i -= 1

                except Exception as e:
                    print(f"Error: {e}")

    # Additionally, if data includes cross-year data, manually download dates from previous year from monthly data
    if cross_year == 1:
        for filename in dir_stations["STATION_ID"]:
            if filename in canadian_owners:
                continue
            else:
                try:
                    # If list spans two years, remove new year months.
                    if "January" in month_list:
                        # Get index of January
                        prior_year_ind = month_list.index("January")
                        prior_year_months = month_list[:prior_year_ind]

                    else:
                        prior_year_months = month_list

                    # Then, download all months in month list
                    for i in prior_year_months:
                        # Get month in number form
                        month_nom = datetime.strptime(i, "%b").month

                        # Last three months have different file format. Reformat to match
                        # If files in folder have just station id and no time in their filename, they are not final data and the script will skip them automatically (line 371)
                        if month_nom > 9:
                            if month_nom == 10:
                                month_nom = "a"
                            if month_nom == 11:
                                month_nom = "b"
                            if month_nom == 12:
                                month_nom = "c"

                        url = f"https://www.ndbc.noaa.gov/data/stdmet/{i}/{filename}{month_nom}{year}.txt.gz"

                        s3_obj = s3.Object(
                            BUCKET_NAME,
                            directory + f"{filename}{month_nom}{year}.txt.gz",
                        )

                        with requests.get(url, stream=True) as r:
                            if r.status_code == 404:
                                pass
                            elif r.status_code == 200:
                                s3_obj.put(Body=r.content)
                                print(
                                    f"Saving data for station {filename} for {i} {year}"
                                )
                            else:
                                errors["Station ID"].append(filename)
                                errors["Time"].append(end_api)
                                errors["Error"].append(r.status_code)
                                print(f"Error: {r.status_code}")

                except Exception as e:
                    print(f"Error: {e}")

    # Write errors to csv with respective networks
    error_buffer = StringIO()
    errors = pd.DataFrame(errors)
    errors.to_csv(error_buffer)
    content = error_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME,
        Body=content,
        Key=directory + f"errors_{network}_{end_api}.csv",
    )

    return None


if __name__ == "__main__":
    # Select either "MARITIME" or "NDBC" as network of choice to download for "network_to_run"
    network_to_run = "MARITIME"
    stations = get_maritime_station_ids(
        WECC_TERR, WECC_MAR, DIRECTORY_MAR, DIRECTORY_NDBC
    )
    get_maritime(stations, network_to_run, years=None)
