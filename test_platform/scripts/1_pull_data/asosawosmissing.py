"""
SCRIPT IS DEPRACATED AS OF NOV 14 2022. USE pull_qa.py INSTEAD.
"""
## Step 0: Environment set-up
# Import libraries
from ftplib import FTP
from datetime import datetime, timezone
import pandas as pd
from shapely.geometry import Point
import pandas as pd
import geopandas as gp
from geopandas.tools import sjoin
import boto3  # For AWS integration.
from io import BytesIO, StringIO
import calc_pull
import requests
import numpy as np
import config

# Set envr variables

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")

bucket_name = "wecc-historical-wx"
directory = "1_raw_wx/ASOSAWOS/"

# Set paths to WECC shapefiles in AWS bucket.
wecc_terr = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Set state shortcodes
states = [
    "AK",
    "AL",
    "AR",
    "AZ",
    "CA",
    "CO",
    "CT",
    "DC",
    "DE",
    "FL",
    "GA",
    "HI",
    "IA",
    "ID",
    "IL",
    "IN",
    "KS",
    "KY",
    "LA",
    "MA",
    "MD",
    "ME",
    "MI",
    "MN",
    "MO",
    "MS",
    "MT",
    "NC",
    "ND",
    "NE",
    "NH",
    "NJ",
    "NM",
    "NV",
    "NY",
    "OH",
    "OK",
    "OR",
    "PA",
    "RI",
    "SC",
    "SD",
    "TN",
    "TX",
    "UT",
    "VA",
    "VT",
    "WA",
    "WI",
    "WV",
    "WY",
]


# Function to write FTP data directly to AWS S3 folder.
# ftp here is the current ftp connection
# file is the filename
# directory is the desired path (set of folders) in AWS
def ftp_to_aws(ftp, file, directory):
    r = BytesIO()
    ftp.retrbinary("RETR " + file, r.write)
    r.seek(0)
    s3.upload_fileobj(r, bucket_name, directory + file)
    print("{} saved".format(file))  # Optional.
    r.close()  # Close file


# Function to download and parse ASOS and AWOS station lists (.txt).
# Source: https://www.ncei.noaa.gov/access/homr/reports/platforms
# Inputs: none.
# Outputs: one asosawos-stations.csv file and a stations object.
# This function is called internally in get_wecc_stations().
def get_asosawos_stations():
    # Get AWOS stations.
    awosurl = "https://www.ncei.noaa.gov/access/homr/file/awos-stations.txt"
    awosr = requests.get(awosurl)
    lines = awosr.content.split(b"\n")  # Split by line

    lines = lines[4:]  # Skip the first 4 lines
    df = []
    for line in lines:  # Parse data manually.
        ncdcid = line[0:8]
        wban = line[9:14]
        coopid = line[15:21]
        call = line[22:26]
        name = line[27:57]
        country = line[58:78]
        st = line[79:81]
        county = line[82:112]
        lat = line[113:122]
        lon = line[123:133]
        elev = line[134:140]
        utc = line[141:146]
        stntype = line[147:197]

        row = [
            ncdcid,
            wban,
            coopid,
            call,
            name,
            country,
            st,
            county,
            lat,
            lon,
            elev,
            utc,
            stntype,
        ]
        row = [x.decode("utf-8") for x in row]  # Convert to string
        row = [x.strip() for x in row]  # Strip whitespace
        df.append(row)

    # # Convert to pandas
    stations = pd.DataFrame(
        df,
        columns=[
            "NCDCID",
            "WBAN",
            "COOPID",
            "CALL",
            "NAME",
            "COUNTRY",
            "ST",
            "COUNTY",
            "LAT",
            "LON",
            "ELEV",
            "UTC",
            "STNTYPE",
        ],
    )

    # Get ASOS stations.
    asosurl = "https://www.ncei.noaa.gov/access/homr/file/asos-stations.txt"
    asosr = requests.get(asosurl)
    lines = asosr.content.split(b"\n")  # Split by line
    # Skip the first 4 lines
    lines = lines[4:]
    df = []
    for line in lines:
        ncdcid = line[0:8]
        wban = line[9:14]
        coopid = line[15:21]
        call = line[22:26]
        name = line[27:57]
        alt_name = line[58:88]
        country = line[89:109]
        st = line[110:112]
        county = line[113:143]
        lat = line[144:153]
        lon = line[154:164]
        elev = line[165:171]
        utc = line[172:177]
        stntype = line[178:228]
        begdt = line[229:237]
        ghcnd = line[238:249]
        elev_p = line[250:256]
        elev_a = line[257:263]

        row = [
            ncdcid,
            wban,
            coopid,
            call,
            name,
            alt_name,
            country,
            st,
            county,
            lat,
            lon,
            elev,
            utc,
            stntype,
            begdt,
            ghcnd,
            elev_p,
            elev_a,
        ]
        row = [x.decode("utf-8") for x in row]  # Convert to string
        row = [x.strip() for x in row]  # Strip whitespace
        df.append(row)

    # # Convert to pandas
    stationsasos = pd.DataFrame(
        df,
        columns=[
            "NCDCID",
            "WBAN",
            "COOPID",
            "CALL",
            "NAME",
            "ALTNAME",
            "COUNTRY",
            "ST",
            "COUNTY",
            "LAT",
            "LON",
            "ELEV",
            "UTC",
            "STNTYPE",
            "STARTDATE",
            "GHCN-DailyID",
            "Barometer_elev",
            "Anemometer_elev",
        ],
    )

    # Now, merge the two dataframes
    asosawosstations = pd.concat([stationsasos, stations], axis=0, ignore_index=True)

    # Fill any blank spaces with NaN
    asosawosstations = asosawosstations.replace(r"^\s*$", np.nan, regex=True)

    # Drop 6 records without a WBAN (none in WECC)
    asosawosstations = asosawosstations.dropna(subset=["WBAN"])
    # Note: Only 1 record appears in both ASOS and AWOS - Rock Springs AP in WY.

    return asosawosstations


# Function to get up to date station list of ASOS AWOS stations in WECC.
# Pulls in ISD station list and ASOSAWOS station list (two separate csvs), joins by WBAN and returns list of station IDs.
# Inputs: path to terrestrial WECC shapefile, path to marine WECC file.
# Both paths given relative to home directory for git project.
# Outputs: saves 2 sets of metadata to AWS, returns filtered ISD station object for use in get_asosawos_data_ftp().
def get_wecc_stations(
    terrpath, marpath
):  # Could alter script to have shapefile as input also, if there's a use for this.
    ## Login.
    ## using ftplib, get list of stations as csv
    filename = "isd-history.csv"
    ftp = FTP("ftp.ncdc.noaa.gov")
    ftp.login()  # user anonymous, password anonymous
    ftp.cwd("pub/data/noaa/")  # Change WD.

    # Read in ISD stations.
    r = BytesIO()
    ftp.retrbinary("RETR " + filename, r.write)
    r.seek(0)

    ## Read in csv and only filter to include US stations.
    stations = pd.read_csv(r)
    weccstations = stations[(stations["CTRY"] == "US")]

    # Use spatial geometry to only keep points in wecc marine / terrestrial areas.
    geometry = [
        Point(xy) for xy in zip(weccstations["LON"], weccstations["LAT"])
    ]  # Zip lat lon coords.
    weccgeo = gp.GeoDataFrame(
        weccstations, crs="EPSG:4326", geometry=geometry
    )  # Convert to geodataframe.

    ## get bbox of WECC to use to filter stations against
    t, m, bbox = calc_pull.get_wecc_poly(terrpath, marpath)  # Call get_wecc_poly.

    # Get terrestrial stations.
    weccgeo = weccgeo.to_crs(t.crs)  # Convert to CRS of terrestrial stations.
    terwecc = sjoin(
        weccgeo.dropna(), t, how="left"
    )  # Only keep stations in terrestrial WECC region.
    terwecc = terwecc.dropna()  # Drop empty rows.

    # Get marine stations.
    marwecc = sjoin(
        weccgeo.dropna(), m, how="left"
    )  # Only keep stations in marine WECC region.
    marwecc = marwecc.dropna()  # Drop empty rows.

    # Join and remove duplicates using USAF and WBAN as combined unique identifier.
    weccstations = pd.concat(
        [terwecc.iloc[:, :11], marwecc.iloc[:, :11]], ignore_index=True, sort=False
    ).drop_duplicates(["USAF", "WBAN"], keep="first")

    # Generate ID from USAF/WBAN combo for API call. This follows the naming convention used by FTP/AWS for file names.
    # Add leading zeros where they are missing from WBAN stations.
    weccstations["ISD-ID"] = (
        weccstations["USAF"]
        + "-"
        + weccstations["WBAN"].astype("str").str.pad(5, side="left", fillchar="0")
    )

    # Reformat time strings for FTP/API call.
    weccstations["start_time"] = [
        datetime.strptime(str(i), "%Y%m%d").strftime("%Y-%m-%d")
        for i in weccstations["BEGIN"]
    ]
    weccstations["end_time"] = [
        datetime.strptime(str(i), "%Y%m%d").strftime("%Y-%m-%d")
        for i in weccstations["END"]
    ]

    # Now, read in ASOS and AWOS station files and use to filter to only keep ASOS/AWOS stations.
    asosawos = get_asosawos_stations()

    # Make columns have the same format
    asosawos["WBAN"] = asosawos["WBAN"].astype(str).str.pad(5, fillchar="0")
    weccstations["WBAN"] = weccstations["WBAN"].astype(str).str.pad(5, fillchar="0")

    # Convert ASOSAWOS elevation to feet
    asosawos["ELEV"] = asosawos["ELEV"].astype(float) * 0.3048

    # Sort both dataframes by start date
    asosawos = asosawos.sort_values("STARTDATE")
    weccstations = weccstations.sort_values("BEGIN")

    # Only keep ASOS-AWOS stations in WECC.
    # Merging creates duplicates but gives the number of stations expected.
    # For this stage, keep the 2 sets of station metadata separate and just filter using WBAN.
    m1 = weccstations.WBAN.isin(asosawos.WBAN)
    m2 = asosawos.WBAN.isin(weccstations.WBAN)

    weccstations = weccstations[m1]
    asosawos = asosawos[m2]

    weccstations.reset_index(inplace=True, drop=True)
    asosawos.reset_index(inplace=True, drop=True)

    # Write ASOS AWOS station list to CSV.
    csv_buffer = StringIO()
    asosawos.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3_cl.put_object(
        Bucket=bucket_name, Body=content, Key=directory + "stationlist_asosawos.csv"
    )

    # Write filtered ISD station list to CSV.
    csv_buffer = StringIO()
    weccstations.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3_cl.put_object(
        Bucket=bucket_name, Body=content, Key=directory + "stationlist_isd_asosawos.csv"
    )

    return weccstations


def asosawos_retry_downloads(token, bucket_name, network):
    # Get list of files in folder
    prefix = "1_raw_wx/" + network
    files = []
    for item in s3.Bucket(bucket_name).objects.filter(Prefix=prefix):
        file = str(item.key)
        files += [file]
    files = list(filter(lambda f: f.endswith(".gz"), files))  # Get list of file names

    files = [file for file in files if "errors" not in file]
    files = [
        file for file in files if "station" not in file
    ]  # Remove error and station list files

    # Get only station IDs from file names
    stations = [file.split("/")[-1] for file in files]
    stations = [file.replace(".gz", "") for file in stations]
    stations = [file[0:-5] for file in stations]
    print(stations)

    # Read in station list
    station_list = s3_cl.get_object(
        Bucket=bucket_name, Key="1_raw_wx/ASOSAWOS/stationlist_isd_asosawos.csv"
    )
    station_list = pd.read_csv(station_list["Body"])

    print(station_list)
    # Get list of IDs not in download folder
    missed_stations = [id for id in station_list["ISD-ID"] if id not in stations]
    missed_ids = station_list[
        ["ISD-ID", "start_time"]
    ]  # Format list in way that MADIS_pull script wants it.
    missed_ids = missed_ids[missed_ids["ISD-ID"].isin(missed_stations)]
    downloaded_ids = station_list[~station_list["ISD-ID"].isin(missed_stations)]
    print(missed_ids)
    print(downloaded_ids["ISD-ID"].tail(20))


# stations = get_wecc_stations(wecc_terr, wecc_mar)
asosawos_retry_downloads(
    token=config.token, bucket_name="wecc-historical-wx", network="ASOSAWOS"
)
# Compare file list to stations
