"""
OtherISD_clean.py

This script performs data cleaning for all non-ASOS/AWOS ISD networks for ingestion into the Historical Observations Platform.

Approach
--------
(1) Read through variables, and calculates derived priority variables if not observed
(2) Drops unnecessary variables
(3) Converts station metadata to standard format, with unique identifier
(4) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(5) Converts missing data to standard format
(6) Tracks existing qa/qc flag for review
(7) Merge files by station, and outputs cleaned variables as a single .nc file for an individual network.

Functions
---------
- clean_otherisd: Clean non-ASOSAWOS ISD stations.

Intended Use
------------
Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.

References
----------
https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf

Notes
-----
Removed variables are not saved as removedvars.csv, due to the large number of variables available. These are documented in the ISD format document.
See ISD format document for list of available variables. The QA/QC flag dictionary has been manually formatted and uploaded to the QAQC folder for ASOS/AWOS data.
"""

import os
import xarray as xr
from datetime import datetime, date, timedelta
import re
import numpy as np
import pandas as pd
import boto3
from io import BytesIO, StringIO
import random
from ftplib import FTP
import gzip
import math
import csv
import traceback
import warnings

# Optional: Silence pandas' future warnings about regex (not relevant here)
warnings.filterwarnings(action="ignore", category=FutureWarning)

from clean_utils import get_file_paths
import calc_clean

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes
BUCKET_NAME = "wecc-historical-wx"
WECC_TERR = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
WECC_MAR = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


# Set up directory to save files, if it doesn't already exist.
try:
    os.mkdir("temp")
except:
    pass


def clean_otherisd(rawdir: str, cleandir: str):
    """
    Clean non-ASOSAWOS ISD stations.

    Parameters
    ----------
    rawdir : str
        path to raw data bucket
    cleandir : str
        path to cleaned data bucket

    Returns
    -------
    None
    """
    network = "OtherISD"

    # Set up error handling.
    errors = {"File": [], "Time": [], "Error": []}
    # Set end time to be current time at beginning of download: for error handling csv.
    end_api = datetime.now().strftime("%Y%m%d%H%M")
    # For attributes of netCDF file.
    timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    try:
        # Get files
        files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=rawdir):
            file = str(item.key)
            files += [file]

        # Set up lat/lon bounds for filtering data
        try:
            t, m, bbox = calc_clean.get_wecc_poly(WECC_TERR, WECC_MAR)
            lonmin, lonmax = float(bbox["minx"]), float(bbox["maxx"])
            latmin, latmax = float(bbox["miny"]), float(bbox["maxy"])
        except:
            # If geospatial call fails, hardcode.
            lonmin, lonmax = -139.047795, -102.03721
            latmin, latmax = 30.142739, 60.003861

        # Get station file and read in metadata.
        station_file = [file for file in files if "stationlist" in file]
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj["Body"].read()))
        stations = station_file["ISD-ID"].dropna().astype(str)

        # Remove error, station files
        files = [file for file in files if ".gz" in file]
        files = [file for file in files if "stationlist" not in file]
        files = [file for file in files if "error" not in file]

    except Exception as e:
        # If unable to read files from rawdir, break function.
        print(e)
        errors["File"].append("Whole network")
        errors["Time"].append(end_api)
        errors["Error"].append(f"Whole network error: {e}")

    else:
        # Use ID to grab all files linked to station.
        for id in stations:
            subfiles = list(filter(lambda f: id in f, files))
            # Sort files by year in order to concatenate in correct order.
            subfiles = sorted(subfiles)
            file_count = len(subfiles)

            station = "OtherISD_" + id.replace("-", "")
            station_metadata = station_file.loc[station_file["ISD-ID"] == id]

            if bool(subfiles) is False:
                # If ID returns no files, go to next ID.
                print(f"No raw data found for {station} on AWS.")
                errors["File"].append(station)
                errors["Time"].append(end_api)
                errors["Error"].append("No raw data found for this station on AWS.")
                continue

            # Initialize list of dictionaries.
            data = {
                "station": [],
                "time": [],
                "lat": [],
                "lon": [],
                "elevation": [],
                "qaqc_process": [],
                "ps": [],
                "ps_qc": [],
                "ps_altimeter": [],
                "ps_altimeter_qc": [],
                "psl": [],
                "psl_qc": [],
                "tas": [],
                "tas_qc": [],
                "tdps": [],
                "tdps_qc": [],
                "pr": [],
                "pr_qc": [],
                "pr_duration": [],
                "pr_depth_qc": [],
                "hurs": [],
                "hurs_qc": [],
                "hurs_flag": [],
                "hurs_duration": [],
                "hurs_temp": [],
                "hurs_temp_qc": [],
                "hurs_temp_flag": [],
                "rsds": [],
                "rsds_duration": [],
                "rsds_qc": [],
                "rsds_flag": [],
                "sfcWind": [],
                "sfcWind_qc": [],
                "sfcWind_dir": [],
                "sfcWind_method": [],
                "sfcWind_dir_qc": [],
            }

            for file in subfiles:
                obj = s3.Object(BUCKET_NAME, file)
                try:
                    with gzip.GzipFile(fileobj=obj.get()["Body"]) as gzipped_csv_file:
                        csv_reader = csv.reader(
                            StringIO(gzipped_csv_file.read().decode())
                        )
                        for row in csv_reader:
                            # Each row is a record
                            # Initialize all variables and set to be NA by default.
                            string = row[0]

                            # Filter station by latitude and longitude. Only keep obs in WECC.
                            # POS 29-34: GEOPHYSICAL-POINT-OBSERVATION latitude coordinate
                            lat = float(string[28:34]) / 1000  # Unit: degree
                            # POS 35-41: GEOPHYSICAL-POINT-OBSERVATION longitude coordinate
                            lon = float(string[34:41]) / 1000  # Unit: degree
                            if (
                                lat > latmax
                                or lat < latmin
                                or lon > lonmax
                                or lon < lonmin
                            ):
                                errors["File"].append(station)
                                errors["Time"].append(end_api)
                                # some years hav valid lat/lon, others do not
                                errors["Error"].append(f"File not in WECC: {file}")
                                break

                            else:
                                # POS 16-23: GEOPHYSICAL-POINT-OBSERVATION date, # POS 24-27: GEOPHYSICAL-POINT-OBSERVATION time (in UTC)
                                time = datetime.strptime(string[15:27], "%Y%m%d%H%M")

                                # POS 47-51: GEOPHYSICAL-POINT-OBSERVATION elevation dimension
                                elevation = float(string[46:51])
                                # POS 57-60: METEOROLOGICAL-POINT-OBSERVATION quality control process name
                                qaqc_process = string[56:60]
                                # V01 = No A or M Quality Control applied
                                # V02 = Automated Quality Control
                                # V03 = subjected to Quality Control

                                # Mandatory data

                                # POS: 61-63: WIND-OBSERVATION direction angle
                                sfcWind_dir = int(string[60:63])  # Units degrees

                                # POS: 64-64: WIND-OBSERVATION direction quality code
                                sfcWind_dir_qc = string[63]

                                # POS: 65-65 WIND-OBSERVATION type code
                                sfcWind_method = string[64]

                                # POS: 66-69: WIND-OBSERVATION speed rate
                                sfcWind = float(string[65:69]) / 10  # Units m/s

                                # POS: 70-70: WIND-OBSERVATION speed quality code
                                sfcWind_qc = string[69]
                                # Note: One row returns A here, this is an error.

                                # POS 88-92: AIR-TEMPERATURE-OBSERVATION air temperature
                                tas = float(string[87:92]) / 10

                                # POS 93: AIR-TEMPERATURE-OBSERVATION air temperature quality code
                                tas_qc = string[92]

                                # POS 94-98: AIR-TEMPERATURE-OBSERVATION dew point temperature
                                tdps = float(string[93:98]) / 10

                                # POS 99-99: AIR-TEMPERATURE-OBSERVATION dew point quality code
                                tdps_qc = string[98]

                                # POS 100-104: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure
                                # In hectopascals, CONVERT.
                                psl = float(string[99:104]) / 10

                                # POS 105-105: ATMOSPHERIC-PRESSURE-OBSERVATION sea level pressure quality code
                                psl_qc = string[104]

                                # Additional data - begins with ADD

                                # Liquid precipitation - begins with AA[1-4]

                                # Figure out length of precip. string.
                                precip_len = re.search(
                                    "(?<=AA1|AA2|AA3|AA4)[\da-zA-Z]{8}", string
                                )
                                if precip_len is not None:
                                    # If precip data exists
                                    # Get everything after AA#
                                    # QA/QC codes include letters. So, use precip_length to pull string of correct length.
                                    precip = precip_len.group()
                                    # Get the 8 strings after AA1/2/3/4
                                    pr_duration = int(precip[0:2])  # Hours
                                    pr = float(precip[2:6]) / 10  # In mm
                                    # QA for precipitation depth
                                    pr_depth_qc = int(precip[6:7])
                                    # QA for precipitation data.
                                    pr_qc = precip[7:8]

                                    # Take first measurement as primary, unless missing.
                                    if float(precip[2:6]) == 9999:
                                        # If precip depth of first report is missing.
                                        precip = re.search(
                                            "(?<=AA1|AA2|AA3|AA4)[\da-zA-Z]{16}", string
                                        )
                                        if precip is not None:
                                            precip = precip.group()
                                            # Get second set of precip values.
                                            if precip[9:10].isnumeric():
                                                print(precip)
                                                pr_duration = int(precip[8:10])  # Hours
                                                pr = float(precip[10:14]) / 10  # In mm
                                                # QA for precipitation depth
                                                pr_depth_qc = int(precip[14:15])
                                                # QA for precipitation data.
                                                pr_qc = precip[15:16]
                                        else:
                                            pr = np.nan
                                            pr_qc = np.nan
                                            pr_duration = np.nan
                                            pr_depth_qc = np.nan

                                else:
                                    pr = np.nan
                                    pr_qc = np.nan
                                    pr_duration = np.nan
                                    pr_depth_qc = np.nan

                                # Relative humidity - shouldn't be in this dataset, but capture if it is.
                                # Section starts with CH
                                hurs_string = re.search(
                                    "(?<=CH1|CH2)[\da-zA-Z]{15}", string
                                )

                                if hurs_string is not None:
                                    # If hurs data exists
                                    hurs_string = hurs_string.group()

                                    hurs_duration = int(hurs_string[0:2])  # Minutes
                                    hurs_temp = int(hurs_string[2:7]) / 10  # In deg. C
                                    hurs_temp_qc = hurs_string[7]
                                    hurs_temp_flag = int(hurs_string[8])
                                    hurs = int(hurs_string[9:13]) / 10  # In percent
                                    hurs_qc = hurs_string[13]
                                    hurs_flag = int(hurs_string[14])

                                else:
                                    hurs = np.nan
                                    hurs_qc = np.nan
                                    hurs_flag = np.nan
                                    hurs_duration = np.nan
                                    hurs_temp = np.nan
                                    hurs_temp_qc = np.nan
                                    hurs_temp_flag = np.nan

                                # Solar radiation - "global irradiance" in our dataset.
                                # Not keeping "direct beam irradiance", "diffuse irradiance" which are also provided in same string.
                                # Section starts with CH
                                rsds_string = re.search(
                                    "(?<=GM1)[\da-zA-Z]{11}", string
                                )

                                if rsds_string is not None:
                                    # If rsds data exists
                                    rsds_string = rsds_string.group()
                                    # Time period over which solar radiation integrated, in minutes
                                    rsds_duration = float(rsds_string[0:4])
                                    rsds = float(rsds_string[4:8])  # In w/m2
                                    rsds_flag = rsds_string[8:10]
                                    rsds_qc = rsds_string[10:12]

                                else:
                                    rsds = np.nan
                                    rsds_duration = np.nan
                                    rsds_qc = np.nan
                                    rsds_flag = np.nan

                                # Station pressure, Section starts with CH
                                ps_string = re.search("(?<=MA1)[\da-zA-Z]{12}", string)

                                if ps_string is not None:
                                    # If pressure exists
                                    ps_string = ps_string.group()
                                    ps_altimeter = float(ps_string[0:5]) / 10  # HPa
                                    ps_altimeter_qc = ps_string[5]
                                    ps = float(ps_string[6:11]) / 10  # HPa
                                    ps_qc = ps_string[11]

                                else:
                                    ps = np.nan
                                    ps_qc = np.nan
                                    ps_altimeter = np.nan
                                    ps_altimeter_qc = np.nan

                                # Standardize NAs
                                try:
                                    if lat == 99.999 or lon == 99.999:
                                        # If lat or lon values are NA, skip observation.
                                        continue
                                    if elevation == 9999:
                                        elevation = np.nan
                                    if ps == 9999.9:
                                        ps = np.nan
                                    if ps_altimeter == 9999.9:
                                        ps_altimeter = np.nan
                                    if psl == 9999.9:
                                        psl = np.nan
                                    if tas == 999.9:
                                        tas = np.nan
                                    if tdps == 999.9:
                                        tdps = np.nan
                                    if pr == 999.9:
                                        pr = np.nan
                                    if pr_duration == 99:
                                        pr_duration = np.nan
                                    if hurs == 999.9:
                                        hurs = np.nan
                                    if hurs_duration == 99:
                                        hurs_duration = np.nan
                                    if hurs_temp == 999.9:
                                        hurs_temp = np.nan
                                    if rsds == 9999:
                                        rsds = np.nan
                                    if rsds_duration == 9999:
                                        rsds_duration = np.nan
                                    if sfcWind_dir == 999:
                                        sfcWind_dir = np.nan
                                    if sfcWind == 999.9:
                                        sfcWind = np.nan
                                except Exception as e:
                                    # Could add more robust error handling here if desired, but should not be necessary.
                                    print(e)

                                # For each row of data, append data to list.
                                columns = [
                                    "station",
                                    "time",
                                    "lat",
                                    "lon",
                                    "elevation",
                                    "qaqc_process",
                                    "ps",
                                    "ps_qc",
                                    "ps_altimeter",
                                    "ps_altimeter_qc",
                                    "psl",
                                    "psl_qc",
                                    "tas",
                                    "tas_qc",
                                    "tdps",
                                    "tdps_qc",
                                    "pr",
                                    "pr_qc",
                                    "pr_duration",
                                    "pr_depth_qc",
                                    "hurs",
                                    "hurs_qc",
                                    "hurs_flag",
                                    "hurs_duration",
                                    "hurs_temp",
                                    "hurs_temp_qc",
                                    "hurs_temp_flag",
                                    "rsds",
                                    "rsds_duration",
                                    "rsds_qc",
                                    "rsds_flag",
                                    "sfcWind",
                                    "sfcWind_qc",
                                    "sfcWind_dir",
                                    "sfcWind_method",
                                    "sfcWind_dir_qc",
                                ]
                                variables = [
                                    station,
                                    time,
                                    lat,
                                    lon,
                                    elevation,
                                    qaqc_process,
                                    ps,
                                    ps_qc,
                                    ps_altimeter,
                                    ps_altimeter_qc,
                                    psl,
                                    psl_qc,
                                    tas,
                                    tas_qc,
                                    tdps,
                                    tdps_qc,
                                    pr,
                                    pr_qc,
                                    pr_duration,
                                    pr_depth_qc,
                                    hurs,
                                    hurs_qc,
                                    hurs_flag,
                                    hurs_duration,
                                    hurs_temp,
                                    hurs_temp_qc,
                                    hurs_temp_flag,
                                    rsds,
                                    rsds_duration,
                                    rsds_qc,
                                    rsds_flag,
                                    sfcWind,
                                    sfcWind_qc,
                                    sfcWind_dir,
                                    sfcWind_method,
                                    sfcWind_dir_qc,
                                ]

                                for i, j in zip(columns, variables):
                                    data[i].append(j)

                                # For testing: progress update. Print status update every 1k rows.
                                # If station reports every 10 min, expect up to 50k observations per year.
                                if len(data) % 1000 == 0:
                                    print(f"{len(data)} observations parsed.")

                except Exception as e:
                    print(file, e)
                    errors["File"].append(file)
                    errors["Time"].append(end_api)
                    errors["Error"].append(f"Error in pandas df set-up: {e}")

            # Take all lists and convert to dataframe.
            # Remove any variables where all the data is nan.
            df = pd.DataFrame.from_dict(data)

            if df.empty is True:
                print(
                    f"Station {station} reports all nan meteorological data - not cleaned."
                )

            if df.empty is False:
                # If there is data in the dataframe, convert to xarray object.
                try:
                    # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022.
                    df = df.loc[
                        (df["time"] < "2022-09-01") & (df["time"] > "1979-12-31")
                    ]

                    ds = df.to_xarray()

                    # Update global attributes
                    ds = ds.assign_attrs(title="ISD (non-ASOS/AWOS networks) cleaned")
                    ds = ds.assign_attrs(institution="Eagle Rock Analytics / Cal Adapt")
                    ds = ds.assign_attrs(source="")
                    ds = ds.assign_attrs(
                        history=f"OtherISD_clean.py script run on {timestamp} UTC"
                    )
                    ds = ds.assign_attrs(
                        comment="Intermediate data product: may not have been subject to any cleaning or QA/QC processing"
                    )
                    ds = ds.assign_attrs(license="")
                    ds = ds.assign_attrs(citation="")
                    ds = ds.assign_attrs(
                        disclaimer="This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document."
                    )

                    # Add station metadata
                    # Station name
                    ds = ds.assign_attrs(
                        station_name=station_metadata["STATION NAME"].values[0]
                    )

                    # Other station IDs - only add if not NA
                    ids = ["USAF", "WBAN", "ICAO"]
                    for i in ids:
                        if isinstance(
                            station_metadata[i].values[0], (int, str, np.integer)
                        ):
                            ds.attrs[i] = station_metadata[i].values[0]

                    # Sensor heights
                    # Could be cross referenced with HOMR data, but currently unknown.
                    ds = ds.assign_attrs(thermometer_height_m=np.nan)
                    # As per p. 99 of ISD manual.
                    ds = ds.assign_attrs(anemometer_height_m=1.5)
                    # Could be cross referenced with HOMR data, but currently unknown.
                    ds = ds.assign_attrs(barometer_elev_m=np.nan)
                    # Keep count of how many files merged per station.
                    ds = ds.assign_attrs(raw_files_merged=file_count)

                    # Update dimensions and coordinates

                    # Add dimensions: station ID and time.
                    ds = ds.set_coords("time").swap_dims({"index": "time"})
                    ds = ds.assign_coords(id=str(station))
                    # Add station as index.
                    ds = ds.expand_dims("id")
                    # Drop station variable and index coordinate.
                    ds = ds.drop_vars(("station", "index"))
                    # Rename id to station.
                    ds = ds.rename({"id": "station"})

                    # Add coordinates: lat and longitude.
                    ds = ds.set_coords(("lat", "lon"))

                    # Update dimension and coordinate attributes.

                    # Time
                    ds["time"].attrs["long_name"] = "time"
                    ds["time"].attrs["standard_name"] = "time"
                    ds["time"].attrs["comment"] = "In UTC."

                    # Station ID
                    ds["station"].attrs["long_name"] = "station_id"
                    ds["station"].attrs[
                        "comment"
                    ] = "Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network."

                    # lat
                    ds["lat"].attrs["long_name"] = "latitude"
                    ds["lat"].attrs["standard_name"] = "latitude"
                    ds["lat"].attrs["units"] = "degrees_north"

                    # Longitude
                    ds["lon"].attrs["long_name"] = "longitude"
                    ds["lon"].attrs["standard_name"] = "longitude"
                    ds["lon"].attrs["units"] = "degrees_east"

                    # Elevation
                    ds["elevation"].attrs[
                        "standard_name"
                    ] = "height_above_mean_sea_level"
                    ds["elevation"].attrs["long_name"] = "station_elevation"
                    ds["elevation"].attrs["units"] = "meter"
                    # Define which direction is positive
                    ds["elevation"].attrs["positive"] = "up"

                    # Update variable attributes and do unit conversions

                    # tas: air surface temperature (K)
                    if "tas" in ds.keys():
                        ds["tas"] = calc_clean._unit_degC_to_K(ds["tas"])

                        ds["tas"].attrs["long_name"] = "air_temperature"
                        ds["tas"].attrs["standard_name"] = "air_temperature"
                        ds["tas"].attrs["units"] = "degree_Kelvin"
                        # List other variables associated with variable (QA/QC)
                        ds["tas"].attrs["ancillary_variables"] = "tas_qc"
                        ds["tas"].attrs["comment"] = "Converted from Celsius."

                    if "tas_qc" in ds.keys():
                        ds["tas_qc"].attrs[
                            "flag_values"
                        ] = "0 1 2 3 4 5 6 7 9 A C I M P R U"
                        ds["tas_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # ps: surface air pressure (Pa)
                    if "ps" in ds.keys():
                        ds["ps"] = calc_clean._unit_pres_hpa_to_pa(ds["ps"])
                        ds["ps"].attrs["long_name"] = "station_air_pressure"
                        ds["ps"].attrs["standard_name"] = "air_pressure"
                        ds["ps"].attrs["units"] = "Pa"
                        # List other variables associated with variable (QA/QC)
                        ds["ps"].attrs[
                            "ancillary_variables"
                        ] = "ps_qc ps_altimeter ps_altimeter_qc"
                        ds["tas"].attrs["comment"] = "Converted from hPa to Pa."

                        # Delete sea level pressure if station pressure included.
                        if "psl" in ds.keys():
                            ds.drop(("psl", "psl_qc"))

                        if "ps_qc" in ds.keys():
                            ds["ps_qc"].attrs["flag_values"] = "0 1 2 3 4 5 6 7 M 9"
                            ds["ps_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."

                        if "ps_altimeter" in ds.keys():
                            ds["ps_altimeter"] = calc_clean._unit_pres_hpa_to_pa(
                                ds["ps_altimeter"]
                            )
                            ds["ps_altimeter"].attrs["long_name"] = "altimeter_setting"
                            ds["ps_altimeter"].attrs["units"] = "Pa"
                            # List other variables associated with variable (QA/QC)
                            ds["ps_altimeter"].attrs[
                                "ancillary_variables"
                            ] = "ps ps_qc ps_altimeter ps_altimeter_qc"
                            ds["ps_altimeter"].attrs[
                                "comment"
                            ] = "Converted from hPa to Pa. The pressure value to which an aircraft altimeter is set so that it will indicate the altitude relative to mean sea level of an aircraft on the ground at the location for which the value was determined."  # Description of variable meaning, as not CF-standard.

                        if "ps_altimeter_qc" in ds.keys():
                            ds["ps_altimeter_qc"].attrs[
                                "flag_values"
                            ] = "0 1 2 3 4 5 6 7 M 9"
                            ds["ps_altimeter_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."

                    # If station air pressure not reported, convert from sea level pressure.
                    # Inputs: sea level pressure (mb/hPa), elevation (m), and air temperature (K)
                    elif "psl" in ds.keys():
                        ds["psl"] = calc_clean._unit_pres_hpa_to_pa(ds["psl"])
                        ds["psl"].attrs["long_name"] = "sea_level_air_pressure"
                        ds["psl"].attrs["standard_name"] = "air_pressure_at_sea_level"
                        ds["psl"].attrs["units"] = "Pa"
                        # List other variables associated with variable (QA/QC)
                        ds["psl"].attrs["ancillary_variables"] = "psl_qc"
                        ds["psl"].attrs["comment"] = "Converted from hPa to Pa."

                    # tdps: dew point temperature (K)
                    # tdps always provided by ISD reports, so no need to calculate here.
                    if "tdps" in ds.keys():  # If variable already exists, rename.
                        ds["tdps"] = calc_clean._unit_degC_to_K(ds["tdps"])

                        ds["tdps"].attrs["long_name"] = "dew_point_temperature"
                        ds["tdps"].attrs["standard_name"] = "dew_point_temperature"
                        ds["tdps"].attrs["units"] = "degree_Kelvin"
                        # List other variables associated with variable (QA/QC)
                        ds["tdps"].attrs["ancillary_variables"] = "tdps_qc"
                        ds["tdps"].attrs["comment"] = "Converted from Celsius."

                    if "tdps_qc" in ds.keys():
                        ds["tdps_qc"].attrs[
                            "flag_values"
                        ] = "0 1 2 3 4 5 6 7 9 A C I M P R U"
                        ds["tdps_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # pr: precipitation (Flag for discussion about CF compliance and units)
                    if "pr" in ds.keys():
                        # Note leave final conversion here to next stage.
                        ds["pr"].attrs["long_name"] = "precipitation_accumuation"
                        ds["pr"].attrs["units"] = "mm/?"
                        # List other variables associated with variable (QA/QC)
                        ds["pr"].attrs[
                            "ancillary_variables"
                        ] = "pr_qc pr_depth_qc pr_duration"
                        ds["pr"].attrs["comment"] = ""  # To be completed.

                        ds["pr_duration"].attrs[
                            "long_name"
                        ] = "precipitation measurement interval"
                        ds["pr_duration"].attrs["units"] = "hours"
                        # List other variables associated with variable (QA/QC)
                        ds["pr_duration"].attrs[
                            "ancillary_variables"
                        ] = "pr pr_qc pr_depth_qc"

                    if "pr_qc" in ds.keys():
                        ds["pr_qc"].attrs[
                            "flag_values"
                        ] = "0 1 2 3 4 5 6 7 9 A I M P R U"
                        ds["pr_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        ds["pr_depth_qc"].attrs[
                            "flag_values"
                        ] = "1 2 3 4 5 6 7 8 E I J 9"
                        ds["pr_depth_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # hurs: relative humidity
                    if "hurs" in ds.keys():
                        ds["hurs"].attrs["long_name"] = "average_relative_humidity"
                        ds["hurs"].attrs["standard_name"] = "relative_humidity"
                        ds["hurs"].attrs["units"] = "percent"
                        # List other variables associated with variable (QA/QC)
                        ds["hurs"].attrs[
                            "ancillary_variables"
                        ] = "hurs_qc hurs_flag hurs_duration hurs_temp hurs_temp_qc hurs_temp_flag"

                        ds["hurs_duration"].attrs[
                            "long_name"
                        ] = "relative_humidity_duration"
                        ds["hurs_duration"].attrs["units"] = "minutes"
                        # List other variables associated with variable (QA/QC)
                        ds["hurs_duration"].attrs[
                            "ancillary_variables"
                        ] = "hurs hurs_qc hurs_flag hurs_temp hurs_temp_qc hurs_temp_flag"

                        ds["hurs_temp"].attrs[
                            "long_name"
                        ] = "average_air_temperature_at_hurs_instrument"
                        ds["hurs_temp"].attrs["units"] = "degree_celsius"
                        # List other variables associated with variable (QA/QC)
                        ds["hurs_temp"].attrs[
                            "ancillary_variables"
                        ] = "hurs hurs_qc hurs_flag hurs_temp hurs_temp_qc hurs_temp_flag"

                    if "hurs_qc" in ds.keys():
                        ds["hurs_qc"].attrs["flag_values"] = "1 3 9"
                        ds["hurs_qc"].attrs[
                            "flag_meanings"
                        ] = "passed_all_qc_checks failed_all_qc_checks missing"

                        ds["hurs_flag"].attrs["flag_values"] = "0 1 2 3 4 5 6 7 8 9"
                        ds["hurs_flag"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                        ds["hurs_temp_qc"].attrs["flag_values"] = "1 3 9"
                        ds["hurs_temp_qc"].attrs[
                            "flag_meanings"
                        ] = "passed_all_qc_checks failed_all_qc_checks missing"

                        ds["hurs_temp_flag"].attrs[
                            "flag_values"
                        ] = "0 1 2 3 4 5 6 7 8 9"
                        ds["hurs_temp_flag"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
                    if "rsds" in ds.keys():
                        ds["rsds"].attrs["long_name"] = "solar_radiation"
                        ds["rsds"].attrs[
                            "standard_name"
                        ] = "surface_downwelling_shortwave_flux_in_air"
                        ds["rsds"].attrs["units"] = "W m-2"
                        # List other variables associated with variable (QA/QC)
                        ds["rsds"].attrs[
                            "ancillary_variables"
                        ] = "rsds_duration rsds_qc rsds_flag"
                        ds["rsds"].attrs[
                            "comment"
                        ] = "Waveband ranges from 0.4 - 2.3 micrometers."

                        ds["rsds_duration"].attrs["long_name"] = "time_period"
                        ds["rsds_duration"].attrs["units"] = "minutes"
                        # List other variables associated with variable (QA/QC)
                        ds["rsds_duration"].attrs[
                            "ancillary_variables"
                        ] = "rsds rsds_qc rsds_flag"
                        ds["rsds_duration"].attrs[
                            "comment"
                        ] = "Time period over which solar radiation is integrated."

                        if "rsds_qc" in ds.keys():
                            ds["rsds_qc"].attrs["flag_values"] = "0 1 2 3 9"
                            ds["rsds_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."

                            ds["rsds_flag"].attrs[
                                "flag_values"
                            ] = "00 01 02 03 04 05 06 07 08 09 10-93 94-97 98 99"
                            ds["rsds_flag"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."

                    # sfcWind : wind speed (m/s) (Method of calculation may vary, standardize during hourly merge or QA/QC process.)
                    if "sfcWind" in ds.keys():
                        # No conversions needed, do not add raw column.
                        ds["sfcWind"].attrs["long_name"] = "wind_speed"
                        ds["sfcWind"].attrs["standard_name"] = "wind_speed"
                        ds["sfcWind"].attrs["units"] = "m s-1"
                        # List other variables associated with variable (QA/QC)
                        ds["sfcWind"].attrs[
                            "ancillary_variables"
                        ] = "sfcWind_qc sfcWind_method"
                        ds["sfcWind"].attrs[
                            "comment"
                        ] = "Method of wind speed calculation varies, see sfcWind_method."

                    if "sfcWind_method" in ds.keys():
                        ds["sfcWind_method"].attrs[
                            "long_name"
                        ] = "wind_speed_calculation_method"
                        ds["sfcWind_method"].attrs[
                            "flag_values"
                        ] = "A B C H N R Q T V 9"
                        ds["sfcWind_method"].attrs[
                            "flag_meanings"
                        ] = "abridged_beaufort beaufort calm 5-minute_average normal 60-minute_average squall 180-minute_average variable missing"

                    if "sfcWind_qc" in ds.keys():
                        ds["sfcWind_method"].attrs["flag_values"] = "0 1 2 3 4 5 6 7 9"
                        ds["sfcWind_method"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # sfcWind_dir: wind direction
                    if "sfcWind_dir" in ds.keys():
                        # No conversions needed, do not make raw column.
                        ds["sfcWind_dir"].attrs["long_name"] = "wind_direction"
                        ds["sfcWind_dir"].attrs["standard_name"] = "wind_from_direction"
                        ds["sfcWind_dir"].attrs[
                            "units"
                        ] = "degrees_clockwise_from_north"
                        # List other variables associated with variable (QA/QC)
                        ds["sfcWind_dir"].attrs[
                            "ancillary_variables"
                        ] = "sfcWind_dir_qc"

                    if "sfcWind_dir_qc" in ds.keys():
                        ds["sfcWind_dir_qc"].attrs["flag_values"] = "0 1 2 3 4 5 6 7 9"
                        ds["sfcWind_dir_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # Update attributes for any non-standard variables

                    # QAQC process
                    if "qaqc_process" in ds.keys():
                        ds["qaqc_process"].attrs["long_name"] = "qaqc_process_type"
                        ds["qaqc_process"].attrs["flag_values"] = "V01 V02 V03"
                        ds["qaqc_process"].attrs[
                            "flag_meanings"
                        ] = "no_qaqc automated_qaqc subjected_to_qaqc"

                    # Data source
                    if "data_source" in ds.keys():
                        ds["data_source"].attrs["long_name"] = "source_of_data"
                        ds["data_source"].attrs[
                            "flag_values"
                        ] = "1 2 3 4 5 6 7 8 A B C D E F G H I J K L M N O 9"
                        ds["data_source"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."

                    # For QA/QC flags, replace np.nan with "nan" to avoid h5netcdf overwrite to blank.
                    for key in ds.keys():
                        if "qc" in key:
                            # Coerce all values in key to string
                            ds[key] = ds[key].astype(str)

                    # drop any column that does not have any valid (non-nan data)
                    # need to keep elevation separate, as it does have "valid" nan value, only drop if all other variables are also nans
                    for key in ds.keys():
                        try:
                            if key != "elevation":
                                if np.isnan(ds[key].values).all():
                                    print(f"Dropping empty var: {key}")
                                    ds = ds.drop(key)

                            # only drop elevation if all other variables are also nans
                            if (key == "elevation") & (len(ds.keys()) == 1):
                                # only elevation remains
                                print(f"Dropping empty var: {key}")
                                # slightly unnecessary since the entire dataset will be empty too
                                ds = ds.drop(key)
                                continue

                        except Exception as e:
                            # Add to handle errors for unsupported data types
                            continue

                    # Reorder variables
                    # In following order:
                    desired_order = [
                        "ps",
                        "tas",
                        "tdps",
                        "pr",
                        "hurs",
                        "rsds",
                        "sfcWind",
                        "sfcWind_dir",
                    ]
                    # Only keep vars which are in ds.
                    desired_order = [i for i in desired_order if i in list(ds.keys())]
                    # Retain rest of variables at the bottom.
                    rest_of_vars = [
                        i for i in list(ds.keys()) if i not in desired_order
                    ]
                    new_index = desired_order + rest_of_vars
                    ds = ds[new_index]

                except Exception as e:
                    # If error in xarray reorganization
                    print(file, e)
                    traceback.print_exc()
                    errors["File"].append(station)
                    errors["Time"].append(end_api)
                    errors["Error"].append(f"Error in ds set-up: {e}")

                # Write station file to netcdf.
                if ds is None:
                    # Should be caught be error handling above, but add in case.
                    print(
                        f"Station {file} does not report any valid meteorological data - not cleaned."
                    )
                    errors["File"].append(file)
                    errors["Time"].append(end_api)
                    errors["Error"].append(
                        "Station reports all nan meteorological data."
                    )
                    continue

                else:
                    try:
                        filename = station + ".nc"
                        filepath = cleandir + filename

                        # Write locally
                        ds.to_netcdf(path="temp/temp.nc", engine="netcdf4")

                        # Push file to AWS with correct file name.
                        s3.Bucket(BUCKET_NAME).upload_file("temp/temp.nc", filepath)

                        print(f"Saving {filename} with dims {ds.dims}")
                        ds.close()
                        continue

                    except Exception as e:
                        print(filename, e)
                        errors["File"].append(filename)
                        errors["Time"].append(end_api)
                        errors["Error"].append(
                            f"Error saving ds as .nc file to AWS bucket: {e}"
                        )
                        continue

    # Write errors to csv
    finally:
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=cleandir + f"errors_otherisd_{end_api}.csv",
        )

    return None


# Run function
if __name__ == "__main__":
    rawdir, cleandir, qaqcdir = get_file_paths("OtherISD")
    clean_otherisd(rawdir, cleandir)
