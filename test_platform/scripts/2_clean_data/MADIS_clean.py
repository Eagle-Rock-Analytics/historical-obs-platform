"""
MADIS_clean.py

This script performs data cleaning for networks pulled through Synoptic API for ingestion into the Historical Observations Platform.

Approach
--------
(1) Read through variables and drop unnecessary variables
(2) Converts station metadata to standard format, with unique identifier
(3) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(4) Converts missing data to standard format
(5) Tracks existing qa/qc flag for review
(6) Merge files by station, and outputs cleaned variables as a single .nc file for each station in an individual network.

Functions
---------
- get_qaqc_flags: Get MADIS QA/QC raw flag data and parse into a csv for reference.
- parse_madis_headers: Parsing the header of MADIS csv files.
- parse_madis_to_pandas: Take csv headers and clean MADIS data.
- clean_madis: Cleans MADIS data for a variety of networks.

Intended Use
------------
Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.

References
----------
https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
"""

import os
import xarray as xr
from datetime import datetime, date
import re
import numpy as np
import pandas as pd
import requests
from collections import Counter
import boto3
from io import BytesIO, StringIO
import smart_open
import traceback
import botocore
import warnings

# Optional: Silence pandas' future warnings about regex (not relevant here)
warnings.filterwarnings(action="ignore", category=FutureWarning)

try:
    import config  # Import API keys.
except:
    print("Missing config.py file with API token. Make file if necessary.")
    exit()

from clean_utils import get_file_paths
import calc_clean

s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes
BUCKET_NAME = "wecc-historical-wx"

# Set up directory to save files temporarily, if it doesn't already exist.
try:
    os.mkdir("temp")
except:
    pass


def get_qaqc_flags(token: str, qaqcdir: str, network: str) -> pd.DataFrame:
    """
    Get MADIS QA/QC raw flag data and parse into a csv for reference.

    Parameters
    ----------
    token : str
        API url for QAQC metadata
    qaqcdir : str
        path to AWS qaqc bucket
    network : str
        name of network

    Returns
    -------
    ids : pd.DataFrame
        qaqc flags from MADIS networks

    Notes
    -----
    1. We are no longer retrieving MADIS data from Synoptic. Will need a reconfig.
    """

    url = f"https://api.synopticdata.com/v2/qctypes?token={token}"
    request = requests.get(url).json()

    ids = []
    for each in request["QCTYPES"]:
        ids.append([each["ID"], each["SHORTNAME"], each["NAME"]])
    # Sort by start date (note some stations return 'None' here)
    ids = pd.DataFrame(ids, columns=["FLAG", "SHORTNAME", "NAME"])

    # Save to AWS bucket.
    print("Saving QAQC flags to csv file.")

    csv_buffer = StringIO()
    ids.to_csv(csv_buffer, index=False)
    content = csv_buffer.getvalue()
    s3_cl.put_object(
        Bucket=BUCKET_NAME, Body=content, Key=qaqcdir + f"qaqcdict_{network}.csv"
    )
    return ids


def parse_madis_headers(file: str) -> dict[str, str]:
    """
    Parsing the header of MADIS csv files.

    Parameters
    ----------
    file : str
        file name

    Returns
    -------
    headers : dict of str
        dictionary of relevant station attributes for metadata from MADIS raw data header
    """

    url = f"s3://{BUCKET_NAME}/{file}"
    index = 0
    unitstocheck = [
        "Pascals",
        "%",
        "m/s",
        "Celsius",
        "QC_type",
        "Degrees",
        "Millimeters",
    ]

    # If no data available for station, both of these variables will remain as NaN
    units = np.nan
    first_row = np.nan
    for line in smart_open.open(url, mode="rb"):
        if index > 10:
            continue
        index += 1

        # row = (line.decode('utf-8'))
        # fix for non-ASCII character, safe for all stations
        row = line.decode(errors="ignore")

        if "STATION:" in row:
            # Skip first row.
            station_id = row.partition(": ")[2].replace(" ", "").replace("\n", "")
            continue

        elif "STATION NAME:" in row:
            station_name = str(
                row.partition(": ")[2].replace("']", "").replace("\n", "")
            )
            station_name = station_name.replace(")", "")
            continue

        elif "LATITUDE" in row:
            latitude = float(row.partition(": ")[2].replace("']", ""))
            continue

        elif "LONGITUDE" in row:
            longitude = float(row.partition(": ")[2].replace("']", ""))
            continue

        elif "ELEVATION" in row:
            if row.partition(": ")[2].replace("\n", "") == "None":
                elevation = np.nan
            else:
                elevation = float(row.partition(": ")[2].replace("']", ""))  # In feet.
            continue

        elif "STATE" in row:
            state = row.partition(": ")[2].replace("']", "").replace("\n", "")
            continue

        elif "Station_ID" in row:
            columns = row
            columns = columns.split(",")
            # Clean: Remove line break from list
            columns = [k.replace("\n", "") for k in columns]
            continue

        elif any(unit in row for unit in unitstocheck):
            units = row
            units = units.replace("\n", "")
            continue

        else:
            # Get row where data begins (station ID first appears)
            if station_id in row:
                first_row = index
                break
            else:
                continue

    # Read in entire csv starting at first_row, using pandas.
    # First, some error handling.

    # Check for duplicated columns. If duplicated names, manually rename as 1 and 2, and then check for duplicates after reading in.
    if set([x for x in columns if columns.count(x) > 1]):
        dup = True
        # Add dup flag.
        dup_col = set([x for x in columns if columns.count(x) > 1])
        # Manually rename second iteration of column.
        d = {
            a: list(range(1, b + 1)) if b > 1 else ""
            for a, b in Counter(columns).items()
        }
        columns = [i + str(d[i].pop(0)) if len(d[i]) else i for i in columns]
    else:
        dup = False
        dup_col = []

    headers = {
        "station_id": station_id,
        "station_name": station_name,
        "latitude": latitude,
        "longitude": longitude,
        "elevation": elevation,
        "state": state,
        "columns": columns,
        "units": units,
        "first_row": first_row,
        "dup": dup,
        "dup_col": dup_col,
    }

    return headers


def parse_madis_to_pandas(
    file: str, headers: dict[str, str], errors: dict[str, str], removedvars: list[str]
) -> pd.DataFrame | None:
    """
    Take csv headers and clean MADIS data.

    Parameters
    ----------
    file : str
        file path from rawdir
    headers : dict of str
        dictionary of relevant station attributes for metadata from MADIS raw data header
    errors : dict of str
        errors dictionary for tracking
    removedvars : list of str
        non-essential variables available in the raw data that are removed

    Returns
    -------
    If success, pd.DataFrame
    If raw file is empty, does not report data within v1 period, or has a mismatched metadata header, return None

    Note
    ----
    Headers is created in parse_madis_headers function for a list of keys, please go check the headers dictionary at the end of that function.
    Errors is created in parse_madis_headers function, please check how it is utilized in the parse_madis_to_pandas function.
    """

    ### TEMPORARY MISMATCH METADATA FIX
    if len(file) > 24:
        # identifies files whose name is longer than 24 characters; for CWOP standard is 23 chars, update pull files are 32
        # skips any file whose name is longer than 24 characters --> targets the CWOP_STID-20YY-MM-DD files
        return None
    ###

    # Read in CSV, removing header.
    obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=file)
    try:
        # Ignore dtype warning here, we resolve this manually below.
        df = pd.read_csv(
            BytesIO(obj["Body"].read()),
            names=headers["columns"],
            header=headers["first_row"] - 1,
            low_memory=False,
        )

    except:
        print(
            f"Raw data file {file} is either empty or does not report data within v1 period (1/1980 - 8/2022). Not cleaned."
        )
        return None

    if len(df.index) == 0:
        print(f"Raw data file for {file} is empty, station not cleaned.")
        return None

    # Handling for timeout errors
    # If a timeout occurs, the last line of valid data is repeated in the next file and is skipped, but if that is the only line of data, script breaks
    # Remove "status" error rows
    if (
        len(
            df[
                df.isin([' "message": "Request could not complete. Timeout."}']).any(
                    axis=1
                )
            ]
        )
        != 0
    ):
        df = df.drop(df.tail(1).index)
    elif len(df[df.isin(['{"status": 408']).any(axis=1)]) != 0:
        df = df.drop(df.tail(1).index)

    # Drop any columns that only contain NAs.
    df = df.dropna(axis=1, how="all")

    # If columns duplicated, check to see if they are identical and drop second if so.
    if headers["dup"] is True:
        for i in headers["dup_col"]:
            # For each duplicated column
            cols = df.filter(like=i).columns
            if df[cols[0]].equals(df[cols[1]]):
                df.drop(cols[1], axis=1, inplace=True)

            else:
                print("Non-identical duplicate columns found.")
                errors["File"].append(file)
                errors["Time"].append(end_api)
                errors["Error"].append(
                    f"Non-identical duplicate columns found. Columns: {i}"
                )
                continue

    # Fix time format issues caused by "Timeout" errors.
    # Convert time to datetime and catch any incorrectly formatted columns.
    df["Date_Time"] = pd.to_datetime(df["Date_Time"], errors="coerce")
    # Remove any rows where time is missing.
    df = df[pd.notnull(df["Date_Time"])]

    # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022.
    df = df.loc[(df["Date_Time"] < "09-01-2022") & (df["Date_Time"] > "12-31-1979")]

    # note this will return an empty df if all data is outside these time bounds
    if len(df.index) == 0:
        print(
            f"No data for {file} during v1 period (1/1980 - 8/2022), station not cleaned."
        )
        return None

    # Remove any non-essential columns.
    coltokeep = [
        "Station_ID",
        "Date_Time",
        "altimeter_set_1",
        "altimeter_set_1_qc",
        "air_temp_set_1",
        "air_temp_set_1_qc",
        "relative_humidity_set_1",
        "relative_humidity_set_1_qc",
        "wind_speed_set_1",
        "wind_speed_set_1_qc",
        "wind_direction_set_1",
        "wind_direction_set_1_qc",
        "precip_accum_since_local_midnight_set_1",
        "precip_accum_since_local_midnight_set_1_qc",
        "precip_accum_24_hour_set_1",
        "precip_accum_24_hour_set_1_qc",
        "dew_point_temperature_set_1d",
        "pressure_set_1d",
        "pressure_set_1",
        "pressure_set_1_qc",
        "dew_point_temperature_set_1",
        "dew_point_temperature_set_1_qc",
        "precip_accum_one_hour_set_1",
        "precip_accum_one_hour_set_1_qc",
        "solar_radiation_set_1",
        "solar_radiation_set_1_qc",
        "precip_accum_set_1",
        "precip_accum_set_1_qc",
        "precip_accum_five_minute_set_1",
        "precip_accum_five_minute_set_1_qc",
    ]

    othercols = [
        col for col in df.columns if col not in coltokeep and col not in removedvars
    ]
    # Add any new columns from drop list to removedvars, to save later.
    removedvars += othercols
    # Drop all columns not in coltokeep list.
    df = df.drop(columns=[col for col in df if col not in coltokeep])

    # Manually convert "None" to np.nan
    df.replace(to_replace="None", value=np.nan, inplace=True)

    return df


def clean_madis(
    rawdir: str, cleandir: str, network: str, cwop_letter: str | None = None
):
    """Cleans MADIS data for a variety of networks.

    Paramters
    ---------
    rawdir : str
        path to raw data bucket
    cleandir : str
        path to cleaned data bucket
    network : str
        name of network
    cwop_letter : str, optional
        option to subset CWOP stations by start letter, recommended!!

    Returns
    -------
    None
    """

    # Ensuring that non-CWOP networks do not accidentally subset
    if network != "CWOP":
        cwop_letter = None

    try:
        files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=rawdir):
            file = str(item.key)
            files += [file]

        # Get list of file names
        files = list(filter(lambda f: f.endswith(".csv"), files))
        # Remove error handling files.
        files = [file for file in files if "error" not in file]
        # Remove error handling files.
        files = [file for file in files if "station" not in file]

        # Set up error handling.
        errors = {"File": [], "Time": [], "Error": []}
        # Set end time to be current time at beginning of download: for error handling csv.
        end_api = datetime.now().strftime("%Y%m%d%H%M")
        # For attributes of netCDF file.
        timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

        # Set up list of variables to be removed
        removedvars = []

        # # Get list of station IDs from filename and clean.
        ids = list()
        for file in files:
            id = file.split("/")[-1]  # Remove leading folders
            id = re.sub(".csv", "", id)  # Remove file extension
            id_mult = id.split("_")  # Split on _ char in two cases
            if len(id_mult[0]) <= 2:  # timeout split files go up to 11_STID so far
                id = id_mult[-1]  # X_STID timeout split files, id is after _
            elif len(id_mult[-1]) == 10:
                id = id_mult[0]  # STID_20YY-MM-DD update pull files, id is before _
            if id not in ids:
                ids.append(id)

        # Get sensor metadata from QA/QC folder
        sensor_filepath = (
            f"s3://wecc-historical-wx/3_qaqc_wx/{network}/sensorlist_{network}.csv"
        )
        sensor_data = pd.read_csv(smart_open.smart_open(sensor_filepath))

    except Exception as e:
        # If unable to read files from cleandir, break function.
        print("Whole network error.")
        errors["File"].append("Whole network")
        errors["Time"].append(end_api)
        errors["Error"].append(f"Whole network error: {e}")

    else:
        # If files read successfully, continue.
        dfs = []  # intialize empty df for appending

        # Procedure for grouping of data in CWOP to split up 7k+ stations by first letter

        # catch-all single letter stations (K, L, M, P, S, T, U, W at present)
        not_ABCDEFG = (
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
        )
        if network == "CWOP" and cwop_letter != None:
            # cwop_letter = "other"
            if "other" in cwop_letter and len(cwop_letter) == 5:
                ids = [id for id in ids if not id.startswith(not_ABCDEFG)]

            elif "other" in cwop_letter and len(cwop_letter) != 5:
                # additional letters + other category called, ex: cwop_letter = "ABC + other"
                letter_to_clean = cwop_letter.replace(" ", "")
                # so it doesn't clean "o t h e r"
                letter_to_clean = letter_to_clean.replace("other", "")
                letter_to_clean = letter_to_clean.replace("+", "")
                letter_ids = tuple(letter_to_clean)
                other_ids = [id for id in ids if not id.startswith(not_ABCDEFG)]
                letter_ids = [id for id in ids if id.startswith(letter_ids)]
                ids = other_ids + letter_ids

            if len(cwop_letter) == 1:
                # single letter cleaning, ex: cwop_letter = "A"
                ids = [id for id in ids if id.startswith(str(cwop_letter))]

            if "other" not in cwop_letter and len(cwop_letter) != 1:
                # more than one letter provided, but not other category, ex: cwop_letter = "ACD"
                letter_ids = tuple(cwop_letter)
                ids = [id for id in ids if id.startswith(letter_ids)]

            print(
                f"CWOP batch cleaning for '{cwop_letter}' stations: batch-size of {len(ids)} stations"
            )

        elif network == "CWOP" and cwop_letter == None:
            # This a full network clean with no batch sub-setting, ex: cwop_letter = None
            print(
                "Warning: Setting cwop_letter = None is for an entire network clean of CWOP, estimated 1 week of continuous runtime to complete."
            )
            ids = ids

        else:
            # network should not be CWOP, and will complete full clean
            ids = ids

        for i in ids:
            try:
                # Get list of files with station ID in them.
                stat_files = [k for k in files if i in k]
                # Save file ID as uppercase always.
                station_id = f"{network}_{i.upper()}"
                # Iterate through files to clean. Each file represents a station's data.
                headers = []

                if not stat_files:
                    # If no files left in list
                    print(f"No raw data found for {station_id} on AWS.")
                    errors["File"].append(station_id)
                    errors["Time"].append(end_api)
                    errors["Error"].append("No raw data found for this station on AWS.")
                    continue  # Skip this station

                for file in stat_files:
                    try:
                        skip = 0
                        header = parse_madis_headers(file)
                        # If units are NaN, this signifies an empty dataframe. Write to errors, but do not clean station.
                        if isinstance(header["units"], float) and np.isnan(
                            header["units"]
                        ):
                            print(
                                f"{station_id} reports empty data file -- not cleaned."
                            )
                            errors["File"].append(file)
                            errors["Time"].append(end_api)
                            errors["Error"].append(
                                "No data available for station. Cleaning stage skipped."
                            )
                            stat_files.remove(file)  # Remove from file list
                            continue
                        headers.append(header)

                    except Exception as e:
                        print(
                            f"Error parsing MADIS headers, please check for {station_id}."
                        )
                        errors["File"].append(file)
                        errors["Time"].append(end_api)
                        errors["Error"].append("Error parsing MADIS headers.")
                        continue

                # If more than one file for station, metadata should be identical. Test this.
                if all(a == headers[0] for a in headers[1:]):
                    headers = headers[0]
                else:
                    # If not, provide user input option to proceed.
                    ### TEMPORARY MISMATCH METADATA FIX
                    headers = headers[0]
                    # resp = input("Station files provide conflicting metadata for station {0}. Which file's metadata would you like to use (0-{1})?".format(station_id, len(headers)-1))
                    # try:
                    #     headers = headers[int(resp)]
                    # except ValueError:
                    #     errors['File'].append(stat_files)
                    #     errors['Time'].append(end_api)
                    #     errors['Error'].append("Invalid response to metadata query. Skipping file.")
                    #     print("Invalid response. Skipping file {}.")
                    #     continue
                    ###

                try:
                    dfs = [
                        parse_madis_to_pandas(file, headers, errors, removedvars)
                        for file in stat_files
                    ]

                except Exception as e:
                    # Note: error handling will be slightly coarse here bc of list comprehension.
                    errors["File"].append(stat_files)
                    errors["Time"].append(end_api)
                    errors["Error"].append(
                        f"Error in parsing MADIS files to pandas dfs: {e}"
                    )
                    continue

            except Exception as e:
                # If stat_files is none, this will default to saving ID of station.
                errors["File"].append(i)
                errors["Time"].append(end_api)
                errors["Error"].append(f"Error in stat_files set-up: {e}")
                continue

            try:
                station_name = headers["station_name"]  # Get station name
                file_count = len(dfs)
                if file_count == 0:
                    # Example: # AP907 has a known timeout issue
                    print(f"{i} is having dfs appending issues, please check")
                    errors["File"].append(i)
                    errors["Time"].append(end_api)
                    errors["Error"].append("Dataframe appending issue, please check")
                    continue  # skip station

                elif all(df is None for df in dfs) == True:
                    # If all files for a station do not have data within time bound range
                    errors["File"].append(stat_files)
                    errors["Time"].append(end_api)
                    errors["Error"].append(
                        "All raw data files do not report data for 1/1980 - 8/2022. Not cleaned"
                    )
                    continue  # skip station

                # fixes mismatch in num. of cols if sensor is added to station in newer raw datafile
                df_stat = pd.concat(dfs, axis=0, ignore_index=True)

                # Deal with units
                units = pd.DataFrame(
                    list(zip(headers["columns"], list(headers["units"].split(",")))),
                    columns=["column", "units"],
                )
                varstokeep = list(df_stat.columns)
                # Only keep non-removed cols
                units = units[units.column.isin(varstokeep)]
                # Only keep non qa-qc data types
                units = units[units["units"] != "QC_type"]
                # Pull code set to get units as metric. Break if this condition is not met.
                if "Fahrenheit" in units["units"]:
                    print("Units not standardized! Fix code")
                    exit()

                # Fix multi-type columns
                # If column has QC in it, force to string.
                for b in df_stat.columns:
                    multitype = set(type(x).__name__ for x in df_stat[b])
                    if len(multitype) > 1:
                        if "qc" in b:
                            # Coerce to string (to handle multiple QA/QC flags)
                            df_stat[b] = df_stat[b].astype(str)
                            # Remove trailing .0 from float conversion.
                            df_stat[b] = df_stat[b].str.replace(".0", "")
                        elif "wind_cardinal_direction" in b:
                            df_stat[b] = df_stat[b].astype(str)  # Coerce to string
                        elif "sea_level_pressure_set" in b:
                            df_stat[b] = df_stat[b].astype(float)  # Coerce to float.
                        elif "wind_gust_set_1" in b:
                            df_stat[b] = df_stat[b].astype(float)  # Coerce to float.
                        elif "heat_index_set_1" in b:
                            df_stat[b] = df_stat[b].astype(float)  # Coerce to float.
                        elif "wind_direction_set_1" in b:
                            df_stat[b] = df_stat[b].astype(float)  # Coerce to float.
                        else:
                            # Code to flag novel exceptions, correct and add explicit handling above.
                            print(
                                f"Multitype error for column {b} with data types {multitype}. Please resolve"
                            )
                            errors["File"].append(file)
                            errors["Time"].append(end_api)
                            errors["Error"].append(
                                f"Multitype error for column {b} with data types {multitype}. Please resolve"
                            )
                            continue
                    else:
                        if "qc" in b:
                            # Coerce QA/QC flag to string in all instances.
                            df_stat[b] = df_stat[b].astype(str)
                            # Remove trailing .0 from float conversion.
                            df_stat[b] = df_stat[b].str.replace(".0", "")

                # Fix issue with "nan" and nan causing comparison errors
                df_stat = df_stat.replace("nan", np.nan)

                # Sort by time and remove any overlapping timestamps.
                df_stat = df_stat.sort_values(by="Date_Time")
                df_stat = df_stat.drop_duplicates()

                # Move df to xarray object.
                ds = df_stat.to_xarray()

                # Update global attributes
                ds = ds.assign_attrs(title=f"{network} cleaned")
                ds = ds.assign_attrs(institution="Eagle Rock Analytics / Cal Adapt")
                ds = ds.assign_attrs(source="")
                ds = ds.assign_attrs(
                    history=f"MADIS_clean.py script run on {timestamp} UTC"
                )
                ds = ds.assign_attrs(
                    comment="Intermediate data product: may not have been subject to any cleaning or QA/QC processing"
                )
                ds = ds.assign_attrs(license="")
                ds = ds.assign_attrs(citation="")
                ds = ds.assign_attrs(
                    disclaimer="This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document."
                )
                ds = ds.assign_attrs(station_name=station_name.replace("\n", ""))
                # Keep count of how many files merged per station.
                ds = ds.assign_attrs(raw_files_merged=file_count)

                # Update dimensions and coordinates

                # Add dimensions: station ID and time.
                ds = ds.rename({"Date_Time": "time"})  # Rename time variable.
                # Swap index with time.
                ds = ds.set_coords("time").swap_dims({"index": "time"})
                ds = ds.assign_coords(id=str(station_id))
                ds = ds.expand_dims("id")  # Add station_id as index.
                # Drop station_id variable and index coordinate.
                ds = ds.drop_vars(("index"))
                # Rename id to station_id.
                ds = ds.rename({"id": "station"})

                # Add coordinates: latitude and longitude.
                lat = np.asarray([headers["latitude"]] * len(ds["time"]))
                lat.shape = (1, len(ds["time"]))

                lon = np.asarray([headers["longitude"]] * len(ds["time"]))
                lon.shape = (1, len(ds["time"]))

                # reassign lat and lon as coordinates
                ds = ds.assign_coords(
                    lat=(["station", "time"], lat), lon=(["station", "time"], lon)
                )

                # If any observation is missing lat or lon coordinates, drop these observations.
                if np.count_nonzero(np.isnan(ds["lat"])) != 0:
                    ds = ds.where(~np.isnan(ds["lat"]))

                if np.count_nonzero(np.isnan(ds["lon"])) != 0:
                    ds = ds.where(~np.isnan(ds["lon"]))

                # Add variable: elevation
                elev = np.asarray([headers["elevation"]] * len(ds["time"]))
                elev.shape = (1, len(ds["time"]))
                ds["elevation"] = (["station", "time"], elev)

                # Update dimension and coordinate attributes.

                # Convert column to datetime (and remove any rows that cannot be coerced).
                ds["time"] = pd.to_datetime(ds["time"], utc=True)
                ds["time"] = pd.to_datetime(ds["time"], unit="ns")  # Fix time format.

                # Update attributes.
                ds["time"].attrs["long_name"] = "time"
                ds["time"].attrs["standard_name"] = "time"
                ds["time"].attrs["comment"] = "In UTC."

                # Station ID
                ds["station"].attrs["long_name"] = "station_id"
                ds["station"].attrs[
                    "comment"
                ] = "Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network."

                # Latitude
                ds["lat"].attrs["long_name"] = "latitude"
                ds["lat"].attrs["standard_name"] = "latitude"
                ds["lat"].attrs["units"] = "degrees_north"

                # Longitude
                ds["lon"].attrs["long_name"] = "longitude"
                ds["lon"].attrs["standard_name"] = "longitude"
                ds["lon"].attrs["units"] = "degrees_east"

                # Elevation
                # Convert from feet to meters.
                ds["elevation"] = calc_clean._unit_elev_ft_to_m(ds["elevation"])

                ds["elevation"].attrs["standard_name"] = "height_above_mean_sea_level"
                ds["elevation"].attrs["long_name"] = "station_elevation"
                ds["elevation"].attrs["units"] = "meters"
                # Define which direction is positive
                ds["elevation"].attrs["positive"] = "up"
                ds["elevation"].attrs["comment"] = "Converted from feet to meters."

                # Update sensor metadata
                # May be multiple rows if sensors added/removed over time.
                station_sensors = sensor_data.loc[sensor_data.STID == i]
                # Get all position columns, dropping duplicate rows
                sensorheights = station_sensors[
                    [x for x in station_sensors.columns if "position" in x]
                ].drop_duplicates()

                # If sensor heights is completely NA, ignore.
                if sensorheights.isnull().all().all():
                    ds = ds.assign_attrs(anemometer_height_m=np.nan)
                    ds = ds.assign_attrs(thermometer_height_m=np.nan)
                    ds = ds.assign_attrs(barometer_elevation_m=np.nan)

                else:
                    # For wind direction, air temp and barometer sensors, keep multiple sensor heights if there are more than one.

                    # Wind speed & direction (m)
                    if (
                        "wind_speed_1_position" in sensorheights.columns
                        and True in sensorheights.wind_speed_1_position.notnull().values
                    ):
                        # If any value not null
                        if (
                            len(sensorheights.wind_speed_1_position) > 1
                            and sensorheights.wind_speed_1_position.nunique() != 1
                        ):
                            print(
                                f"{station_id} has more than one wind sensor height. Check these are saved as expected."
                            )
                            # Testing, to flag cases of this.
                            # If more than one anemometer sensor w/ diff heights
                            # Use time column to name attributes

                            # Get start date for each sensor height
                            start_dates = (
                                station_sensors.sort_values("wind_speed_1_start")
                                .groupby("wind_speed_1_position")
                                .head(1)
                            )
                            start_dates = start_dates[
                                ["wind_speed_1_position", "wind_speed_1_start"]
                            ]

                            # Get end date for each sensor height
                            end_dates = (
                                station_sensors.sort_values("wind_speed_1_end")
                                .groupby("wind_speed_1_position")
                                .tail(1)
                            )
                            end_dates = end_dates[
                                ["wind_speed_1_position", "wind_speed_1_end"]
                            ]

                            # Join by sensor heights
                            dates = start_dates.merge(
                                end_dates, on="wind_speed_1_position"
                            )
                            row["names"] = np.nan  # Add names column

                            # generate attribute names
                            for index, row in dates.iterrows():
                                row["names"] = (
                                    f"anemometer_height_m_{row.wind_speed_1_start[0:10]}_{row.wind_speed_1_end[0:10]}"
                                )
                                ds.attrs[row["names"]] = float(
                                    row["wind_speed_1_position"]
                                )

                        else:
                            sensor_height = (
                                sensorheights.wind_speed_1_position.dropna().unique()
                            )
                            # Get all unique values from column.
                            ds = ds.assign_attrs(
                                anemometer_height_m=float(sensor_height[0])
                            )
                    else:
                        ds = ds.assign_attrs(anemometer_height_m=np.nan)

                    # Air temperature (m)
                    if (
                        "air_temp_1_position" in sensorheights.columns
                        and True in sensorheights.air_temp_1_position.notnull().values
                    ):
                        # If any value not null
                        if (
                            len(sensorheights.air_temp_1_position) > 1
                            and sensorheights.air_temp_1_position.nunique() != 1
                        ):
                            print(
                                f"{station_id} has more than one air temp sensor height. Check these are saved as expected."
                            )
                            # Testing, to flag cases of this.
                            # If more than one anemometer sensor w/ diff heights
                            # Use time column to name attributes

                            # Get start date for each sensor height
                            start_dates = (
                                station_sensors.sort_values("air_temp_1_start")
                                .groupby("air_temp_1_position")
                                .head(1)
                            )
                            start_dates = start_dates[
                                ["air_temp_1_position", "air_temp_1_start"]
                            ]

                            # Get end date for each sensor height
                            end_dates = (
                                station_sensors.sort_values("air_temp_1_end")
                                .groupby("air_temp_1_position")
                                .tail(1)
                            )
                            end_dates = end_dates[
                                ["air_temp_1_position", "air_temp_1_end"]
                            ]

                            # Join by sensor heights
                            dates = start_dates.merge(
                                end_dates, on="air_temp_1_position"
                            )
                            row["names"] = np.nan  # Add names column

                            # generate attribute names
                            for index, row in dates.iterrows():
                                row["names"] = (
                                    f"thermometer_height_m_{row.air_temp_1_start[0:10]}_{row.air_temp_1_end[0:10]}"
                                )
                                ds.attrs[row["names"]] = float(
                                    row["air_temp_1_position"]
                                )
                        else:
                            sensor_height = (
                                sensorheights.air_temp_1_position.dropna().unique()
                            )
                            # Get all unique values from column.
                            ds = ds.assign_attrs(
                                thermometer_height_m=float(sensor_height[0])
                            )
                    else:
                        ds = ds.assign_attrs(thermometer_height_m=np.nan)

                    # Barometer elevation (convert from height)
                    if (
                        "pressure_1_position" in sensorheights.columns
                        and True in sensorheights.pressure_1_position.notnull().values
                    ):
                        # If any value not null
                        if (
                            len(sensorheights.pressure_1_position) > 1
                            and sensorheights.pressure_1_position.nunique() != 1
                        ):
                            print(
                                f"{station_id} has more than one pressure sensor height. Check these are saved as expected."
                            )
                            # Testing, to flag cases of this.
                            # If more than one anemometer sensor w/ diff heights
                            # Use time column to name attributes

                            # Get start date for each sensor height
                            start_dates = (
                                station_sensors.sort_values("pressure_1_start")
                                .groupby("pressure_1_position")
                                .head(1)
                            )
                            start_dates = start_dates[
                                ["pressure_1_position", "pressure_1_start"]
                            ]

                            # Get end date for each sensor height
                            end_dates = (
                                station_sensors.sort_values("pressure_1_end")
                                .groupby("pressure_1_position")
                                .tail(1)
                            )
                            end_dates = end_dates[
                                ["pressure_1_position", "pressure_1_end"]
                            ]

                            # Join by sensor heights
                            dates = start_dates.merge(
                                end_dates, on="pressure_1_position"
                            )
                            row["names"] = np.nan  # Add names column

                            if pd.notnull(ds["elevation"].values[0]):
                                # generate attribute names
                                for index, row in dates.iterrows():
                                    row["names"] = (
                                        f"barometer_elevation_m_{row.pressure_1_start[0:10]}_{row.pressure_1_end[0:10]}"
                                    )
                                    ds.attrs[row["names"]] = float(
                                        row["pressure_1_position"]
                                    ) + float(ds["elevation"].values[0])
                            else:
                                for index, row in dates.iterrows():
                                    row["names"] = (
                                        f"barometer_height_m_{row.pressure_1_start[0:10]}_{row.pressure_1_end[0:10]}"
                                    )
                                    ds.attrs[row["names"]] = float(
                                        row["pressure_1_position"]
                                    )
                                    ds = ds.assign_attrs(barometer_elevation_m=np.nan)

                        else:
                            if pd.notnull(ds["elevation"].values[0]):
                                # If station has elevation
                                sensor_height = (
                                    sensorheights.pressure_1_position.dropna().unique()
                                )
                                # Get all unique values from column.
                                barometer_elev = sensor_height[0] + float(
                                    ds["elevation"].values[0]
                                )
                                ds = ds.assign_attrs(
                                    barometer_elevation_m=barometer_elev
                                )
                            else:
                                # if no station elevation, keep barometer height and record barometer elevation as NaN
                                sensor_height = (
                                    sensorheights.pressure_1_position.dropna().unique()
                                )
                                # Get all unique values from column.
                                ds = ds.assign_attrs(
                                    barometer_height_m=sensor_height[0]
                                )
                                ds = ds.assign_attrs(barometer_elevation_m=np.nan)

                    else:
                        ds = ds.assign_attrs(barometer_elevation_m=np.nan)

                    # Any other sensors with values
                    sensorheights = sensorheights.dropna(axis="columns", how="all")

                    # For precip columns, get single value
                    precip_heights = sensorheights[
                        [col for col in sensorheights.columns if "precip" in col]
                    ]
                    if not precip_heights.isnull().values.all():
                        if len(precip_heights.columns) > 1:
                            # If multiple rainfall columns
                            if precip_heights.nunique(axis=1).eq(1).all():
                                # If all values the same in row
                                ds = ds.assign_attrs(
                                    rain_gauge_height_m=float(
                                        precip_heights.values[0][0]
                                    )
                                )
                            else:
                                print("Competing rain gauge height values")
                                exit()
                        else:
                            ds = ds.assign_attrs(
                                rain_gauge_height_m=float(precip_heights.values[0])
                            )

                    # Currently, these only take one sensor height if there is more than one provided.
                    for col in sensorheights.columns:
                        if col not in [
                            "altimeter_1_position",
                            "wind_speed_1_position",
                            "pressure_1_position",
                        ]:
                            sensor_height = sensorheights[col].dropna().unique()
                            if col == "relative_humidity_1_position":
                                ds = ds.assign_attrs(
                                    humidity_height_m=float(sensor_height[0])
                                )
                            if col == "dew_point_temperature_1_position":
                                ds = ds.assign_attrs(
                                    dew_point_temperature_height_m=float(
                                        sensor_height[0]
                                    )
                                )
                            if col == "solar_radiation_1_position":
                                ds = ds.assign_attrs(
                                    pyranometer_height_m=float(sensor_height[0])
                                )
                            if col == "wind_direction_1_position":
                                ds = ds.assign_attrs(
                                    wind_vane_height_m=float(sensor_height[0])
                                )

                # Update variable attributes and do unit conversions

                # tas: air surface temperature (K)
                if "air_temp_set_1" in ds.keys():
                    ds["tas"] = calc_clean._unit_degC_to_K(ds["air_temp_set_1"])
                    ds = ds.drop("air_temp_set_1")

                    if "air_temp_set_1_qc" in ds.keys():
                        # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
                        flagvals = ds["air_temp_set_1_qc"].values.tolist()[0]
                        flagvals = [x for x in flagvals if pd.isnull(x) == False]
                        flagvals = list(np.unique(flagvals))  # Get unique values
                        # Split any rows with multiple flags and run unique again.
                        flagvals = list(
                            np.unique(
                                [
                                    split_item
                                    for item in flagvals
                                    for split_item in str(item).split(";")
                                ]
                            )
                        )

                        if flagvals == ["nan"]:
                            # This should not occur, but leave in here for additional robustness.
                            ds = ds.drop("air_temp_set_1_qc")
                            ds["tas"].attrs["long_name"] = "air_temperature"
                            ds["tas"].attrs["standard_name"] = "air_temperature"
                            ds["tas"].attrs["units"] = "degree_Kelvin"
                            ds["tas"].attrs[
                                "comment"
                            ] = "Converted from Celsius to Kelvin."

                        else:
                            ds = ds.rename({"air_temp_set_1_qc": "tas_qc"})
                            # Reformat to match csv.
                            flagvals = [ele.replace(".0", "") for ele in flagvals]
                            ds["tas_qc"].attrs["flag_values"] = flagvals
                            # Generate values from unique values from dataset.
                            ds["tas_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."

                            ds["tas"].attrs["long_name"] = "air_temperature"
                            ds["tas"].attrs["standard_name"] = "air_temperature"
                            ds["tas"].attrs["units"] = "degree_Kelvin"
                            ds["tas"].attrs["ancillary_variables"] = "tas_qc"
                            # List other variables associated with variable (QA/QC)
                            ds["tas"].attrs[
                                "comment"
                            ] = "Converted from Celsius to Kelvin."

                    else:
                        ds["tas"].attrs["long_name"] = "air_temperature"
                        ds["tas"].attrs["standard_name"] = "air_temperature"
                        ds["tas"].attrs["units"] = "degree_Kelvin"
                        ds["tas"].attrs["comment"] = "Converted from Celsius to Kelvin."

                # ps: surface air pressure (Pa)
                # Note here that if "pressure_set_1" has values this is a direct station observation reading.
                # Otherwise, if "pressure_set_1d" has values this is a derived value calculated from altimeter and elevation.
                # We will manually recalculate this here.
                if "pressure_set_1" in ds.keys():
                    # If station pressure directly observed
                    if not np.isnan(ds["pressure_set_1"].values).all():
                        # If station pressure directly observed
                        ds = ds.rename({"pressure_set_1": "ps"})

                        # Set attributes
                        ds["ps"].attrs["long_name"] = "station_air_pressure"
                        ds["ps"].attrs["standard_name"] = "air_pressure"
                        ds["ps"].attrs["units"] = "Pa"

                        if "sea_level_pressure_set_1" in ds.keys():
                            ds = ds.drop("sea_level_pressure_set_1")
                            # Drop psl if station pressure available

                if "ps" not in ds.keys():
                    # If this didn't work, look for sea level pressure
                    if "sea_level_pressure_set_1" in ds.keys():
                        ds = ds.rename({"sea_level_pressure_set_1": "psl"})
                        ds["psl"].attrs["long_name"] = "sea_level_air_pressure"
                        ds["psl"].attrs["standard_name"] = "air_pressure"
                        ds["psl"].attrs["units"] = "Pa"

                if "pressure_set_1_qc" in ds.keys():
                    # If QA/QC exists.
                    flagvals = ds["pressure_set_1_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in here for additional check.
                        # print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pressure_set_1_qc")
                    else:
                        ds = ds.rename({"pressure_set_1_qc": "ps_qc"})
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        ds["ps_qc"].attrs["flag_values"] = flagvals
                        # Generate values from unique values from dataset.
                        ds["ps_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        if "ps" in ds.keys():
                            ds["ps"].attrs["ancillary_variables"] = "ps_qc"

                if "sea_level_pressure_set_1_qc" in ds.keys():
                    # If QA/QC exists.
                    if "psl" in ds.keys():
                        flagvals = ds["sea_level_pressure_set_1_qc"].values.tolist()[0]
                        flagvals = [x for x in flagvals if pd.isnull(x) == False]
                        flagvals = list(np.unique(flagvals))  # Get unique values
                        flagvals = list(
                            np.unique(
                                [
                                    split_item
                                    for item in flagvals
                                    for split_item in str(item).split(";")
                                ]
                            )
                        )
                        # Split any rows with multiple flags and run unique again.

                        if flagvals == ["nan"]:
                            # This should not occur, but leave in as additional check.
                            ds = ds.drop("sea_level_pressure_set_1_qc")
                        else:
                            ds = ds.rename({"sea_level_pressure_set_1_qc": "psl_qc"})
                            # Reformat to match csv.
                            flagvals = [ele.replace(".0", "") for ele in flagvals]
                            # Generate values from unique values from dataset.
                            ds["psl_qc"].attrs["flag_values"] = flagvals
                            ds["psl_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."
                            ds["psl"].attrs["ancillary_variables"] = "psl_qc"
                    else:
                        # if no sea level pressure, drop QC flags.
                        ds = ds.drop("sea_level_pressure_set_1_qc")

                # tdps: dew point temperature (K)
                # if raw dew point temperature observed, use that.
                if "dew_point_temperature_set_1" in ds.keys():
                    ds["tdps"] = calc_clean._unit_degC_to_K(
                        ds["dew_point_temperature_set_1"]
                    )
                    ds = ds.drop("dew_point_temperature_set_1")

                    # Set attributes for conversion.
                    ds["tdps"].attrs["long_name"] = "dew_point_temperature"
                    ds["tdps"].attrs["standard_name"] = "dew_point_temperature"
                    ds["tdps"].attrs["units"] = "degree_Kelvin"
                    ds["tdps"].attrs["comment"] = "Converted from Celsius to Kelvin."

                    # QAQC flag
                    if "dew_point_temperature_set_1_qc" in ds.keys():
                        # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
                        flagvals = ds["dew_point_temperature_set_1_qc"].values.tolist()[
                            0
                        ]
                        flagvals = [x for x in flagvals if pd.isnull(x) == False]
                        flagvals = list(np.unique(flagvals))  # Get unique values
                        flagvals = list(
                            np.unique(
                                [
                                    split_item
                                    for item in flagvals
                                    for split_item in str(item).split(";")
                                ]
                            )
                        )
                        # Split any rows with multiple flags and run unique again.

                        if flagvals == ["nan"]:
                            # This should not occur, but leave in as additional check.
                            # print("Flag value is {}".format(flagvals)) # For testing.
                            ds = ds.drop("dew_point_temperature_set_1_qc")
                        else:
                            ds = ds.rename(
                                {"dew_point_temperature_set_1_qc": "tdps_qc"}
                            )
                            # Reformat to match csv.
                            flagvals = [ele.replace(".0", "") for ele in flagvals]
                            # Generate values from unique values from dataset.
                            ds["tdps_qc"].attrs["flag_values"] = flagvals
                            ds["tdps_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."

                            ds["tdps"].attrs["ancillary_variables"] = "tdps_qc"
                            # List other variables associated with variable (QA/QC)

                # pr: precipitation
                # We have 4 different raw precipitation variables for precip.
                # precip_accum_24_hour_set_1 # Precipitation from last 24 hours.
                # precip_accum_since_local_midnight_set_1 # Precipitation from local midnight.
                # precip_accum_set_1 # Precipitation since last record.
                # precip_accum_one_hour_set_1 # Precipitation in last hour.

                # At this stage, no infilling. So we will keep all columns with data and simply rename them.

                # Reformat remaining columns
                if "precip_accum_24_hour_set_1" in ds.keys():
                    ds = ds.rename({"precip_accum_24_hour_set_1": "pr_24h"})
                    ds["pr_24h"].attrs["long_name"] = "24_hr_precipitation_amount"
                    ds["pr_24h"].attrs["units"] = "mm/24hr"
                    ds["pr_24h"].attrs[
                        "comment"
                    ] = "Precipitation accumulated in previous 24 hour period."

                    if "precip_accum_24_hour_set_1_qc" in ds.keys():
                        ds = ds.rename({"precip_accum_24_hour_set_1_qc": "pr_24h_qc"})

                if "precip_accum_since_local_midnight_set_1" in ds.keys():
                    ds = ds.rename(
                        {"precip_accum_since_local_midnight_set_1": "pr_localmid"}
                    )
                    ds["pr_localmid"].attrs[
                        "long_name"
                    ] = "precipitation_since_local_midnight"
                    ds["pr_localmid"].attrs["units"] = "mm"
                    # since local midnight, including "midnight" in unit will break xr read in for time unit
                    ds["pr_localmid"].attrs[
                        "comment"
                    ] = "Precipitation accumulated since local midnight."

                    if "precip_accum_since_local_midnight_set_1_qc" in ds.keys():
                        ds = ds.rename(
                            {
                                "precip_accum_since_local_midnight_set_1_qc": "pr_localmid_qc"
                            }
                        )

                if "precip_accum_set_1" in ds.keys():
                    ds = ds.rename({"precip_accum_set_1": "pr"})
                    ds["pr"].attrs["long_name"] = "precipitation_amount"
                    ds["pr"].attrs["units"] = "mm/interval"
                    ds["pr"].attrs[
                        "comment"
                    ] = "Precipitation accumulated since previous measurement."

                    if "precip_accum_set_1_qc" in ds.keys():
                        ds = ds.rename({"precip_accum_set_1_qc": "pr_qc"})

                if "precip_accum_one_hour_set_1" in ds.keys():
                    ds = ds.rename({"precip_accum_one_hour_set_1": "pr_1h"})
                    ds["pr_1h"].attrs["long_name"] = "hourly_precipitation_amount"
                    ds["pr_1h"].attrs["units"] = "mm/hr"
                    ds["pr_1h"].attrs[
                        "comment"
                    ] = "Precipitation accumulated in previous hour."

                    if "precip_accum_one_hour_set_1_qc" in ds.keys():
                        ds = ds.rename({"precip_accum_one_hour_set_1_qc": "pr_1h_qc"})

                if "precip_accum_five_minute_set_1" in ds.keys():
                    ds = ds.rename({"precip_accum_five_minute_set_1": "pr_5min"})
                    ds["pr_5min"].attrs["long_name"] = "5_minute_precipitation_amount"
                    ds["pr_5min"].attrs["units"] = "mm/5 min"
                    ds["pr_5min"].attrs[
                        "comment"
                    ] = "Precipitation accumulated in previous 5 minutes."

                    if "precip_accum_five_minute_set_1_qc" in ds.keys():
                        ds = ds.rename(
                            {"precip_accum_five_minute_set_1_qc": "pr_5min_qc"}
                        )

                # Reformat qc columns
                if "pr_24h_qc" in ds.keys():
                    flagvals = ds["pr_24h_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        # print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pr_24h_qc")
                    else:
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        # Generate values from unique values from dataset.
                        ds["pr_24h_qc"].attrs["flag_values"] = flagvals
                        ds["pr_24h_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["pr_24h"].attrs["ancillary_variables"] = "pr_24h_qc"

                if "pr_localmid_qc" in ds.keys():
                    flagvals = ds["pr_localmid_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("pr_localmid_qc")
                    else:
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        # Generate values from unique values from dataset.
                        ds["pr_localmid_qc"].attrs["flag_values"] = flagvals
                        ds["pr_localmid_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["pr_localmid"].attrs[
                            "ancillary_variables"
                        ] = "pr_localmid_qc"

                if "pr_qc" in ds.keys():
                    flagvals = ds["pr_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("pr_qc")
                    else:
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        # Generate values from unique values from dataset.
                        ds["pr_qc"].attrs["flag_values"] = flagvals
                        ds["pr_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["pr"].attrs["ancillary_variables"] = "pr_qc"

                if "pr_1h_qc" in ds.keys():
                    flagvals = ds["pr_1h_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("pr_1h_qc")
                    else:
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        # Generate values from unique values from dataset.
                        ds["pr_1h_qc"].attrs["flag_values"] = flagvals
                        ds["pr_1h_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["pr_1h"].attrs["ancillary_variables"] = "pr_1h_qc"

                if "pr_5min_qc" in ds.keys():
                    flagvals = ds["pr_5min_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("pr_5min_qc")
                    else:
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        # Generate values from unique values from dataset.
                        ds["pr_5min_qc"].attrs["flag_values"] = flagvals
                        ds["pr_5min_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["pr_5min"].attrs["ancillary_variables"] = "pr_5min_qc"

                # Set ancillary variables based on other precip cols
                precip_cols = [elem for elem in ds.keys() if "pr" in elem]
                precip_cols = [elem for elem in precip_cols if "pressure" not in elem]

                if precip_cols:
                    # if list not empty
                    if len(precip_cols) > 1:
                        # If list has more than one var in it.
                        for col in precip_cols:
                            relcols = [k for k in precip_cols if col not in k]
                            # Get list of all vars other than col.
                            ds[col].attrs["ancillary_variables"] = " ".join(relcols)

                # hurs: relative humidity
                if "relative_humidity_set_1" in ds.keys():
                    ds = ds.rename({"relative_humidity_set_1": "hurs"})
                    # Set attributes
                    ds["hurs"].attrs["long_name"] = "relative_humidity"
                    ds["hurs"].attrs["standard_name"] = "relative_humidity"
                    ds["hurs"].attrs["units"] = "percent"

                    # If QA/QC column exists
                    if "relative_humidity_set_1_qc" in ds.keys():
                        flagvals = ds["relative_humidity_set_1_qc"].values.tolist()[0]
                        flagvals = [x for x in flagvals if pd.isnull(x) == False]
                        flagvals = list(np.unique(flagvals))  # Get unique values
                        flagvals = list(
                            np.unique(
                                [
                                    split_item
                                    for item in flagvals
                                    for split_item in str(item).split(";")
                                ]
                            )
                        )
                        # Split any rows with multiple flags and run unique again.

                        if flagvals == ["nan"]:
                            # This should not occur, but leave in as additional check.
                            ds = ds.drop("relative_humidity_set_1_qc")
                        else:
                            ds = ds.rename({"relative_humidity_set_1_qc": "hurs_qc"})
                            # Reformat to match csv.
                            flagvals = [ele.replace(".0", "") for ele in flagvals]
                            # Generate values from unique values from dataset.
                            ds["hurs_qc"].attrs["flag_values"] = flagvals
                            ds["hurs_qc"].attrs[
                                "flag_meanings"
                            ] = "See QA/QC csv for network."
                            # List other variables associated with variable (QA/QC)
                            ds["hurs"].attrs["ancillary_variables"] = "hurs_qc"

                # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)

                if "solar_radiation_set_1" in ds.keys():
                    # Already in w/m2, no need to convert units.
                    # If column exists, rename.
                    ds = ds.rename({"solar_radiation_set_1": "rsds"})

                    # Set attributes
                    ds["rsds"].attrs["long_name"] = "solar_radiation"
                    ds["rsds"].attrs[
                        "standard_name"
                    ] = "surface_downwelling_shortwave_flux_in_air"
                    ds["rsds"].attrs["units"] = "W m-2"

                # rsds: QA/QC flags
                if "solar_radiation_set_1_qc" in ds.keys():
                    flagvals = ds["solar_radiation_set_1_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("solar_radiation_set_1_qc")
                    else:
                        ds = ds.rename({"solar_radiation_set_1_qc": "rsds_qc"})
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        # Generate values from unique values from dataset.
                        ds["rsds_qc"].attrs["flag_values"] = flagvals
                        ds["rsds_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["rsds"].attrs["ancillary_variables"] = "rsds_qc"

                # sfcWind : wind speed (m/s)
                if "wind_speed_set_1" in ds.keys():
                    # Data already in m/s.
                    ds = ds.rename({"wind_speed_set_1": "sfcWind"})
                    ds["sfcWind"].attrs["long_name"] = "wind_speed"
                    ds["sfcWind"].attrs["standard_name"] = "wind_speed"
                    ds["sfcWind"].attrs["units"] = "m s-1"
                    ds["sfcWind"].attrs[
                        "comment"
                    ] = "Method of wind speed calculation varies within network, with 2-minute mean as CWOP sampling standard."
                    # (Method of calculation may vary and is unknown source by source.)
                    # See: https://weather.gladstonefamily.net/CWOP_Guide.pdf

                if "wind_speed_set_1_qc" in ds.keys():
                    flagvals = ds["wind_speed_set_1_qc"].values.tolist()[0]
                    flagvals = [
                        x for x in flagvals if pd.isnull(x) == False
                    ]  # Remove nas
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("wind_speed_set_1_qc")
                    else:
                        # Otherwise, rename and reformat.
                        ds = ds.rename({"wind_speed_set_1_qc": "sfcWind_qc"})
                        # Reformat to match csv.
                        flagvals = [ele.replace(".0", "") for ele in flagvals]

                        ds["sfcWind_qc"].attrs["flag_values"] = flagvals
                        ds["sfcWind_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["sfcWind"].attrs["ancillary_variables"] = "sfcWind_qc"

                # sfcWind_dir: wind direction
                if "wind_direction_set_1" in ds.keys():
                    # No conversions needed, do not make raw column.
                    ds = ds.rename({"wind_direction_set_1": "sfcWind_dir"})
                    ds["sfcWind_dir"].attrs["long_name"] = "wind_direction"
                    ds["sfcWind_dir"].attrs["standard_name"] = "wind_from_direction"
                    ds["sfcWind_dir"].attrs["units"] = "degrees_clockwise_from_north"
                    ds["sfcWind_dir"].attrs[
                        "comment"
                    ] = "Wind direction is defined by the direction that the wind is coming from (i.e., a northerly wind originates in the north and blows towards the south)."

                if "wind_direction_set_1_qc" in ds.keys():
                    flagvals = ds["wind_direction_set_1_qc"].values.tolist()[0]
                    flagvals = [x for x in flagvals if pd.isnull(x) == False]
                    flagvals = list(np.unique(flagvals))  # Get unique values
                    flagvals = list(
                        np.unique(
                            [
                                split_item
                                for item in flagvals
                                for split_item in str(item).split(";")
                            ]
                        )
                    )
                    # Split any rows with multiple flags and run unique again.

                    if flagvals == ["nan"]:
                        # This should not occur, but leave in as additional check.
                        ds = ds.drop("wind_direction_set_1_qc")
                    else:
                        # Otherwise, rename and reformat.
                        ds = ds.rename({"wind_direction_set_1_qc": "sfcWind_dir_qc"})
                        # Reformat to match csv.

                        flagvals = [ele.replace(".0", "") for ele in flagvals]
                        ds["sfcWind_dir_qc"].attrs["flag_values"] = flagvals
                        ds["sfcWind_dir_qc"].attrs[
                            "flag_meanings"
                        ] = "See QA/QC csv for network."
                        # List other variables associated with variable (QA/QC)
                        ds["sfcWind_dir"].attrs[
                            "ancillary_variables"
                        ] = "sfcWind_dir_qc"

                # Other variables: rename to match format
                if "altimeter_set_1" in ds.keys():
                    ds = ds.rename({"altimeter_set_1": "ps_altimeter"})
                if "altimeter_set_1_qc" in ds.keys():
                    ds = ds.rename({"altimeter_set_1_qc": "ps_altimeter_qc"})
                if "dew_point_temperature_set_1d" in ds.keys():
                    ds = ds.rename({"dew_point_temperature_set_1d": "tdps_derived"})
                if "pressure_set_1d" in ds.keys():
                    ds = ds.rename({"pressure_set_1d": "ps_derived"})

                if "ps_altimeter" in ds.keys():
                    ds["ps_altimeter"].attrs["long_name"] = "altimeter"
                    ds["ps_altimeter"].attrs["units"] = "Pa"
                    ds["ps_altimeter"].attrs["ancillary_variables"] = "ps"

                if "tdps_derived" in ds.keys():
                    ds["tdps_derived"] = calc_clean._unit_degC_to_K(ds["tdps_derived"])
                    ds["tdps_derived"].attrs[
                        "long_name"
                    ] = "derived_dew_point_temperature"
                    ds["tdps_derived"].attrs["units"] = "degree_Kelvin"
                    ds["tdps_derived"].attrs[
                        "comment"
                    ] = "Derived by Synoptic. Converted from Celsius to Kelvin."

                if "ps_derived" in ds.keys():
                    ds["ps_derived"].attrs["long_name"] = "derived_station_pressure"
                    ds["ps_derived"].attrs["units"] = "Pa"
                    ds["ps_derived"].attrs["comment"] = "Derived by Synoptic."

                # Drop repeated Station_ID column.
                ds = ds.drop("Station_ID")

                # Quality control: if any variable is completely empty, drop it.
                for key in ds.keys():
                    try:
                        if np.isnan(ds[key].values).all():
                            if "elevation" not in key:
                                # Don't drop elevation if NaN
                                print(f"Dropping {key}")
                                ds = ds.drop(key)
                    except Exception as e:
                        # Add to handle errors for unsupported data types
                        continue

                # For QA/QC flags, replace np.nan with "nan" to avoid h5netcdf overwrite to blank.
                for key in ds.keys():
                    if "qc" in key:
                        # Coerce all values in key to string.
                        ds[key] = ds[key].astype(str)

                # Reorder variables
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
                actual_order = [i for i in desired_order if i in list(ds.keys())]
                rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order]
                # Retain rest of variables at the bottom.
                new_index = actual_order + rest_of_vars
                ds = ds[new_index]

            except Exception as e:
                print(traceback.format_exc())
                errors["File"].append(stat_files)
                errors["Time"].append(end_api)
                errors["Error"].append(f"Error in ds set-up: {e}")
                continue

            # Write station file to netcdf.
            if len(ds.keys()) == 0:
                # skip station if the entire dataset will be empty because no data is observed
                print(
                    f"{station_id} has no data for all meteorological variables of interest throughout its current reporting; station not cleaned."
                )
                errors["File"].append(station_id)
                errors["Time"].append(end_api)
                errors["Error"].append("Station reports all nan meteorological data.")
                continue

            elif len(ds.time) == 0:
                # this should not be necessary, but placing for testing
                print(
                    f"No data for {station_id} during v1 period (1/1980 - 8/2022), station not cleaned."
                )
                errors["File"].append(stat_files)
                errors["Time"].append(end_api)
                errors["Error"].append(
                    "No data for this station during v1 period (1/1980 - 8/2022)"
                )
                continue

            else:
                try:
                    filename = station_id + ".nc"
                    filepath = cleandir + filename  # Write file path

                    # Write locally
                    ds.to_netcdf(path="temp/temp.nc", engine="netcdf4")
                    # Previously engine="h5netcdf" but was breaking on macOS

                    # Push file to AWS with correct file name.
                    s3.Bucket(BUCKET_NAME).upload_file("temp/temp.nc", filepath)

                    print(f"Saving {filename} with dims {ds.dims}")
                    ds.close()  # Close dataframe.

                except Exception as e:
                    print(e)
                    errors["File"].append(stat_files)
                    errors["Time"].append(end_api)
                    errors["Error"].append(
                        f"Error saving ds as .nc file to AWS bucket: {e}"
                    )
                    continue

        # Write the list of removed variables to csv for future reference. Keep these in one centralized "removedvars.csv" file that gets appended to
        # Read in existing removedvars.csv from AWS folder, if it already exists.
        # If removedvars.csv already in file list, save it to append to.
        key = cleandir + "removedvars.csv"
        removedvars = pd.DataFrame(removedvars, columns=["Variable"])
        try:
            test = s3.Bucket(BUCKET_NAME).Object(key).get()
            removedvarscsv = pd.read_csv(test["Body"])
            removedvars = pd.concat(
                [removedvarscsv, removedvars], axis=0, ignore_index=True
            )
            removedvars = removedvars[["Variable"]].drop_duplicates(ignore_index=True)

            csv_buffer = StringIO()
            removedvars.to_csv(csv_buffer)
            content = csv_buffer.getvalue()
            s3_cl.put_object(
                Bucket=BUCKET_NAME, Body=content, Key=cleandir + "removedvars.csv"
            )

        except botocore.exceptions.ClientError as ex:
            if ex.response["Error"]["Code"] == "NoSuchKey":
                # If removed vars doesn't already exist, save.
                # Save to AWS
                csv_buffer = StringIO()
                removedvars.to_csv(csv_buffer)
                content = csv_buffer.getvalue()
                s3_cl.put_object(
                    Bucket=BUCKET_NAME, Body=content, Key=cleandir + "removedvars.csv"
                )

            else:
                # If some other error occurs.
                errors["File"].append("Whole network error")
                errors["Time"].append(end_api)
                errors["Error"].append(ex.response)

    # Write errors.csv
    finally:
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(
            Bucket=BUCKET_NAME,
            Body=content,
            Key=cleandir + f"errors_{network}_{end_api}.csv",
        )

    return None


# Run functions
if __name__ == "__main__":
    network = "CWOP"
    rawdir, cleandir, qaqcdir = get_file_paths(network)
    get_qaqc_flags(token=config.token, qaqcdir=qaqcdir, network=network)
    clean_madis(
        rawdir, cleandir, network, cwop_letter=None
    )  # if cwop_letter is not None, argument must be passed as a string


# List of MADIS network names
# 'CAHYDRO', 'CDEC', 'CNRFC', 'CRN', 'CWOP', 'HNXWFO', 'HOLFUY', 'HPWREN',
# 'LOXWFO', 'MAP', 'MTRWFO', 'NCAWOS', 'NOS-NWLON', 'NOS-PORTS',
# 'RAWS', 'SGXWFO', 'SHASAVAL', 'VCAPCD', 'HADS'

# Note: CWOP, RAWS, and HADS will take a long time to run to complete full network clean
