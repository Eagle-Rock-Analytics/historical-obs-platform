"""
CW3E_clean.py

This script performs data cleaning for CW3E stations pulled through Scripps Institution for Oceanography for
ingestion into the Historical Observations Platform.

Approach
--------
(1) Read through variables and drop unnecessary variables
(2) Converts station metadata to standard format, with unique identifier
(3) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(4) Converts missing data to standard format
(5) Merge files by station, and outputs cleaned variables as a single .nc file for each station in an individual network.

Functions
---------
- clean_cw3e: Cleans CW3E data. 

Intended Use
------------
Cleaned data for an individual network, priority variables, all times. Organized by station per year as .nc file.

Notes
-----
1. One station, Lower Bath House, has an entirely different data structure and file format and is not included in MADIS metadata.
Data is available from 04/23/2021-11/18/2021. At this time this script does not clean this file, but it could be included in future rounds of cleaning.
2. This dataset does not have QA/QC flags.
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd
import boto3
from io import BytesIO, StringIO
import dask.dataframe as dd
import traceback
import warnings

# Optional: Silence pandas' future warnings about regex (not relevant here)
warnings.filterwarnings(action="ignore", category=FutureWarning)

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


def clean_cw3e(rawdir: str, cleandir: str):
    """Clean CW3E data.

    Paramters
    ---------
    rawdir : str
        path to raw data bucket
    cleandir : str
        path to cleaned data bucket

    Returns
    -------
    None
    """

    network = "CW3E"

    # Set up error handling.
    errors = {"File": [], "Time": [], "Error": []}
    # Set end time to be current time at beginning of download: for error handling csv.
    end_api = datetime.now().strftime("%Y%m%d%H%M")
    # For attributes of netCDF file.
    timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    # Specify columns to remove
    removecols = [
        "Datalogger ID",
        "Wind Direction Standard Deviation (degrees)",
        "Vector Wind Speed (m/s)",
        "Battery Voltage (volts)",
        "Maximum Wind Speed (m/s)",
        "Soil Temperature (C) 5cm",
        "Soil Temperature (C) 10cm",
        "Soil Temperature (C) 15cm",
        "Soil Temperature (C) 20cm",
        "Soil Temperature (C) 50cm",
        "Soil Temperature (C) 100cm",
        "Soil Reflectometer Output Period (usec) 5cm",
        "Soil Reflectometer Output Period (usec) 10cm",
        "Soil Reflectometer Output Period (usec) 15cm",
        "Soil Reflectometer Output Period (usec) 20cm",
        "Soil Reflectometer Output Period (usec) 50cm",
        "Soil Reflectometer Output Period (usec) 100cm",
        "Soil Reflectometer Output Period (%) 5cm",
        "Soil Reflectometer Output Period (%) 10cm",
        "Soil Reflectometer Output Period (%) 15cm",
        "Soil Reflectometer Output Period (%) 20cm",
        "Soil Reflectometer Output Period (%) 50cm",
        "Soil Reflectometer Output Period (%) 100cm",
    ]

    # Specify columns to use if there is no DataFormat file
    default_cols = [
        "Datalogger ID",
        "Year (end time of average)",
        "Julian Day (end time of average)",
        "HoursMinutes (end time of average)",
        "Pressure (mb)",
        "Temperature (C)",
        "Relative Humidity (%)",
        "Scalar Wind Speed (m/s)",
        "Vector Wind Speed (m/s)",
        "Wind Direction (degrees)",
        "Wind Direction Standard Deviation (degrees)",
        "Solar Radiation (W/m^2)",
        "Battery Voltage (volts)",
        "Precipitation (mm)",
        "Maximum Wind Speed (m/s)",
        "Soil Temperature (C) 5cm",
        "Soil Temperature (C) 10cm",
        "Soil Temperature (C) 15cm",
        "Soil Temperature (C) 20cm",
        "Soil Temperature (C) 50cm",
        "Soil Temperature (C) 100cm",
        "Soil Reflectometer Output Period (usec) 5cm",
        "Soil Reflectometer Output Period (usec) 10cm",
        "Soil Reflectometer Output Period (usec) 15cm",
        "Soil Reflectometer Output Period (usec) 20cm",
        "Soil Reflectometer Output Period (usec) 50cm",
        "Soil Reflectometer Output Period (usec) 100cm",
    ]

    date_parser = lambda x, y, z: datetime.strptime(f"{x}.{y}.{z}", "%Y.%j.%H%M")

    try:
        # Get files
        files = []
        for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=rawdir):
            file = str(item.key)
            files += [file]

        # Get station file and read in metadata.
        station_file = [file for file in files if "stationlist_" in file]
        obj = s3_cl.get_object(Bucket=BUCKET_NAME, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj["Body"].read()))
        stations = station_file["STID"].dropna().tolist()
        stations = [station.replace("C3", "") for station in stations]

        # Remove error, station files, only some stations have a data format file
        format_files = [file for file in files if "_DataFormat.txt" in file]
        # all valid stations have a readme file
        readme_files = [file for file in files if "_README.txt" in file]
        files = [file for file in files if "README" not in file]
        files = [file for file in files if "stationlist" not in file]
        files = [file for file in files if "error" not in file]

    except Exception as e:
        # If unable to read files from rawdir, break function.
        print(f"Error: {e.args[0]}. Code line: {str(traceback.extract_stack()[-1][1])}")
        errors["File"].append("Whole network")
        errors["Time"].append(end_api)
        errors["Error"].append(f"Whole network error: {e}")

    else:
        # If files read successfully, continue.
        for station in stations:
            # for station in ["CAT"]: # uncomment for a specific station
            # Full network clean
            print(f"Parsing: {station}")
            try:
                # If station does not have a README file, automatically skip - there will be no data on AWS
                read_stn = f"{rawdir}{station}_README.txt"
                if read_stn not in readme_files:
                    print(f"No raw data found for {station} on AWS.")
                    errors["File"].append(station)
                    errors["Time"].append(end_api)
                    errors["Error"].append("No raw data found on AWS -- not cleaned")
                    continue

                # Get metadata from .txt file
                obj = s3_cl.get_object(
                    Bucket=BUCKET_NAME, Key=rawdir + f"{station}_README.txt"
                )
                with StringIO(obj["Body"].read().decode()) as f:
                    for line in f.readlines():
                        if "Lat" in line:
                            latval = float(line.split(": ")[1])
                        if "Lon" in line:
                            lonval = float(line.split(": ")[1])
                        if "Elev" in line:
                            elevval = line.split(": ")[1]
                            if "m" in elevval:
                                elevval = float(elevval.replace("m", ""))
                            elif "f" in elevval:
                                elevval = float(elevval.replace("f", ""))
                                elevval = calc_clean._unit_elev_ft_to_m(elevval)

                # Get column names from DataFormat.txt file if available, if not use default columns list
                if station == "DLA" or station == "CAT" or station == "FRC":
                    colnames = default_cols
                else:
                    obj = s3_cl.get_object(
                        Bucket=BUCKET_NAME, Key=rawdir + f"{station}_DataFormat.txt"
                    )
                    dataformat = pd.read_csv(
                        BytesIO(obj["Body"].read()),
                        sep=":",
                        skipinitialspace=True,
                        names=["No", "ColName"],
                    )
                    colnames = dataformat["ColName"].tolist()[1:]

                # Station full name is not in README so we pull it in from MADIS metadata.
                station_metadata = station_file.loc[
                    station_file["STID"].str.replace("C3", "") == station
                ]
                station_id = "CW3E_" + str(station)

                # Get files per station
                station_files = [
                    file
                    for file in files
                    if file.startswith(f"1_raw_wx/CW3E/{station.lower()}")
                ]

                # Get years - note this is hardcoded and should be updated for the update script.
                years = ["19", "20", "21", "22"]

                for year in years:
                    # Get files for station for year.
                    year_files = [file for file in station_files if year == file[-9:-7]]
                    file_count = len(year_files)

                    if file_count == 0:
                        print(
                            f"{station} does not have data for 20{year}, moving to next available year of data."
                        )
                        continue
                    else:
                        # Read in data.
                        try:
                            print(
                                f"{station} has {file_count} files for 20{year}, cleaning in progress..."
                            )
                            df_stat = dd.read_csv(
                                f"s3://wecc-historical-wx/1_raw_wx/CW3E/{station.lower()}{year}*m",
                                names=colnames,
                                na_values=[-99999],
                                parse_dates={
                                    "time": [
                                        "Year (end time of average)",
                                        "Julian Day (end time of average)",
                                        "HoursMinutes (end time of average)",
                                    ]
                                },
                                date_parser=date_parser,
                                dtype={
                                    "Temperature (C)": "float64",
                                    "Pressure (mb)": "float64",
                                    "Solar Radiation (W/m^2)": "float64",
                                    "Relative Humidity (%)": "float64",
                                    "Precipitation (mm)": "float64",
                                    "Scalar Wind Speed (m/s)": "float64",
                                    "Wind Direction (degrees)": "float64",
                                },
                                assume_missing=True,
                            )

                        except OSError as e:
                            # If year has data, but fails to read-in
                            errors["File"].append(station_id)
                            errors["Time"].append(end_api)
                            errors["Error"].append(f"Error in df set-up: {e}")
                            continue

                        try:
                            # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022.
                            df_stat = df_stat.loc[
                                (df_stat["time"] < "2022-09-01")
                                & (df_stat["time"] > "1979-12-31")
                            ]

                            df_stat = df_stat.compute()

                            # Drop any columns that only contain NAs.
                            df_stat = df_stat.dropna(axis=1, how="all")

                            # Drop unnecessary columns
                            try:
                                cols_to_keep = [
                                    col
                                    for col in df_stat.columns
                                    if col not in removecols
                                ]
                                df_stat = df_stat[
                                    df_stat.columns.intersection(cols_to_keep)
                                ]
                            except Exception as e:
                                print(
                                    f"Dropping unnecessary column error for {station} -- check."
                                )
                                errors["File"].append(station)
                                errors["Time"].append(end_api)
                                errors["Error"].append("Unnecessary column drop error.")
                                continue

                            # Sort by time and remove any overlapping timestamps.
                            df_stat = df_stat.sort_values(by="time")
                            df_stat = df_stat.drop_duplicates()

                            ds = df_stat.to_xarray()

                            # Update global attributes
                            ds = ds.assign_attrs(title=f"{network} cleaned")
                            ds = ds.assign_attrs(
                                institution="Eagle Rock Analytics / Cal Adapt"
                            )
                            ds = ds.assign_attrs(source="")
                            ds = ds.assign_attrs(
                                history=f"CW3E_clean.py script run on {timestamp} UTC"
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
                                station_name=station_metadata["NAME"].values[0]
                            )

                            # Sensor heights - Final confirmation is still waiting on CW3E response.
                            ds = ds.assign_attrs(barometer_elev_m=np.nan)
                            ds = ds.assign_attrs(pyranometer_height_m=np.nan)
                            ds = ds.assign_attrs(wind_vane_height_m=np.nan)
                            ds = ds.assign_attrs(anemometer_height_m=np.nan)
                            ds = ds.assign_attrs(thermometer_height_m=np.nan)
                            ds = ds.assign_attrs(humidity_height_m=np.nan)
                            ds = ds.assign_attrs(rain_gauge_height_m=np.nan)

                            # Keep count of how many files merged per station.
                            ds = ds.assign_attrs(raw_files_merged=file_count)

                            # Add dimensions: station ID and time.
                            # Swap index with time.
                            ds = ds.set_coords("time").swap_dims({"index": "time"})
                            ds = ds.assign_coords(id=str(station_id))
                            # Add station_id as index.
                            ds = ds.expand_dims("id")
                            # Drop station_id variable and index coordinate.
                            ds = ds.drop_vars(("index"))
                            # Rename id to station_id.
                            ds = ds.rename({"id": "station"})

                            # Update dimensions and coordinates

                            # Add coordinates: latitude and longitude.
                            lat = np.asarray([latval] * len(ds["time"]))
                            lat.shape = (1, len(ds["time"]))

                            lon = np.asarray([lonval] * len(ds["time"]))
                            lon.shape = (1, len(ds["time"]))

                            # reassign lat and lon as coordinates
                            ds = ds.assign_coords(
                                lat=(["station", "time"], lat),
                                lon=(["station", "time"], lon),
                            )

                            # If any observation is missing lat or lon coordinates, drop these observations.
                            if np.count_nonzero(np.isnan(ds["lat"])) != 0:
                                ds = ds.where(~np.isnan(ds["lat"]))

                            if np.count_nonzero(np.isnan(ds["lon"])) != 0:
                                ds = ds.where(~np.isnan(ds["lon"]))

                            # Add variable: elevation (in meters)
                            elev = np.asarray([elevval] * len(ds["time"]))
                            elev.shape = (1, len(ds["time"]))
                            ds["elevation"] = (["station", "time"], elev)

                            # Update dimension and coordinate attributes.

                            # Update attributes.
                            # Remove timezone data from string (to match other networks)
                            ds["time"] = pd.to_datetime(ds["time"])
                            ds["time"] = pd.to_datetime(ds["time"], unit="ns")
                            ds["time"].attrs["long_name"] = "time"
                            ds["time"].attrs["standard_name"] = "time"
                            ds["time"].attrs["comment"] = "Converted from PST to UTC."

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
                            ds["elevation"].attrs[
                                "standard_name"
                            ] = "height_above_mean_sea_level"
                            ds["elevation"].attrs["long_name"] = "station_elevation"
                            ds["elevation"].attrs["units"] = "meters"
                            # Define which direction is positive
                            ds["elevation"].attrs["positive"] = "up"

                            # Update variable attributes and do unit conversions

                            # tas: air surface temperature (K)
                            if "Temperature (C)" in ds.keys():
                                ds["tas"] = calc_clean._unit_degC_to_K(
                                    ds["Temperature (C)"]
                                )
                                ds = ds.drop("Temperature (C)")

                                ds["tas"].attrs["long_name"] = "air_temperature"
                                ds["tas"].attrs["standard_name"] = "air_temperature"
                                ds["tas"].attrs["units"] = "degree_Kelvin"
                                ds["tas"].attrs[
                                    "comment"
                                ] = "Converted from Celsius to Kelvin."

                            # ps: surface air pressure (Pa)
                            if "Pressure (mb)" in ds.keys():
                                # If barometric pressure available
                                # Convert from inHg to PA
                                ds["psl"] = calc_clean._unit_pres_hpa_to_pa(
                                    ds["Pressure (mb)"]
                                )
                                ds = ds.drop("Pressure (mb)")

                                # Set attributes
                                ds["psl"].attrs["long_name"] = "station_air_pressure"
                                ds["psl"].attrs["standard_name"] = "air_pressure"
                                ds["psl"].attrs["units"] = "Pa"
                                ds["psl"].attrs["comment"] = "Converted from mb."

                            # tdps: dew point temperature (K)
                            # Not available in this dataset.

                            # pr: precipitation
                            if "Precipitation (mm)" in ds.keys():
                                ds = ds.rename({"Precipitation (mm)": "pr"})
                                ds["pr"].attrs[
                                    "long_name"
                                ] = "precipitation_accumulation"
                                ds["pr"].attrs["units"] = "mm/2 minutes"
                                ds["pr"].attrs["comment"] = "Accumulated precipitation."

                            # hurs: relative humidity (%)
                            if "Relative Humidity (%)" in ds.keys():
                                # Already in %, no need to convert units.
                                ds = ds.rename({"Relative Humidity (%)": "hurs"})
                                # Set attributes
                                ds["hurs"].attrs["long_name"] = "relative_humidity"
                                ds["hurs"].attrs["standard_name"] = "relative_humidity"
                                ds["hurs"].attrs["units"] = "percent"

                            # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
                            if "Solar Radiation (W/m^2)" in ds.keys():
                                # Already in w/m2, no need to convert units.
                                # If column exists, rename.
                                ds = ds.rename({"Solar Radiation (W/m^2)": "rsds"})

                                # Set attributes
                                ds["rsds"].attrs["long_name"] = "solar_radiation"
                                ds["rsds"].attrs[
                                    "standard_name"
                                ] = "surface_downwelling_shortwave_flux_in_air"
                                ds["rsds"].attrs["units"] = "W m-2"

                            # sfcWind : wind speed (m/s)
                            if "Scalar Wind Speed (m/s)" in ds.keys():
                                # Data originally in mph.
                                ds = ds.rename({"Scalar Wind Speed (m/s)": "sfcWind"})
                                ds["sfcWind"].attrs["long_name"] = "wind_speed"
                                ds["sfcWind"].attrs["standard_name"] = "wind_speed"
                                ds["sfcWind"].attrs["units"] = "m s-1"

                            # sfcWind_dir: wind direction
                            if "Wind Direction (degrees)" in ds.keys():
                                # No conversions needed, do not make raw column.
                                ds = ds.rename(
                                    {"Wind Direction (degrees)": "sfcWind_dir"}
                                )
                                ds["sfcWind_dir"].attrs["long_name"] = "wind_direction"
                                ds["sfcWind_dir"].attrs[
                                    "standard_name"
                                ] = "wind_from_direction"
                                ds["sfcWind_dir"].attrs[
                                    "units"
                                ] = "degrees_clockwise_from_north"
                                ds["sfcWind_dir"].attrs[
                                    "comment"
                                ] = "Wind direction is defined by the direction that the wind is coming from (i.e., a northerly wind originates in the north and blows towards the south)."

                            # Quality control: if any variable is completely empty, drop it.
                            for key in ds.keys():
                                try:
                                    if np.isnan(ds[key].values).all():
                                        if "elevation" not in key:  # Exclude elevation
                                            print(f"Dropping {key}")
                                            ds = ds.drop(key)
                                except:
                                    # Add to handle errors for unsupported data types
                                    continue

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
                            actual_order = [
                                i for i in desired_order if i in list(ds.keys())
                            ]
                            # Retain rest of variables at the bottom.
                            rest_of_vars = [
                                i for i in list(ds.keys()) if i not in desired_order
                            ]
                            new_index = actual_order + rest_of_vars
                            ds = ds[new_index]

                            # Write station file to netcdf.
                            if ds is None:
                                # Should be caught in error handling above, but add in case.
                                print(f"{station + year} not saved.")
                                errors["File"].append(station + year)
                                errors["Time"].append(end_api)
                                errors["Error"].append("File has no data.")
                                continue

                            else:
                                try:
                                    filename = station_id + "_" + year + ".nc"
                                    # Make file name (stationID+year)
                                    filepath = cleandir + filename  # Write file path

                                    # Write locally
                                    ds.to_netcdf(path="temp/temp.nc", engine="netcdf4")

                                    # Push file to AWS with correct file name.
                                    s3.Bucket(BUCKET_NAME).upload_file(
                                        "temp/temp.nc", filepath
                                    )

                                    print(f"Saving {filename} with dims {ds.dims}")
                                    ds.close()  # Close dataframe.

                                except Exception as e:
                                    print(
                                        f"Error: {e.args[0]}. Code line: {str(traceback.extract_stack()[-1][1])}"
                                    )
                                    errors["File"].append(filename)
                                    errors["Time"].append(end_api)
                                    errors["Error"].append(
                                        f"Error in saving ds as .nc file to AWS: {e}"
                                    )
                                    continue

                        except Exception as e:
                            print(
                                f"Error: {e.args[0]}. Code line: {str(traceback.extract_stack()[-1][1])}"
                            )
                            errors["File"].append(station + year)  # Save ID of station.
                            errors["Time"].append(end_api)
                            errors["Error"].append(f"Error in ds set-up: {e}")
                            continue

            except Exception as e:
                print(
                    f"Error: {e.args[0]}. Code line: {str(traceback.extract_stack()[-1][1])}"
                )
                errors["File"].append(station)  # Save ID of station.
                errors["Time"].append(end_api)
                errors["Error"].append(
                    f"Error in full clean set-up/GetObject set-up for {station}: {e}"
                )
                continue

        # Save the list of removed variables to AWS
        removedvars = pd.DataFrame(removecols, columns=["Variable"])
        csv_buffer = StringIO()
        removedvars.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(
            Bucket=BUCKET_NAME, Body=content, Key=cleandir + "removedvars.csv"
        )

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
    rawdir, cleandir, qaqcdir = get_file_paths("CW3E")
    clean_cw3e(rawdir, cleandir)
