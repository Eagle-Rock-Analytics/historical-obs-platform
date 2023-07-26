"""
This script performs data cleaning for networks pulled through Synoptic API for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables and drop unnecessary variables
(2) Converts station metadata to standard format, with unique identifier
(3) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(4) Converts missing data to standard format
(5) Tracks existing qa/qc flag for review
(6) Merge files by station, and outputs cleaned variables as a single .nc file for each station in an individual network.
Inputs: Raw data for the network's stations, with each csv file representing a station.
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
Reference: https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf

Note: QAQC flags and removed variable lists both formatted and uploaded manually. Last update Nov 9 2022.
"""

# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime, date, timedelta
import re
import numpy as np
import warnings
warnings.filterwarnings(action = 'ignore', category = FutureWarning) # Optional: Silence pandas' future warnings about regex (not relevant here)
import pandas as pd
import boto3
from io import BytesIO, StringIO
import random
from cleaning_helpers import get_file_paths
# To be able to open xarray files from S3, h5netcdf must also be installed, but doesn't need to be imported.


# Import calc_clean.py.
try:
    import calc_clean
except:
    print("Error importing calc_clean.py")
    pass

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

# Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

# Set up directory to save files temporarily, if it doesn't already exist.
try:
    os.mkdir('temp') # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass


## FUNCTION: Clean SCAN and SNOTEL data.
## Input:
## bucket_name: name of AWS bucket.
## rawdir: path to where raw data is saved as .csv files, with each file representing a station's records from download start date to present (by default 01-01-1980).
## cleandir: path to where cleaned files should be saved.
## Output: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
## Note: files are already filtered by bbox in the pull step, so additional geospatial filtering is not required here.

## Function: take heads, read csv into pandas db and clean.
def clean_scansnotel(rawdir, cleandir):

    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
    timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.

    try:
        # Get files
        files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir):
            file = str(item.key)
            files += [file]

        # Get station file and read in metadata.
        station_file = [file for file in files if 'station' in file]
        if len(station_file)>1: # If more than one file returned
            station_file = [file for file in station_file if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj['Body'].read()))

        # Remove error, station files
        files = [file for file in files if '.csv' in file]
        files = [file for file in files if 'station' not in file]
        files = [file for file in files if 'error' not in file]

        vars_to_remove = ['TAVG', 'RHUMV', 'SRADV', 'SRADT', 'WDIRV', 'WSPDV'] # Variables to remove
        cols_to_remove = ['{}_flag'.format(elem) for elem in vars_to_remove] + ['{}_value'.format(elem) for elem in vars_to_remove] + ['{}_time'.format(elem) for elem in vars_to_remove] # Associated columns for each variable

        # Get list of station IDs from filename and clean.
        ids = list()
        for file in files:
            id = file.split("/")[-1] # Remove leading folders
            id = id.split("_")[-1] # Remove leading prefixes (for mult files)
            id = re.sub('.csv', "", id)
            if id not in ids:
                ids.append(id)

    except Exception as e: # If unable to read files from rawdir, break function.
        # print(e)
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append("Whole network error: {}".format(e))

    else: # If files read successfully, continue.
        for i in ids: # For each station (full run)
        # for i in random.sample(ids,2): # Subsample for testing errors
            df_stat = None # Initialize merged df.
            try:
                stat_files = [k for k in files if i in k] # Get list of files with station ID in them.

                station_id = "{}_".format(network)+i.split(":")[0] # Save file ID, keeping STID from triplet.

                if not stat_files: # If ID has no files downloaded
                    print('No raw data found for {} on AWS.'.format(station_id))
                    errors['File'].append(station_id)
                    errors['Time'].append(end_api)
                    errors['Error'].append("No raw data found for station on AWS.")
                    continue # Skip ID.

                # Get metadata attributes from station list.
                station_metadata = station_file.loc[station_file['stationTriplet'] == i] # Get metadata for station.
                station_metadata = station_metadata.iloc[0]

                for file in stat_files: # For each station
                # Read in CSV, removing header.
                    try:
                        obj = s3_cl.get_object(Bucket=bucket_name, Key=file)
                        df = pd.read_csv(BytesIO(obj['Body'].read()), low_memory=False) # here we can either set engine = 'python' (slower) to suppress dtype warning, or remove but ignore it - it gets fixed down the line either way.

                        # Fix any NA mixed types
                        df.replace("NaN", np.nan, inplace=True)

                        # Drop any columns that only contain NAs.
                        df = df.dropna(axis = 1, how = 'all')

                        # Drop any unnecessary columns.
                        for col in cols_to_remove:
                            if col in df.columns:
                                df = df.drop(col, axis = 1)

                        # Resolve time differences between time and variable time observations.
                        # For each column, check if any value different than original time column
                        # Do so by setting to NA if identical, and then deleting all columns with all NAs.
                        timecols = [col for col in df.columns if 'time' in col]
                        timecols = [col for col in timecols if col != 'time'] # Remove 'time' column

                        for col in timecols:
                            df[col] = np.where(df['time']==df[col], np.nan, df[col]) # Assign nan to any value identical to the time column.

                        df = df.dropna(axis = 1, how = 'all') # if all columns are dropped, handled at last step, but could happen here too

                        if len([col for col in df.columns if 'time' in col]) > 1: # If more than one time column remains,
                            print("Conflicting time values: {}".format(list([col for col in df.columns if 'time' in col]))) # print warning
                            exit() # And kill code.
 
                        # Fix time format issues caused by "Timeout" errors.
                        df['time'] = pd.to_datetime(df['time'], errors = 'coerce') # Convert time to datetime and catch any incorrectly formatted columns.
                        df = df[pd.notnull(df['time'])] # Remove any rows where time is missing.

                        # Adjust station time by local time offset.
                        hours_adj = float(station_metadata['stationDataTimeZone']) # Get UTC offset from metadata
                        df['time'] = df['time'] - pd.Timedelta(hours = hours_adj)

                        # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022.
                        df = df.loc[(df['time']<'2022-09-01') & (df['time']>'1979-12-31')]

                        # If more than one file for station, merge files.
                        # will not be triggered for full clean, but set up to work if slices of data added in future.
                        # will require additional testing at this point.
                        if df_stat is None:
                            df_stat = df
                            del(df) # For memory.
                        else:
                            if len(stat_files)>1: # If there is more than one file per station
                                df_stat = pd.concat([df_stat, df], axis = 0, ignore_index= True)

                    except Exception as e: # Exceptions thrown during individual file read in.
                        print('Error in pandas df set-up')
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append('Error in pandas df set-up: {}'.format(e))
                        continue


                # Format joined station file.
                file_count = len(stat_files)

                # Fix multi-type columns
                # If column has QC in it, force to string.
                for b in df_stat.columns:
                    multitype = set(type(x).__name__ for x in df_stat[b])
                    if len(multitype)>1:
                        if 'flag' in b: # QC columns
                            df_stat[b] = df_stat[b].astype(str) # Coerce to string (to handle multiple QA/QC flags)
                        else:
                            print("Multitype error for column {} with data types {}. Please resolve".format(b, multitype)) # Code to flag novel exceptions, correct and add explicit handling above.
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("Multitype error for column {} with data types {}. Please resolve".format(b, multitype))
                            continue
                    else:
                        if 'flag' in b:
                            df_stat[b] = df_stat[b].astype(str) # Coerce QA/QC flag to string in all instances.


                # Fix issue with "nan" and nan causing comparison errors
                df_stat.replace("nan", np.nan, inplace=True)

                # Sort by time and remove any overlapping timestamps.
                df_stat = df_stat.sort_values(by = "time")
                df_stat = df_stat.drop_duplicates()

                # Move df to xarray object.
                ds = df_stat.to_xarray()
                del(df_stat)

                # Update global attributes
                ds = ds.assign_attrs(title = "{} cleaned".format(network))
                ds = ds.assign_attrs(institution = "Eagle Rock Analytics / Cal Adapt")
                ds = ds.assign_attrs(source = "")
                ds = ds.assign_attrs(history = "SCANSNOTEL_clean.py script run on {} UTC".format(timestamp))
                ds = ds.assign_attrs(comment = "Intermediate data product: may not have been subject to any cleaning or QA/QC processing")
                ds = ds.assign_attrs(license = "")
                ds = ds.assign_attrs(citation = "")
                ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")

                # Sensor heights - Setting thermometer & barometer as NaN pending response from USDA
                ds = ds.assign_attrs(thermometer_height_m = np.nan)
                ds = ds.assign_attrs(anemometer_height_m = 1.5)
                ds = ds.assign_attrs(barometer_elev_m = np.nan)

                # Add station metadata
                # Station name
                ds = ds.assign_attrs(station_name = station_metadata['name'])

                # Station sub-network
                if "SNTLT" in i:
                    subnetwork = "SNOTEL Lite"
                elif "SNTL" in i:
                    subnetwork = "SNOTEL"
                elif "CSCAN" in i:
                    subnetwork = "Tribal SCAN"
                elif "SCAN" in i:
                    subnetwork = "SCAN"
                ds = ds.assign_attrs(subnetwork = subnetwork)

                # Other station IDs - only keep if not NA
                if not np.isnan(station_metadata['huc']):
                    ds = ds.assign_attrs(huc = int(station_metadata['huc']))
                if not np.isnan(station_metadata['hud']):
                    ds = ds.assign_attrs(hud = int(station_metadata['hud']))
                if isinstance(station_metadata['shefId'], str):
                    ds = ds.assign_attrs(shefId = station_metadata['shefId'])

                ## commenting out the following section - it breaks 4 stations in terms of opening (sometimes) but it appears SCAN is missing actonId anyways?
                # if network == "SCAN":
                #     if isinstance(station_metadata['actonId'], str):
                #         ds = ds.assign_attrs(actonId = station_metadata['actonId'])

                ds = ds.assign_attrs(raw_files_merged = file_count) # Keep count of how many files merged per station.

                # Add dimensions: station ID and time.
                ds = ds.set_coords('time').swap_dims({'index': 'time'}) # Swap index with time.
                ds = ds.assign_coords(id = str(station_id))
                ds = ds.expand_dims("id") # Add station_id as index.
                ds = ds.drop_vars(("index")) # Drop station_id variable and index coordinate.
                ds = ds.rename({'id': 'station'}) # Rename id to station_id.

                # Update dimensions and coordinates

                # Add coordinates: latitude and longitude.
                lat = np.asarray([station_metadata['latitude']]*len(ds['time']))
                lat.shape = (1, len(ds['time']))

                lon = np.asarray([station_metadata['longitude']]*len(ds['time']))
                lon.shape = (1, len(ds['time']))

                # reassign lat and lon as coordinates
                ds = ds.assign_coords(lat = (["station","time"],lat), lon = (["station","time"], lon))

                # If any observation is missing lat or lon coordinates, drop these observations.
                if np.count_nonzero(np.isnan(ds['lat'])) != 0:
                    #print("Dropping missing lat values") # For testing.
                    ds = ds.where(~np.isnan(ds['lat']))

                if np.count_nonzero(np.isnan(ds['lon'])) != 0:
                    #print("Dropping missing lon values") # For testing.
                    ds = ds.where(~np.isnan(ds['lon']))


                # Add variable: elevation (in feet)
                elev = np.asarray([station_metadata['elevation']]*len(ds['time']))
                elev.shape = (1, len(ds['time']))
                ds['elevation'] = (['station', 'time'], elev)

                # Update dimension and coordinate attributes.

                # Update attributes.
                ds['time'].attrs['long_name'] = "time"
                ds['time'].attrs['standard_name'] = "time"
                ds['time'].attrs['comment'] = "Converted from local time zone UTC{} to UTC.".format(str(int(hours_adj)))

                # Station ID
                ds['station'].attrs['long_name'] = "station_id"
                ds['station'].attrs['comment'] = "Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network."

                # Latitude
                ds['lat'].attrs['long_name'] = "latitude"
                ds['lat'].attrs['standard_name'] = "latitude"
                ds['lat'].attrs['units'] = "degrees_north"

                # Longitude
                ds['lon'].attrs['long_name'] = "longitude"
                ds['lon'].attrs['standard_name'] = "longitude"
                ds['lon'].attrs['units'] = "degrees_east"

                # Elevation
                # Convert from feet to meters.
                ds['elevation'] = calc_clean._unit_elev_ft_to_m(ds['elevation']) # Convert feet to meters.
                ds['elevation'].attrs['standard_name'] = "height_above_mean_sea_level"
                ds['elevation'].attrs['long_name'] = "station_elevation"
                ds['elevation'].attrs['units'] = "meters"
                ds['elevation'].attrs['positive'] = "up" # Define which direction is positive
                ds['elevation'].attrs['comment'] = "Converted from feet to meters."

                # Update variable attributes and do unit conversions

                #Testing: Manually check values to see that they seem correctly scaled, no unexpected NAs.

                #tas: air surface temperature (K)
                if "TOBS_value" in ds.keys():
                    ds['tas'] = calc_clean._unit_degF_to_K(ds['TOBS_value'])
                    ds = ds.drop("TOBS_value")

                    ds['tas'].attrs['long_name'] = "air_temperature"
                    ds['tas'].attrs['standard_name'] = "air_temperature"
                    ds['tas'].attrs['units'] = "degree_Kelvin"

                    if "TOBS_flag" in ds.keys():
                        # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
                        ds = ds.rename({'TOBS_flag': 'tas_qc'})
                        ds['tas_qc'].attrs['flag_values'] = "V S E"
                        ds['tas_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['tas'].attrs['ancillary_variables'] = "tas_qc" # List other variables associated with variable (QA/QC)

                    ds['tas'].attrs['comment'] = "Converted from Fahrenheit to Kelvin."

                # ps: surface air pressure (Pa)

                # Here we only have barometric pressure (sea-level).
                if 'PRES_value' in ds.keys(): # If barometric pressure available
                    # Convert from inHg to PA
                    ds['psl'] = calc_clean._unit_pres_inHg_to_pa(ds['PRES_value'])
                    ds = ds.drop("PRES_value")

                    # Set attributes
                    ds['psl'].attrs['long_name'] = "barometric_air_pressure"
                    ds['psl'].attrs['standard_name'] = "air_pressure"
                    ds['psl'].attrs['units'] = "Pa"

                    if 'PRES_flag' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'PRES_flag': 'psl_qc'}) # this was previously set to TOBS_flag and tas_qc, in case this errors in the future
                        ds['psl_qc'].attrs['flag_values'] = "V S E"
                        ds['psl_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['psl'].attrs['ancillary_variables'] = "psl_qc" # List other variables associated with variable (QA/QC)

                    ds['psl'].attrs['comment'] = "Converted from inHg."


                # tdps: dew point temperature (K)
                # if raw dew point temperature observed, use that.
                if 'DPTP_value' in ds.keys():
                    ds['tdps'] = calc_clean._unit_degF_to_K(ds['DPTP_value'])
                    ds = ds.drop("DPTP_value")

                    # Set attributes for conversion.
                    ds['tdps'].attrs['long_name'] = "dew_point_temperature"
                    ds['tdps'].attrs['standard_name'] = "dew_point_temperature"
                    ds['tdps'].attrs['units'] = "degree_Kelvin"

                    if 'DPTP_flag' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'DPTP_flag': 'tdps_qc'})
                        ds['tdps_qc'].attrs['flag_values'] = "V S E"
                        ds['tdps_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['tdps'].attrs['ancillary_variables'] = "tdps_qc" # List other variables associated with variable (QA/QC)

                    ds['tdps'].attrs['comment'] = "Converted from Fahrenheit to Kelvin."


                # pr: precipitation
                # We have  different raw precipitation variables for precip.
                # 'PREC', # Precipitation accumulation (in)
                # 'PRCP', # Precipitation increment (in)
                # 'PRCPSA', # Precipitation increment snow-adjusted (in)

                # At this stage, no infilling. So we will keep all columns and simply rename them.

                # Precipitation accumulation
                if "PREC_value" in ds.keys():
                    ds['pr'] = calc_clean._unit_precip_in_to_mm(ds['PREC_value'])
                    ds = ds.drop("PREC_value")
                    ds['pr'].attrs['long_name'] = "precipitation_accumulation"
                    ds['pr'].attrs['units'] = "mm/?"

                    if 'PREC_flag' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'PREC_flag': 'pr_qc'})
                        ds['pr_qc'].attrs['flag_values'] = "V S E"
                        ds['pr_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['pr'].attrs['ancillary_variables'] = "pr_qc" # List other variables associated with variable (QA/QC)

                    ds['pr'].attrs['comment'] = "Accumulated precipitation. Converted from inches to mm."

                # Precipitation increment
                if 'PRCP_value' in ds.keys():
                    ds['pr_inc'] = calc_clean._unit_precip_in_to_mm(ds['PRCP_value'])
                    ds = ds.drop("PRCP_value")
                    ds['pr_inc'].attrs['long_name'] = "precipitation_increment"
                    ds['pr_inc'].attrs['units'] = "mm/?"

                    if 'PRCP_flag' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'PRCP_flag': 'pr_inc_qc'})
                        ds['pr_inc_qc'].attrs['flag_values'] = "V S E"
                        ds['pr_inc_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['pr_inc'].attrs['ancillary_variables'] = "pr_inc_qc" # List other variables associated with variable (QA/QC)

                    ds['pr_inc'].attrs['comment'] = "Precipitation increment. Converted from inches to mm."


                if 'PRCPSA_value' in ds.keys():
                    ds['pr_incsa'] = calc_clean._unit_precip_in_to_mm(ds['PRCPSA_value'])
                    ds = ds.drop("PRCPSA_value")
                    ds['pr_incsa'].attrs['long_name'] = "precipitation_increment_snow_adjusted"
                    ds['pr_incsa'].attrs['units'] = "mm/?"

                    if 'PRCPSA_flag' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'PRCPSA_flag': 'pr_incsa_qc'})
                        ds['pr_incsa_qc'].attrs['flag_values'] = "V S E"
                        ds['pr_incsa_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['pr_incsa'].attrs['ancillary_variables'] = "pr_incsa_qc" # List other variables associated with variable (QA/QC)

                    ds['pr_incsa'].attrs['comment'] = "Precipitation increment (snow-adjusted). Converted from inches to mm."


                # hurs: relative humidity (%)
                if 'RHUM_value' in ds.keys(): # Already in %, no need to convert units.
                    ds = ds.rename({'RHUM_value': 'hurs'})
                    # Set attributes
                    ds['hurs'].attrs['long_name'] = "relative_humidity"
                    ds['hurs'].attrs['standard_name'] = "relative_humidity"
                    ds['hurs'].attrs['units'] = "percent"

                    # If QA/QC column exists
                    if 'RHUM_flag' in ds.keys():
                        ds = ds.rename({'RHUM_flag': 'hurs_qc'})
                        ds['hurs_qc'].attrs['flag_values'] = "V S E"
                        ds['hurs_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['hurs'].attrs['ancillary_variables'] = "hurs_qc" # List other variables associated with variable (QA/QC)


                #rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
                if 'SRAD_value' in ds.keys(): # Already in w/m2, no need to convert units.
                    # If column exists, rename.
                    ds = ds.rename({'SRAD_value': 'rsds'})

                    # Set attributes
                    ds['rsds'].attrs['long_name'] = "solar_radiation"
                    ds['rsds'].attrs['standard_name'] = "surface_downwelling_shortwave_flux_in_air"
                    ds['rsds'].attrs['units'] = "W m-2"

                    #rsds: QA/QC flags
                    if 'SRAD_flag' in ds.keys():
                        ds = ds.rename({'SRAD_flag': 'rsds_qc'})
                        ds['rsds_qc'].attrs['flag_values'] = "V S E"
                        ds['rsds_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['rsds'].attrs['ancillary_variables'] = "rsds_qc" # List other variables associated with variable (QA/QC)


                # sfcWind : wind speed (m/s)
                if "WSPD_value" in ds.keys(): # Data originally in mph.
                    ds['sfcWind'] = calc_clean._unit_windspd_mph_to_ms(ds['WSPD_value'])
                    ds = ds.drop("WSPD_value")
                    ds['sfcWind'].attrs['long_name'] = "wind_speed"
                    ds['sfcWind'].attrs['standard_name'] = "wind_speed"
                    ds['sfcWind'].attrs['units'] = "m s-1"

                    #rsds: QA/QC flags
                    if 'WSPD_flag' in ds.keys():
                        ds = ds.rename({'WSPD_flag': 'sfcWind_qc'})
                        ds['sfcWind_qc'].attrs['flag_values'] = "V S E"
                        ds['sfcWind_qc'].attrs['flag_meanings'] =  "valid suspect edited"
                        ds['sfcWind'].attrs['ancillary_variables'] = "sfcWind_qc" # List other variables associated with variable (QA/QC)

                    ds['sfcWind'].attrs['comment'] = "Converted from mph to m/s."


                # sfcWind_dir: wind direction
                if "WDIR_value" in ds.keys(): # No conversions needed, do not make raw column.
                    ds = ds.rename({'WDIR_value': 'sfcWind_dir'})
                    ds['sfcWind_dir'].attrs['long_name'] = "wind_direction"
                    ds['sfcWind_dir'].attrs['standard_name'] = "wind_from_direction"
                    ds['sfcWind_dir'].attrs['units'] = "degrees_clockwise_from_north"

                    if 'WDIR_flag' in ds.keys():
                        ds = ds.rename({'WDIR_flag': 'sfcWind_dir_qc'})
                        ds['sfcWind_dir_qc'].attrs['flag_values'] = "V S E"
                        ds['sfcWind_dir_qc'].attrs['flag_meanings'] = "valid suspect edited"
                        ds['sfcWind_dir'].attrs['ancillary_variables'] = "sfcWind_dir_qc" # List other variables associated with variable (QA/QC)

                    ds['sfcWind_dir'].attrs['comment'] = "Wind direction is defined by the direction that the wind is coming from (i.e., a northerly wind originates in the north and blows towards the south)."


                # Other variables: rename to match format

                # Partial vapor pressure (kPa -> Pa) ## No CMIP standard name for this var, just CF.
                if 'PVPV_value' in ds.keys():
                    ds['pvp'] = calc_clean._unit_pres_kpa_to_pa(ds['PVPV_value'])
                    ds = ds.drop("PVPV_value")

                    ds['pvp'].attrs['long_name'] = "partial_vapor_pressure"
                    ds['pvp'].attrs['standard_name'] = "water_vapor_partial_pressure_in_air"
                    ds['pvp'].attrs['units'] = "Pa"

                    if "PVPV_flag" in ds.keys():
                        ds = ds.rename({'PVPV_flag': 'pvp_qc'})
                        ds['pvp_qc'].attrs['flag_values'] = "V S E"
                        ds['pvp_qc'].attrs['flag_meanings'] = "valid suspect edited"
                        ds['pvp'].attrs['ancillary_variables'] = "pvp_qc" # List other variables associated with variable (QA/QC)

                    ds['pvp'].attrs['comment'] = "Converted from kPa to Pa."

                # Saturated vapor pressure (kPa -> Pa)
                if 'SVPV_value' in ds.keys():
                    ds['svp'] = calc_clean._unit_pres_kpa_to_pa(ds['SVPV_value'])
                    ds = ds.drop("SVPV_value")

                    ds['svp'].attrs['long_name'] = "saturated_vapor_pressure"
                    ds['svp'].attrs['units'] = "Pa"

                    if "SVPV_flag" in ds.keys():
                        ds = ds.rename({'SVPV_flag': 'svp_qc'})
                        ds['svp_qc'].attrs['flag_values'] = "V S E"
                        ds['svp_qc'].attrs['flag_meanings'] = "valid suspect edited"
                        ds['svp'].attrs['ancillary_variables'] = "svp_qc" # List other variables associated with variable (QA/QC)

                    ds['svp'].attrs['comment'] = "Converted from kPa to Pa."

                #Testing: Manually check values to see that they seem correctly scaled, no unexpected NAs.
                #for var in ds.variables:
                    # try:
                    #     print([var, float(ds[var].min()), float(ds[var].max())])
                    # except:
                    #     next

                # Quality control: if any variable is completely empty, drop it.
                # drop any column that does not have any valid (non-nan) data
                # need to keep elevation separate, as it does have "valid" nan value, only drop if all other variables are also nans
                for key in ds.keys():
                    try:
                        if (key != 'elevation'):
                            if np.isnan(ds[key].values).all():
                                print("Dropping empty var: {}".format(key))
                                ds = ds.drop(key)

                        # only drop elevation if all other variables are also nans
                        if (key == 'elevation') & (len(ds.keys())==1): # only elevation remains
                            print("Dropping empty var: {}".format(key)) # slightly unnecessary since the entire dataset will be empty too
                            ds = ds.drop(key)
                            continue
                    except Exception as e: # Add to handle errors for unsupported data types
                        next

                # For QA/QC flags, replace np.nan with "nan" to avoid h5netcdf overwrite to blank.
                for key in ds.keys():
                    if 'qc' in key:
                        ds[key] = ds[key].astype(str) # Coerce all values in key to string.

                # Reorder variables
                desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir', 'pvp', 'svp']
                actual_order = [i for i in desired_order if i in list(ds.keys())]
                rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom.
                new_index = actual_order + rest_of_vars
                ds = ds[new_index]

            except Exception as e:
                # print(e)
                errors['File'].append(station_id) # If stat_files is none, this will default to saving ID of station.
                errors['Time'].append(end_api)
                errors['Error'].append('Error in ds set-up: {}'.format(e))
                continue

            #Write station file to netcdf.
            if len(ds.keys())==0:   # skip station if the entire dataset will be empty because no data is observed (as in only ocean obs are recorded, but not needed)
                print("{} has no data for all meteorological variables of interest throughout its current reporting; station not cleaned.".format(station_id))
                errors['File'].append(station_id)
                errors['Time'].append(end_api)
                errors['Error'].append('Station reports all nan meteorological data.')
                continue

            elif len(ds.time) == 0: # this should not be necessary, but placing for testing
                print('No data for this station during v1 period (1/1980 - 8/2022), station not cleaned.')
                errors['File'].append(stat_files)
                errors['Time'].append(end_api)
                errors['Error'].append("No data for this station during v1 period (1/1980 - 8/2022") 
                continue

            else:
                try:
                    filename = station_id+".nc" # Make file name
                    filepath = cleandir+filename # Write file path
                    #print(filepath) # For testing

                    # Write locally
                    ds.to_netcdf(path = 'temp/temp.nc', engine = 'netcdf4') # Save station file.

                    # Push file to AWS with correct file name.
                    s3.Bucket(bucket_name).upload_file('temp/temp.nc', filepath)

                    print("Saving {} with dims {}".format(filename, ds.dims))
                    ds.close() # Close dataframe.
                    continue

                except Exception as e:
                    # print(e)
                    errors['File'].append(stat_files)
                    errors['Time'].append(end_api)
                    errors['Error'].append('Error saving ds as .nc file to AWS bucket: {}'.format(e))
                    continue


    # The list of removed variables is manually saved to AWS, as it includes many variables not initially downloaded.
    # Obtained from: https://www.nrcs.usda.gov/wps/portal/wcc/home/dataAccessHelp/webService/webServiceReference

    # Write errors.csv
    finally: # Always execute this.
        # print(errors)
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"errors_{}_{}.csv".format(network, end_api))


# # Run functions
if __name__ == "__main__":
    network = "SCAN"
    rawdir, cleandir, qaqcdir = get_file_paths(network)
    print(rawdir, cleandir, qaqcdir)
    clean_scansnotel(rawdir, cleandir)

#----------------------------------------------------------------------------------------
    ## Testing:
    # import random # To get random subsample
    # import s3fs # To read in .nc files

    # # # ## Import file.
    # files = []
    # for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
    #     file = str(item.key)
    #     files += [file]

    # files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names
    # files = [file for file in files if "error" not in file] # Remove error handling files.
    # files = [file for file in files if "station" not in file] # Remove error handling files.
    # files = random.sample(files, 4)

    # # File 1:
    # fs = s3fs.S3FileSystem()
    # aws_urls = ["s3://wecc-historical-wx/"+file for file in files]

    # with fs.open(aws_urls[0]) as fileObj:
    #     test = xr.open_dataset(fileObj)
    #     print(test)
    #     for var in test.keys():
    #         print(var)
    #         print(test[var])
    #     test.close()

    # # File 2:
    # # Test: multi-year merges work as expected.
    # with fs.open(aws_urls[1]) as fileObj:
    #     test = xr.open_dataset(fileObj, engine='h5netcdf')
    #     print(str(test['time'].min())) # Get start time
    #     print(str(test['time'].max())) # Get end time
    #     test.close()


    # # File 3:
    # # Test: Inspect vars and attributes
    # with fs.open(aws_urls[2]) as fileObj:
    #     test = xr.open_dataset(fileObj, engine='h5netcdf')
    #     for var in test.variables:
    #         try:
    #             print([var, float(test[var].min()), float(test[var].max())])
    #         except:
    #             continue
    #     test.close()
