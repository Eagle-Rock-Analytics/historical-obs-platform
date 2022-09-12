"""
This script performs data cleaning for the CWOP network for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables, and calculates derived priority variables if not observed
(2) Drops unnecessary variables
(3) Converts station metadata to standard format, with unique identifier
(4) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(5) Converts missing data to standard format
(6) Tracks existing qa/qc flag for review
(7) Merge files by station, and outputs cleaned variables as a single .nc file for an individual network.
Inputs: Raw data for CWOP stations, with each csv file representing a station.
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
Reference: https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf
"""
# To do:
#  fix dividing by 0 error



# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime
import re
import numpy as np
import csv
import gzip
import pandas as pd
import math
import requests
from collections import Counter

# Import calc_clean.py.
try:
    import calc_clean
except:
    print("Error importing calc_clean.py")
    pass

# Import config.py
try:
    import config # Import API keys.
except:
    print("Missing config.py file with API token. Make file if necessary.")
    exit()

# Set envr variables and calculate any needed variables

# Set path to head of git repository.
homedir = os.getcwd() # Get current working directory.
if "historical-obs-platform" in homedir: # If git folder in path
    homedir = homedir[0:homedir.index("historical-obs-platform")]+"historical-obs-platform" # Set path to top folder.
    os.chdir(homedir) # Change directory.
else:
    print("Error: Set current working directory to the git repository or a subfolder, and then rerun script.")
    exit()

# Set relative paths to other folders and objects in repository.
workdir = "test_platform/data/1_raw_wx/CWOP/"
savedir = "test_platform/data/2_clean_wx/CWOP/"
wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 

# Set up directory to save files, if it doesn't already exist.
try:
    os.mkdir(savedir) # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass

## FUNCTION: Get CWOP QA/QC flag data and parse into csv.
# Input: API url for QAQC metadata.
# Output: parsed table of QAQC codes and flag values, saved as csv.
# Note: this function does not need to be run daily, just occasionally.
def get_qaqc_flags(token, savedir):
    url = "https://api.synopticdata.com/v2/qctypes?token={}".format(token)
    print(url)
    request = requests.get(url).json()
    ids = []
    for each in request['QCTYPES']:
        ids.append([each['ID'], each['SHORTNAME'], each['NAME']])
    ids = pd.DataFrame(ids, columns = ['ID', 'SHORTNAME', 'NAME']) # Sort by start date (note some stations return 'None' here)
    
    # Write to csv.
    filepath = savedir+'qaqc_flags.csv' # Set path to desired folder. # Write file to name of station ID in synoptic-- Change file name to reflect dates?? 
    print(ids)
    print("Saving QAQC flags to csv file.")
    ids.to_csv(filepath, index = False)

#get_qaqc_flags(token = config.token, savedir = savedir)

## FUNCTION: Clean CWOP data.
# Input: 
# homedir: path to git repository.
# workdir: path to where raw data is saved as .csv files, with each file representing a station's records from download start date to present (by default 01-01-1980).
# Some stations have more than one csv file, with the second file named [STID]_2.csv.
# savedir: path to where cleaned files should be saved.
# Output: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
# Note: files are already filtered by bbox in the pull step, so additional geospatial filtering is not required here.

def clean_cwop(homedir, workdir, savedir):

    os.chdir(homedir)
    os.chdir(workdir)
    
    files = os.listdir() # Gets list of files in directory to work with
    files = list(filter(lambda f: f.endswith(".csv"), files)) # Get list of file names
    files = list(filter(lambda f: not f.startswith("error"), files)) # Remove error handling files.
    
    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
    timestamp = datetime.now().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.

    # Set up list of variables to be removed
    removedvars = []

    # # Get list of station IDs from filename and clean.
    ids = list()
    for file in files:
        if 'error' not in file: # Ignore all errors csv
            # Remove errors_[....].csv file names - TO DO.
            id = re.sub('_2', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
            id = re.sub('.csv', "", id) # Remove leading NDBC
            if id not in ids:
                ids.append(id)
    
    #for i in ids: # For each station
    for i in ids[0:5]: # Subsample for testing
    #for i in ['AP795']: # For testing merge (on AWS)    
        ds_stat = None # Initialize ds_stat
        file_count = 1 # Initialize file counter
        stat_files = [k for k in files if i in k] # Get list of files with station ID in them.
        station_id = "CWOP_"+i
        print(stat_files)
        
        # Iterate through files to clean. Each file represents a station's data.
        for file in stat_files:
            try:
                with open(file, "r") as csv_file:
                    csv_reader = csv.reader(csv_file)
                    for index, row in enumerate(csv_reader):
                        
                        # Parse header information from CSV.
                        if len(row)==1:
                            row = str(row)
                            
                            if "STATION NAME" in row:
                                station_name = row.partition(": ")[2].replace("']", "")
                                station_name = station_name.replace(")", "")
                            elif "LATITUDE" in row:
                                latitude = float(row.partition(": ")[2].replace("']", ""))
                            elif "LONGITUDE" in row:
                                longitude = float(row.partition(": ")[2].replace("']", ""))
                            elif "ELEVATION" in row:
                                elevation = float(row.partition(": ")[2].replace("']", "")) # In feet.
                            elif "STATE" in row:
                                state = row.partition(": ")[2].replace("']", "")
                            continue

                        else:
                            # Get column names
                            if "Station_ID" in row:
                                columns = row
                                continue

                            # Get units
                            unitstocheck = ['Pascals', '%', 'm/s', 'Celsius', 'QC_Type', 'Degrees']
                            if any(unit in row for unit in unitstocheck):
                                units_index = index
                                continue
                                
                            # Get row where data begins
                            if file.replace(".csv", "") in row:
                                first_row = index
                                break
                        
                # Read in entire csv starting at first_row, using pandas.
                
                # First, some error handling.
                
                    #Check for duplicated columns. If duplicated names, manually rename as 1 and 2, and then check for duplicates after reading in.
                    
                    if set([x for x in columns if columns.count(x) > 1]):
                        dup = True # Add dup flag.
                        dup_col = set([x for x in columns if columns.count(x) > 1])
                        # Manually rename second iteration of column.
                        d = {a:list(range(1, b+1)) if b>1 else '' for a,b in Counter(columns).items()}
                        columns = [i+str(d[i].pop(0)) if len(d[i]) else i for i in columns]
                    else:
                        dup = False
                
                
                df = pd.read_csv(file, names = columns, header = first_row-1, low_memory = False)

                # If columns duplicated, check to see if they are identical and drop second if so.
                if dup is True:
                    for i in dup_col: # For each duplicated column
                        cols = df.filter(like=i).columns
                        if df[cols[0]].equals(df[cols[1]]):
                            #print("Dropping duplicate column {}".format(cols[1]))
                            df.drop(cols[1], axis = 1, inplace = True)
                            #print(df.columns)
                        else:
                            print("None-identical duplicate columns found.") # For testing
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("Non-identical duplicate columns found. Columns: {}".format(i))
                            continue

                # Fix time format issues caused by "Timeout" errors.
                df['Date_Time'] = pd.to_datetime(df['Date_Time'], errors = 'coerce') # Convert time to datetime and catch any incorrectly formatted columns.
                df = df[pd.notnull(df['Date_Time'])] # Remove any rows where time is missing.
                
                # Remove any non-essential columns.
                coltokeep = ['Station_ID', 'Date_Time', 'altimeter_set_1', 'altimeter_set_1_qc',
                            'air_temp_set_1', 'air_temp_set_1_qc', 'relative_humidity_set_1',
                            'relative_humidity_set_1_qc', 'wind_speed_set_1', 'wind_speed_set_1_qc',
                            'wind_direction_set_1', 'wind_direction_set_1_qc', 'precip_accum_since_local_midnight_set_1',
                            'precip_accum_since_local_midnight_set_1_qc',
                            'precip_accum_24_hour_set_1', 'precip_accum_24_hour_set_1_qc',
                            'dew_point_temperature_set_1d', 'pressure_set_1d',
                            'pressure_set_1', 'pressure_set_1_qc', 'dew_point_temperature_set_1', 'dew_point_temperature_set_1_qc', 
                            'precip_accum_one_hour_set_1', 'precip_accum_one_hour_set_1_qc', 'solar_radiation_set_1', 'solar_radiation_set_1_qc', 
                            'precip_accum_set_1', 'precip_accum_set_1_qc', 'precip_accum_five_minute_set_1', 'precip_accum_five_minute_set_1_qc']
                
                othercols = [col for col in df.columns if col not in coltokeep and col not in removedvars] 
                removedvars += othercols # Add any new columns from drop list to removedvars, to save later.
                df = df.drop(columns=[col for col in df if col not in coltokeep]) # Drop all columns not in coltokeep list.
                
                # # Manually convert "None" to np.nan
                df.replace(to_replace="None", value=np.nan, inplace=True)

                # Fix multi-type columns
                # If column has QC in it, force to string.
                for i in df.columns:
                    multitype = set(type(x).__name__ for x in df[i])
                    if len(multitype)>1:
                        if 'qc' in i:
                            df[i] = df[i].astype(str) # Coerce to string (to handle multiple QA/QC flags)
                        elif 'wind_cardinal_direction' in i:
                            df[i] = df[i].astype(str) # Coerce to string
                        elif 'sea_level_pressure_set' in i:
                            df[i] = df[i].astype(float) # Coerce to float.
                        else:
                            print("Multitype error for column {} with data types {}. Please resolve".format(i, multitype)) # Manually 
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("Multitype error for column {} with data types {}. Please resolve".format(i, multitype))
                            continue
                                
                # Move df to xarray object.
                ds = df.to_xarray()
                
                # Update dimensions and coordinates

                # Add dimensions: station ID and time.
                ds = ds.rename({'Date_Time': 'time'}) # Rename time variable.
                ds = ds.set_coords('time').swap_dims({'index': 'time'}) # Swap index with time.
                ds = ds.assign_coords(id = str(station_id))
                ds = ds.expand_dims("id") # Add station_id as index.
                ds = ds.drop_vars(("index")) # Drop station_id variable and index coordinate.
                ds = ds.rename({'id': 'station'}) # Rename id to station_id.
                
                    
                # Add coordinates: latitude and longitude.
                lat = np.asarray([latitude]*len(ds['time']))
                lat.shape = (1, len(ds['time']))
                
                lon = np.asarray([longitude]*len(ds['time']))
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
                
                # Add variable: elevation
                elev = np.asarray([elevation]*len(ds['time']))
                elev.shape = (1, len(ds['time']))
                ds['elevation'] = (['station', 'time'], elev)
                
                # Update dimension and coordinate attributes.
                # Time
                
                # Convert column to datetime (and remove any rows that cannot be coerced).
                ds['time'] = pd.to_datetime(ds['time'].values, utc = True) # Convert from string to incorrect time format.
                ds['time'] = pd.to_datetime(ds['time'].values, unit = 'ns') # Fix time format.
                
                # Update attributes.
                ds['time'].attrs['long_name'] = "time"
                ds['time'].attrs['standard_name'] = "time"
                ds['time'].attrs['comment'] = "In UTC."
                
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
                ds['elevation_raw'] = ds['elevation']
                ds['elevation'] = calc_clean._unit_elev_ft_to_m(ds['elevation']) # Convert feet to meters.

                ds['elevation'].attrs['standard_name'] = "height_above_mean_sea_level"
                ds['elevation'].attrs['long_name'] = "station_elevation"
                ds['elevation'].attrs['units'] = "meters"
                ds['elevation'].attrs['positive'] = "up" # Define which direction is positive
                ds['elevation'].attrs['comment'] = "Converted from feet"
                
                ds['elevation_raw'].attrs['standard_name'] = "height_above_mean_sea_level"
                ds['elevation_raw'].attrs['long_name'] = "station_elevation"
                ds['elevation_raw'].attrs['units'] = "feet"
                ds['elevation_raw'].attrs['positive'] = "up" # Define which direction is positive
                
                # Update global attributes
                ds = ds.assign_attrs(title = "CWOP cleaned")
                ds = ds.assign_attrs(institution = "Eagle Rock Analytics / Cal Adapt")
                ds = ds.assign_attrs(source = "")
                ds = ds.assign_attrs(history = "CWOP_clean.py script run on {} UTC".format(timestamp))
                ds = ds.assign_attrs(comment = "Intermediate data product: may not have been subject to any cleaning or QA/QC processing")
                ds = ds.assign_attrs(license = "")
                ds = ds.assign_attrs(citation = "")
                ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")
                ds = ds.assign_attrs(station_name = station_name)

                # Update variable attributes and do unit conversions
                
                #tas: air surface temperature (K)
                if "air_temp_set_1" in ds.keys():
                    ds = ds.rename({'air_temp_set_1': 'tas_raw'})
                    ds['tas'] = calc_clean._unit_degC_to_K(ds['tas_raw'])
                    
                    ds['tas'].attrs['long_name'] = "air_temperature"
                    ds['tas'].attrs['standard_name'] = "air_temperature"
                    ds['tas'].attrs['units'] = "degree_Kelvin"
                    ds['tas'].attrs['ancillary_variables'] = "tas_raw" # List other variables associated with variable (QA/QC)
                    ds['tas'].attrs['comment'] = "Converted from Celsius."
                    
                    ds['tas_raw'].attrs['long_name'] = "air_temperature"
                    ds['tas_raw'].attrs['standard_name'] = "air_temperature"
                    ds['tas_raw'].attrs['units'] = "degree_Celsius"
                    ds['tas_raw'].attrs['ancillary_variables'] = "tas" # List other variables associated with variable (QA/QC)
                
                    if "air_temp_set_1_qc" in ds.keys():
                        # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
                        
                        flagvals = list(np.unique(ds['air_temp_set_1_qc'].values))
                        flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                        
                        if flagvals == ['nan']:
                            #print("Flag value is {}".format(flagvals)) # For testing.
                            ds = ds.drop("air_temp_set_1_qc")
                        else:
                            ds = ds.rename({'air_temp_set_1_qc': 'tas_qc'})
                            flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                            ds['tas_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                            ds['tas_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                            
                            ds['tas'].attrs['ancillary_variables'] = "tas_raw tas_qc" # List other variables associated with variable (QA/QC)
                            ds['tas_raw'].attrs['ancillary_variables'] = "tas tas_qc"
                
                else: # If station doesn't have air temperature, add nan variable with same attributes.
                    tas = np.full((1, len(ds['time'])), np.nan)
                    ds['tas'] = (['station', 'time'], tas)
                
                    ds['tas'].attrs['long_name'] = "air_temperature"
                    ds['tas'].attrs['standard_name'] = "air_temperature"
                    ds['tas'].attrs['units'] = "degree_Kelvin"

                # ps: surface air pressure (Pa)
                # Note here that if "pressure_set_1" has values this is a direct station observation reading.
                # Otherwise, if "pressure_set_1d" has values this is a derived value calculated from altimeter and elevation.
                # We will manually recalculate this here.

                # Set up column.
                ps = np.full((1, len(ds['time'])), np.nan)
                ds['ps'] = (['station', 'time'], ps)
                ds['ps'].attrs['long_name'] = "station_air_pressure"
                ds['ps'].attrs['standard_name'] = "air_pressure"
                ds['ps'].attrs['units'] = "Pa"

                if 'pressure_set_1' in ds.keys():
                    if not np.isnan(ds['pressure_set_1'].values).all(): # If station pressure directly observed
                        ds['ps'].values = ds['pressure_set_1'].values
                        ds = ds.drop("pressure_set_1")
                        #print("Renamed column")

                    elif np.isnan(ds['pressure_set_1'].values).all(): # If all direct station observations are NA
                        if 'altimeter_set_1' in ds.keys(): 
                            if not np.isnan(ds['altimeter_set_1'].values).all(): # If altimeter readings observed.
                                # Calculate pressure from altimeter and elevation
                                ds['ps'].values = calc_clean._calc_ps_alt(ds['altimeter_set_1'], ds['elevation']) 
                                ds['ps'].attrs['comment'] = "Calculated from altimeter setting, station elevation using calc_clean.py."
                                ds['ps'].attrs['ancillary_variables'] = "altimeter_set_1 elevation" # List other variables associated with variable (QA/QC)
                                #print("Calculated ps from altimeter.")

                            else: # If altimeter readings all NA
                                if not np.isnan(ds['sea_level_pressure_set_1'].values).all():
                                    ds['ps'].values = calc_clean._calc_ps(ds['sea_level_pressure_set_1'], ds['elevation'], ds['tas']) 
                                    ds['ps'].attrs['comment'] = "Calculated from sea level pressure, station elevation, air temperature using calc_clean.py."
                                    ds = ds.rename({'sea_level_pressure_set_1_qc': 'psl'})
                                    ds['ps'].attrs['ancillary_variables'] = "psl elevation tas" # List other variables associated with variable (QA/QC)
                                
                        else: # If no altimeter readings
                            if not np.isnan(ds['sea_level_pressure_set_1'].values).all():
                                ds['ps'].values = calc_clean._calc_ps(ds['sea_level_pressure_set_1'], ds['elevation'], ds['tas']) 
                                ds['ps'].attrs['comment'] = "Calculated from sea level pressure, station elevation, air temperature using calc_clean.py."
                                ds = ds.rename({'sea_level_pressure_set_1_qc': 'psl'})
                                ds['ps'].attrs['ancillary_variables'] = "psl elevation tas" # List other variables associated with variable (QA/QC)
                            
                else: # If no pressure column exists.
                    if 'altimeter_set_1' in ds.keys(): 
                        if not np.isnan(ds['altimeter_set_1'].values).all(): # If altimeter readings observed.
                                
                                # Calculate pressure from altimeter and elevation
                                ds['ps'].values = calc_clean._calc_ps_alt(ds['altimeter_set_1'], ds['elevation']) 
                                ds['ps'].attrs['comment'] = "Calculated from altimeter setting, station elevation using calc_clean.py."
                                ds['ps'].attrs['ancillary_variables'] = "altimeter_set_1 elevation" # List other variables associated with variable (QA/QC)

                        elif np.isnan(ds['altimeter_set_1'].values).all(): # If altimeter readings not observed:
                            if not np.isnan(ds['sea_level_pressure_set_1'].values).all():
                                ds['ps'].values = calc_clean._calc_ps(ds['sea_level_pressure_set_1'], ds['elevation'], ds['tas']) 
                                ds['ps'].attrs['comment'] = "Calculated from sea level pressure, station elevation, air temperature using calc_clean.py."
                                ds = ds.rename({'sea_level_pressure_set_1_qc': 'psl'})
                                ds['ps'].attrs['ancillary_variables'] = "psl elevation tas" # List other variables associated with variable (QA/QC)
                        
                if 'pressure_set_1_qc' in ds.keys(): # If QA/QC exists.
                    flagvals = list(np.unique(ds['pressure_set_1_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pressure_set_1_qc")
                    else:
                        ds = ds.rename({'pressure_set_1_qc': 'ps_qc'})
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                        ds['ps_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                        ds['ps_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        
                        if 'ancillary_variables' in ds['ps'].attrs:
                            ds['ps'].attrs['ancillary_variables'] = ds['ps'].attrs['ancillary_variables'] + " ps_qc" # List other variables associated with variable (QA/QC)
                        else:
                            ds['ps'].attrs['ancillary_variables'] = "ps_qc"

                            
                # tdps: dew point temperature (K)
                # if raw dew point temperature observed, use that.
                if 'dew_point_temperature_set_1' in ds.keys():
                    ds = ds.rename({'dew_point_temperature_set_1': 'tdps_raw'})
                    ds['tdps'] = calc_clean._unit_degC_to_K(ds['tdps_raw'])
                    
                    # Set attributes for conversion.
                    ds['tdps'].attrs['ancillary_variables'] = "tdps_raw" # List other variables associated with variable (QA/QC)
                    ds['tdps'].attrs['comment'] = "Converted from Celsius."
                
                    # Set attributes for raw column.
                    ds['tdps_raw'].attrs['long_name'] = "dew_point_temperature"
                    ds['tdps_raw'].attrs['standard_name'] = "dew_point_temperature"
                    ds['tdps_raw'].attrs['units'] = "degree_Celsius"
                    ds['tdps_raw'].attrs['ancillary_variables'] = "tdps" # List other variables associated with variable (QA/QC)

                    # QAQC flag
                    if "dew_point_temperature_set_1_qc" in ds.keys():
                    # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
                    
                        flagvals = list(np.unique(ds['dew_point_temperature_set_1_qc'].values))
                        flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                        
                        if flagvals == ['nan']:
                            #print("Flag value is {}".format(flagvals)) # For testing.
                            ds = ds.drop("dew_point_temperature_set_1_qc")
                        else:
                            ds = ds.rename({'dew_point_temperature_set_1_qc': 'tdps_qc'})
                            flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                            ds['tdps_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                            ds['tdps_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                            
                            ds['tdps'].attrs['ancillary_variables'] = "tdps_raw tdps_qc" # List other variables associated with variable (QA/QC)
                            ds['tdps_raw'].attrs['ancillary_variables'] = "tdps tdps_qc"
                

                # If vapor pressure doesn't exist but tas and hurs do, calculate using opt 1.
                elif 'tas' in ds.keys() and 'relative_humidity_set_1' in ds.keys():
                    # Note we don't check for non-NA values here because there's no alternative calculation method.
                    # If calculation returns all NAs, this is equivalent to manually setting the column to be NA.
                    # Inputs: tas (K) and relative humidity (%)
                    ds['tdps'] = calc_clean._calc_dewpointtemp_opt1(ds['tas'], ds['relative_humidity_set_1'])

                    # Set attributes for calculation.
                    ds['tdps'].attrs['ancillary_variables'] = "tas hurs" # List other variables associated with variable (QA/QC)
                    ds['tdps'].attrs['comment'] = "Calculated from air temperature and relative humidity using calc_clean.py."
                
                # Vapor pressure not collected by CWOP, so can't calculate tdps using opt2.
                # If none of these variables available, make NA column.
                else:
                    tdps = np.full((1, len(ds['time'])), np.nan)
                    ds['tdps'] = (['station', 'time'], tdps)
                
                # Set attributes
                ds['tdps'].attrs['long_name'] = "dew_point_temperature"
                ds['tdps'].attrs['standard_name'] = "dew_point_temperature"
                ds['tdps'].attrs['units'] = "degree_Kelvin"
                
                
                # pr: precipitation
                # We have 4 different raw precipitation variables for precip.
                # precip_accum_24_hour_set_1 # Precipitation from last 24 hours.
                # precip_accum_since_local_midnight_set_1 # Precipitation from local midnight.
                # precip_accum_set_1 # Precipitation since last record.
                # precip_accum_one_hour_set_1 # Precipitation in last hour.

                # At this stage, no infilling. So we will keep all columns and simply rename them.
                # Remove any column if completely empty.
                precip_cols = [elem for elem in ds.keys() if "precip" in elem]
                for col in precip_cols:
                    if np.isnan(ds[col].values).all():
                        #print("Dropping {}".format(col))
                        ds = ds.drop(col)

                # Reformat remaining columns
                if "precip_accum_24_hour_set_1" in ds.keys():
                    ds = ds.rename({'precip_accum_24_hour_set_1': 'pr_24h'})
                    ds['pr_24h'].attrs['long_name'] = "24_hr_precipitation"
                    ds['pr_24h'].attrs['units'] = "mm"
                    ds['pr_24h'].attrs['comment'] = "Precipitation accumulated in previous 24 hour period."

                    if 'precip_accum_24_hour_set_1_qc' in ds.keys():
                        ds = ds.rename({'precip_accum_24_hour_set_1_qc': 'pr_24h_qc'})
                
                if 'precip_accum_since_local_midnight_set_1' in ds.keys():
                    ds = ds.rename({'precip_accum_since_local_midnight_set_1': 'pr_localmid'})
                    ds['pr_localmid'].attrs['long_name'] = "precipitation_since_local_midnight"
                    ds['pr_localmid'].attrs['units'] = "mm"
                    ds['pr_localmid'].attrs['comment'] = "Precipitation accumulated since local midnight."

                    if 'precip_accum_since_local_midnight_set_1_qc' in ds.keys():
                        ds = ds.rename({'precip_accum_since_local_midnight_set_1_qc': 'pr_localmid_qc'})
                
                if "precip_accum_set_1" in ds.keys():
                    ds = ds.rename({'precip_accum_set_1': 'pr_raw'})
                    ds['pr_raw'].attrs['long_name'] = "precipitation"
                    ds['pr_raw'].attrs['units'] = "mm"
                    ds['pr_raw'].attrs['comment'] = "Precipitation accumulated since previous measurement."

                    if 'precip_accum_set_1_qc' in ds.keys():
                        ds = ds.rename({'precip_accum_set_1_qc': 'pr_raw_qc'})

                if "precip_accum_one_hour_set_1" in ds.keys():
                    ds = ds.rename({'precip_accum_one_hour_set_1': 'pr_1h'})
                    ds['pr_1h'].attrs['long_name'] = "hourly_precipitation"
                    ds['pr_1h'].attrs['units'] = "mm"
                    ds['pr_1h'].attrs['comment'] = "Precipitation accumulated in previous hour."

                    if 'precip_accum_one_hour_set_1_qc' in ds.keys():
                        ds = ds.rename({'precip_accum_one_hour_set_1_qc': 'pr_1h_qc'})

                # Reformat qc columns
                if "pr_24h_qc" in ds.keys():
                    flagvals = list(np.unique(ds['pr_24h_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pr_24h_qc")
                    else:
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                        ds['pr_24h_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                        ds['pr_24h_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        
                        ds['pr_24h'].attrs['ancillary_variables'] = "pr_24h_qc" # List other variables associated with variable (QA/QC)

                if "pr_localmid_qc" in ds.keys():
                    flagvals = list(np.unique(ds['pr_localmid_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pr_localmid_qc")
                    else:
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                        ds['pr_localmid_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                        ds['pr_localmid_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        
                        ds['pr_localmid'].attrs['ancillary_variables'] = "pr_localmid_qc" # List other variables associated with variable (QA/QC)
                
                if "pr_raw_qc" in ds.keys():
                    flagvals = list(np.unique(ds['pr_raw_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pr_raw_qc")
                    else:
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                        ds['pr_raw_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                        ds['pr_raw_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        
                        ds['pr_raw'].attrs['ancillary_variables'] = "pr_raw_qc" # List other variables associated with variable (QA/QC)
                        
                if "pr_1h_qc" in ds.keys():
                    flagvals = list(np.unique(ds['pr_1h_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("pr_1h_qc")
                    else:
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                        ds['pr_1h_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                        ds['pr_1h_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        
                        ds['pr_1h'].attrs['ancillary_variables'] = "pr_1h_qc" # List other variables associated with variable (QA/QC)
                    
                # Make pr variable (as NA variable)
                pr = np.full((1, len(ds['time'])), np.nan)
                precip_cols = str([elem for elem in ds.keys() if "pr_" in elem])
                
                ds['pr'] = (['station', 'time'], pr)
                # Set attributes
                ds['pr'].attrs['long_name'] = "precipitation"
                ds['pr'].attrs['standard_name'] = ""
                ds['pr'].attrs['units'] = "mm/hr"
                if precip_cols is not None:
                    ds['pr'].attrs['ancillary_variables'] = precip_cols

                # hurs: relative humidity
                if 'relative_humidity_set_1' in ds.keys():
                    ds = ds.rename({'relative_humidity_set_1': 'hurs'})

                    # If QA/QC column exists
                    if 'relative_humidity_set_1_qc' in ds.keys():
                        flagvals = list(np.unique(ds['relative_humidity_set_1_qc'].values))
                        flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                        if flagvals == ['nan']:
                            #print("Flag value is {}".format(flagvals)) # For testing.
                            ds = ds.drop("relative_humidity_set_1_qc")
                        else:
                            ds = ds.rename({'relative_humidity_set_1_qc': 'hurs_qc'})
                            flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                            ds['hurs_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                            ds['hurs_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                            
                            ds['hurs'].attrs['ancillary_variables'] = "hurs_qc" # List other variables associated with variable (QA/QC)
                    
                    # If column exists and is entirely NA
                    if np.isnan(ds['hurs'].values).all():
                        # Attempt to calculate value using temperature and relative humidity.
                        ds['hurs'] = calc_clean._calc_relhumid(ds['tas'], ds['tdps'])
                        ds['hurs'].attrs['ancillary_variables'] = "tas tdps" # List other variables associated with variable
                        ds['hurs'].attrs['comment'] = "Calculated from air temperature and dew point temperature using calc_clean.py."

                else:
                    hurs = np.full((1, len(ds['time'])), np.nan)
                    ds['hurs'] = (['station', 'time'], hurs)
            
                ds['hurs'].attrs['long_name'] = "relative_humidity"
                ds['hurs'].attrs['standard_name'] = "relative_humidity" 
                ds['hurs'].attrs['units'] = "percent" 
                
                #rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
                
                if 'solar_radiation_set_1' in ds.keys(): # Already in w/m2, no need to convert units.
                    # If column exists, rename.
                    ds = ds.rename({'solar_radiation_set_1': 'rsds'})
                    
                else: # Otherwise, generate column of NAs.
                    rsds = np.full((1, len(ds['time'])), np.nan)
                    ds['rsds'] = (['station', 'time'], rsds)
                
                # Set attributes
                ds['rsds'].attrs['long_name'] = "solar_radiation"
                ds['rsds'].attrs['standard_name'] = "surface_downwelling_shortwave_flux_in_air"
                ds['rsds'].attrs['units'] = "W m-2"

                #rsds: QA/QC flags
                if 'solar_radiation_set_1_qc' in ds.keys():
                    flagvals = list(np.unique(ds['solar_radiation_set_1_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("solar_radiation_set_1_qc")
                    else:
                        ds = ds.rename({'solar_radiation_set_1_qc': 'rsds_qc'})
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                        ds['rsds_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                        ds['rsds_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        
                        ds['rsds'].attrs['ancillary_variables'] = "rsds_qc" # List other variables associated with variable (QA/QC)
                        

                # sfcWind : wind speed (m/s)
                # (Method of calculation may vary and is unknown source by source.)
                # See: https://weather.gladstonefamily.net/CWOP_Guide.pdf
                if "wind_speed_set_1" in ds.keys(): # Data already in m/s.
                    ds = ds.rename({'wind_speed_set_1': 'sfcWind'})
                    ds['sfcWind'].attrs['long_name'] = "wind_speed"
                    ds['sfcWind'].attrs['standard_name'] = "wind_speed"
                    ds['sfcWind'].attrs['units'] = "m s-1"
                    ds['sfcWind'].attrs['comment'] = "Method of wind speed calculation varies, with 2-minute mean as sampling standard."
                
                if 'wind_speed_set_1_qc' in ds.keys():
                    flagvals = list(np.unique(ds['wind_speed_set_1_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    # If only value is NA, drop QAQC column:
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("wind_speed_set_1_qc")
                    else: # Otherwise, rename and reformat.
                        ds = ds.rename({'wind_speed_set_1_qc': 'sfcWind_qc'})
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                    
                        ds['sfcWind_qc'].attrs['flag_values'] = flagvals
                        ds['sfcWind_qc'].attrs['flag_meanings'] = "See QA/QC csv for network." 
                        ds['sfcWind'].attrs['ancillary_variables'] = "sfcWind_qc" # List other variables associated with variable (QA/QC)
        
                
                # sfcWind_dir: wind direction
                if "wind_direction_set_1" in ds.keys(): # No conversions needed, do not make raw column.
                    ds = ds.rename({'wind_direction_set_1': 'sfcWind_dir'})
                    ds['sfcWind_dir'].attrs['long_name'] = "wind_direction"
                    ds['sfcWind_dir'].attrs['standard_name'] = "wind_from_direction"
                    ds['sfcWind_dir'].attrs['units'] = "degrees_clockwise_from_north"
                    
                if 'wind_direction_set_1_qc' in ds.keys():
                    flagvals = list(np.unique(ds['wind_direction_set_1_qc'].values))
                    flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
                    
                    # If only value is NA, drop QAQC column:
                    if flagvals == ['nan']:
                        #print("Flag value is {}".format(flagvals)) # For testing.
                        ds = ds.drop("wind_direction_set_1_qc")
                    else: # Otherwise, rename and reformat.
                        ds = ds.rename({'wind_direction_set_1_qc': 'sfcWind_dir_qc'})
                        flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                    
                        ds['sfcWind_dir_qc'].attrs['flag_values'] = flagvals
                        ds['sfcWind_dir_qc'].attrs['flag_meanings'] = "See QA/QC csv for network." 
                        ds['sfcWind_dir'].attrs['ancillary_variables'] = "sfcWind_dir_qc" # List other variables associated with variable (QA/QC)

                # Other variables: rename to match format
                if 'altimeter_set_1' in ds.keys():
                    ds = ds.rename({'altimeter_set_1': 'ps_altimeter',
                                'altimeter_set_1_qc': 'ps_altimeter_qc'})
                if 'dew_point_temperature_set_1d' in ds.keys():
                    ds = ds.rename({'dew_point_temperature_set_1d':'tdps_derived'})
                if 'pressure_set_1d' in ds.keys():
                    ds = ds.rename({'pressure_set_1d':'ps_derived'})

                if 'ps_altimeter' in ds.keys():
                    ds['ps_altimeter'].attrs['long_name'] = "altimeter"
                    ds['ps_altimeter'].attrs['units'] = "Pa"
                    ds['ps_altimeter'].attrs['ancillary_variables'] = "ps"

                if 'tdps_derived' in ds.keys():
                    ds['tdps_derived'].attrs['long_name'] = "derived_dew_point_temperature"
                    ds['tdps_derived'].attrs['units'] = "Celsius"
                    ds['tdps_derived'].attrs['comment'] = "Derived by Synoptic."

                if 'ps_derived' in ds.keys():
                    ds['ps_derived'].attrs['long_name'] = "derived_station_pressure"
                    ds['ps_derived'].attrs['units'] = "Pa"
                    ds['ps_derived'].attrs['comment'] = "Derived by Synoptic."

                ds = ds.drop("Station_ID") # Drop repeated Station_ID column.

                # Reorder variables
                desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
                rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom.
                new_index = desired_order + rest_of_vars
                ds = ds[new_index]
                

                #Testing: Manually check values to see that they seem correctly scaled, no unexpected NAs.
                # for var in ds.variables:
                #     try:
                #         print([var, float(ds[var].min()), float(ds[var].max())]) 
                #     except:
                #         next
                # Note: here we have tdps and hurs values of 0. This is what causes the "dividing by 0" errors, likely.
                # Address in QA/QC stage (and recalc tdps if needed).

                # Merge file to previous time period data for same station.
                if ds_stat is None:
                    ds_stat = ds # For first file in station, set ds_stat to be original ds.

                else:
                    ds_stat = xr.merge([ds_stat,ds], compat = 'no-conflicts') # Otherwise, merge new records to first ds.
                    # TO BE DONE in AWS: run and test that compat mode correctly joins files?
                    print("Merging records for station {}".format(i)) # Note that we don't add QC here to remove duplicated stations, as this will occur in next step.
                    file_count +=1
                    ds_stat.attrs['raw_files_merged'] = file_count # Keep count of how many files merged per station.

                
            except Exception as e:
                print(e)
                errors['File'].append(station_id)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                  

            # Write station file to netcdf.
            if ds_stat is None: # Should be caught be error handling above, but add in case.
                print("{} not saved.".format(file))
                continue
            
            else:
                try:
                    filename = station_id+".nc" # Make file name
                    filepath = homedir+"/"+savedir+filename # Write file path
                    
                    ds_stat.to_netcdf(path = filepath) # Save station file.
                    print("Saving {} with dims {}".format(filename, ds_stat.dims))
                    ds_stat.close() # Close dataframe.
                except Exception as e:
                    print(filename, e)
                    errors['File'].append(filename)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)
                    continue  

    # Write the list of removed variables to csv for future reference. Keep these in one centralized "removedvars.csv" file that gets appended to
    vars = []
    try:
        os.chdir(homedir+'/'+savedir)
        # Do the reading
        if os.path.isfile('removedvars.csv'): # If removedvars.csv already exists, read in data to vars. Otherwise, vars is empty.
            with open('removedvars.csv','r') as f: 
                rows = list(csv.DictReader(f, delimiter=','))
                for row in rows:
                    vars.append(row['Variable'])
        
        for i in removedvars:
            if i not in vars:
                vars.append(i)
        vars = sorted(vars, key = str.casefold) # Sort and ignore case

        dict = {'Variable': vars}
        
        with open('removedvars.csv', 'w') as f: 
            # pass the csv file to csv.writer function.
            writer = csv.writer(f)

            # pass the dictionary keys to writerow
            # function to frame the columns of the csv file
            writer.writerow(dict.keys())

            # make use of writerows function to append
            # the remaining values to the corresponding
            # columns using zip function.
            writer.writerows(zip(*dict.values()))

    except Exception as e:
        errors['File'].append("N/A")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    #Write errors to csv
    filepath = homedir+"/"+savedir+"errors_cwop_{}.csv".format(end_api) # Set path to save error file.
    #print(errors)
    with open(filepath, "w") as outfile:
        # pass the csv file to csv.writer function.
        writer = csv.writer(outfile)

        # pass the dictionary keys to writerow
        # function to frame the columns of the csv file
        writer.writerow(errors.keys())

        # make use of writerows function to append
        # the remaining values to the corresponding
        # columns using zip function.
        writer.writerows(zip(*errors.values()))
    
# Run function
#clean_cwop(homedir, workdir, savedir)

# Testing:
## Import file.
os.chdir(savedir)
test = xr.open_dataset("CWOP_AP741.nc")
print(test.head())

for var in test.keys():
    print(var)
    print(test[var])

# # # # Test 1: multi-year merges work as expected.
print(str(test['time'].min())) # Get start time
print(str(test['time'].max())) # Get end time


# # ## Test 2: Inspect vars and attributes
# # ## 
for var in test.variables: 
    try:
        print([var, float(test[var].min()), float(test[var].max())]) 
    except:
        continue

# # # Test 3: Get one month's data and test subsetting.
print(test.sel(time = "2015-05"))

# # # Next file.
test = xr.open_dataset("CWOP_D3845.nc") 
print(test)
print(str(test['time'].min())) # Get start time
print(str(test['time'].max())) # Get end time

## Inspect vars and attributes
## 
for var in test.variables: 
    try:
        print([var, float(test[var].min()), float(test[var].max())]) 
    except:
        continue

# # Get a few rows
print(test.sel(time = "2010-06")) 