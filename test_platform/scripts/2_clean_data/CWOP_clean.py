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
# savedir: path to where cleaned files should be saved.
# Output: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
# Note: files are already filtered by bbox in the pull step, so additional geospatial filtering is not required here.

def clean_cwop(homedir, workdir, savedir):

    os.chdir(homedir)
    os.chdir(workdir)
    
    files = os.listdir() # Gets list of files in directory to work with
    files = list(filter(lambda f: f.endswith(".csv"), files)) # Get list of file names
    
    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
    timestamp = datetime.now().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.

    # Iterate through files to clean. Each file represents a station's data.
    for file in files:
    #for id in ['C4951', 'C2512']: # For testing, pick a few stations.
        station_id = "CWOP_"+file.replace(".csv", "")
        
        #try:
        with open(file, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            for index, row in enumerate(csv_reader):
                
                # Parse header information from CSV.
                if len(row)==1:
                    row = str(row)
                    
                    if "STATION NAME" in row:
                        station_name = row.partition(": ")[2].replace("']", "")
                    elif "LATITUDE" in row:
                        latitude = float(row.partition(": ")[2].replace("']", ""))
                    elif "LONGITUDE" in row:
                        longitude = float(row.partition(": ")[2].replace("']", ""))
                    elif "ELEVATION" in row:
                        elevation = float(row.partition(": ")[2].replace("']", "")) # In feet.
                    elif "STATE" in row:
                        state = row.partition(": ")[2].replace("']", "")
                else:
                    # Get column names
                    if "Station_ID" in row:
                        columns = row
                        
                    # Get units
                    unitstocheck = ['Pascals', '%', 'm/s', 'Celsius', 'QC_Type', 'Degrees']
                    if any(unit in row for unit in unitstocheck):
                        units_index = index
                        
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
                    print("Dropping duplicate column {}".format(cols[1]))
                    df.drop(cols[1], axis = 1, inplace = True)
                    #print(df.columns)
                else:
                    print("Non-identical duplicate columns found.")
            
        
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

        # Add variable: elevation
        elev = np.asarray([elevation]*len(ds['time']))
        elev.shape = (1, len(ds['time']))
        ds['elevation'] = (['station', 'time'], elev)
        
        # Update dimension and coordinate attributes.
        # Time
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
            ds['tas'].attrs['ancillary_variables'] = "tas_raw tas_qc" # List other variables associated with variable (QA/QC)
            ds['tas'].attrs['comment'] = "Converted from Celsius."
            
            ds['tas_raw'].attrs['long_name'] = "air_temperature"
            ds['tas_raw'].attrs['standard_name'] = "air_temperature"
            ds['tas_raw'].attrs['units'] = "degree_Celsius"
            ds['tas_raw'].attrs['ancillary_variables'] = "tas tas_qc" # List other variables associated with variable (QA/QC)
        
        if "air_temp_set_1_qc" in ds.keys():
            # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
            
            flagvals = list(np.unique(ds['air_temp_set_1_qc'].values))
            flagvals = list(np.unique([split_item for item in flagvals for split_item in str(item).split(";")])) # Split any rows with multiple flags and run unique again.
            
            if len(flagvals)==1:
                print("Flag value is {}".format(flagvals)) # For testing.
                df.drop("air_temp_set_1_qc")
            else:
                ds = ds.rename({'air_temp_set_1_qc': 'tas_qc'})
                flagvals = [ele.replace(".0", "") for ele in flagvals] # Reformat to match csv.
                ds['tas_qc'].attrs['flag_values'] = flagvals # Generate values from unique values from dataset.
                ds['tas_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                
        # ps: surface air pressure (Pa)
        # Note here that if "pressure_set_1" has values this is a direct station observation reading.
        # Otherwise, if "pressure_set_1d" has values this is a derived value calculated from altimeter and elevation.
        # We will manually recalculate this here.

        ## Left off here on 9.05.

        if '' in ds.keys():
            ds['ps_raw'] = ds['ps'] # Move original data to raw column.
            ds['ps'].attrs['long_name'] = "station_air_pressure"
            ds['ps'].attrs['standard_name'] = "air_pressure"
            ds['ps'].attrs['units'] = "hPa"
            ds['ps'].attrs['ancillary_variables'] = "ps_raw ps_qc" # List other variables associated with variable (QA/QC)
    
            ds['ps_raw'].attrs['long_name'] = "station_air_pressure"
            ds['ps_raw'].attrs['standard_name'] = "air_pressure"
            ds['ps_raw'].attrs['units'] = "hPa"
            ds['ps_raw'].attrs['ancillary_variables'] = "ps ps_qc ps_altimeter ps_altimeter_qc" # List other variables associated with variable (QA/QC)

            # Delete sea level pressure if station pressure included.
            ds.drop(("psl", "psl_qc"))
        
        next

        

        #except Exception as e:
        #    print(e)

    
            # Update variable attributes and do unit conversions

    #             
    #             

    #             if "ps_qc" in ds.keys():
    #                 ds['ps_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, M, 9]"
    #                 ds['ps_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
        
    #             if "ps_altimeter" in ds.keys():
    #                 ds['ps_altimeter'].attrs['long_name'] = "altimeter_setting"
    #                 ds['ps_altimeter'].attrs['units'] = "hPa"
    #                 ds['ps_altimeter'].attrs['ancillary_variables'] = "ps ps_qc ps_altimeter ps_altimeter_qc" # List other variables associated with variable (QA/QC)
    #                 ds['ps_altimeter'].attrs['comment'] = "The pressure value to which an aircraft altimeter is set so that it will indicate the altitude relative to mean sea level of an aircraft on the ground at the location for which the value was determined." # Description of variable meaning, as not CF-standard.
                
    #             if "ps_altimeter_qc" in ds.keys():
    #                 ds['ps_altimeter_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, M, 9]"
    #                 ds['ps_altimeter_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."

    #             # If station air pressure not reported, convert from sea level pressure.
    #             # Inputs: sea level pressure (mb/hPa), elevation (m), and air temperature (K)
    #             elif 'psl' in ds.keys():
    #                 ds['psl_raw'] = ds['ps'] # Move original data to raw column.
                    
    #                 # Convert from sea level pressure. TO DO. (Need to know elev of barometer.)
    #                 try:
    #                     ds['ps'] = calc_clean._calc_ps(ds['psl'], ds['elevation'], ds['tas'])
    #                 except:
    #                     print("ps: calc_clean.py not working.") 
    #                     ds['ps'] = ds['psl']*math.e**(-ds['elevation']/(ds['tas']*29.263)) # Backup method, remove once calc_clean works consistently?

    #                 ds['ps'].attrs['long_name'] = "station_air_pressure"
    #                 ds['ps'].attrs['standard_name'] = "air_pressure"
    #                 ds['ps'].attrs['units'] = "hPa"
    #                 ds['ps'].attrs['ancillary_variables'] = "psl_raw psl_qc elevation" # List other variables associated with variable (QA/QC)
    #                 ds['ps'].attrs['comment'] = "Converted from sea level pressure."
                    
    #                 ds['psl_raw'].attrs['long_name'] = "sea_level_air_pressure"
    #                 ds['psl_raw'].attrs['standard_name'] = "air_pressure_at_sea_level"
    #                 ds['psl_raw'].attrs['units'] = "hPa"
    #                 ds['psl_raw'].attrs['ancillary_variables'] = "ps psl_qc" # List other variables associated with variable (QA/QC)
                
    #             # tdps: dew point temperature (K)
    #             # tdps always provided by ISD reports, so no need to calculate here.
    #             if "tdps" in ds.keys(): # If variable already exists, rename.
    #                 ds['tdps_raw'] = ds['tdps'] # Move original data to raw column.

    #                 try:
    #                     ds['tdps'] = calc_clean._unit_degC_to_K(ds['tdps'])
    #                 except:
    #                     print("tdps: calc_clean.py not working.") 
    #                     ds['tdps'] = ds["tdps"]+273.15 # Convert from C to K
    #                 ds['tdps'].attrs['long_name'] = "dew_point_temperature"
    #                 ds['tdps'].attrs['standard_name'] = "dew_point_temperature"
    #                 ds['tdps'].attrs['units'] = "degree_Kelvin"
    #                 ds['tdps'].attrs['ancillary_variables'] = "tdps_raw tdps_qc" # List other variables associated with variable (QA/QC)
    #                 ds['tdps'].attrs['comment'] = "Converted from Celsius."
                    
    #                 ds['tdps_raw'].attrs['long_name'] = "dew_point_temperature"
    #                 ds['tdps_raw'].attrs['standard_name'] = "dew_point_temperature"
    #                 ds['tdps_raw'].attrs['units'] = "degree_Celsius"
    #                 ds['tdps_raw'].attrs['ancillary_variables'] = "tdps tdps_qc" # List other variables associated with variable (QA/QC)

    #             if 'tdps_qc' in ds.keys():
    #                 ds['tdps_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9, A, C, I, M, P, R, U]"
    #                 ds['tdps_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
    #             # pr: precipitation (Flag for discussion about CF compliance and units) 
    #             if 'pr' in ds.keys():
    #                 ds['pr_raw'] = ds['pr'] # Move original data to raw column.
    #                 # Note leave final conversion here to next stage.
    #                 ds['pr'].attrs['long_name'] = "precipitation"
    #                 ds['pr'].attrs['standard_name'] = "" # Flag for discussion.
    #                 ds['pr'].attrs['units'] = "" # Flag for discussion.
    #                 ds['pr'].attrs['ancillary_variables'] = "pr_raw pr_qc pr_depth_qc pr_duration" # List other variables associated with variable (QA/QC)
    #                 ds['pr'].attrs['comment'] = "" # To be completed.

    #                 ds['pr_raw'].attrs['long_name'] = "precipitation"
    #                 ds['pr_raw'].attrs['standard_name'] = "" # TBD. Maybe delete, not standard.
    #                 ds['pr_raw'].attrs['units'] = "mm"
    #                 ds['pr_raw'].attrs['ancillary_variables'] = "pr pr_qc pr_depth_qc pr_duration" # List other variables associated with variable (QA/QC)

    #                 ds['pr_duration'].attrs['long_name'] = "precipitation measurement interval"
    #                 ds['pr_duration'].attrs['units'] = "hours"
    #                 ds['pr_duration'].attrs['ancillary_variables'] = "pr pr_raw pr_qc pr_depth_qc" # List other variables associated with variable (QA/QC)

    #             if 'pr_qc' in ds.keys():
    #                 ds['pr_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9, A, I, M, P, R, U]"
    #                 ds['pr_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
    #                 ds['pr_depth_qc'].attrs['flag_values'] = "[1, 2, 3, 4, 5, 6, 7, 8, E, I, J, 9]"
    #                 ds['pr_depth_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                
    #             # hurs: relative humidity
    #             if 'hurs' in ds.keys():
    #                 ds['hurs'].attrs['long_name'] = "average_relative_humidity"
    #                 ds['hurs'].attrs['standard_name'] = "relative_humidity" 
    #                 ds['hurs'].attrs['units'] = "percent" 
    #                 ds['hurs'].attrs['ancillary_variables'] = "hurs_qc hurs_flag hurs_duration hurs_temp hurs_temp_qc hurs_temp_flag" # List other variables associated with variable (QA/QC)

    #                 ds['hurs_duration'].attrs['long_name'] = "relative_humidity_duration"
    #                 ds['hurs_duration'].attrs['units'] = "minutes" 
    #                 ds['hurs_duration'].attrs['ancillary_variables'] = "hurs hurs_qc hurs_flag hurs_temp hurs_temp_qc hurs_temp_flag" # List other variables associated with variable (QA/QC)

    #                 ds['hurs_temp'].attrs['long_name'] = "average_air_temperature_at_hurs_instrument"
    #                 ds['hurs_temp'].attrs['units'] = "degree_celsius" 
    #                 ds['hurs_temp'].attrs['ancillary_variables'] = "hurs hurs_qc hurs_flag hurs_temp hurs_temp_qc hurs_temp_flag" # List other variables associated with variable (QA/QC)
                
    #             elif 'tas' and 'tdps' in ds.keys():
    #                 try:
    #                     ds['hurs'] = calc_clean._calc_relhumid(ds['tas'], ds['tdps'])
    #                 except:
    #                     print("hurs: calc_clean.py not working.")
    #                     es = 0.611 * np.exp(5423 * ((1/273) - (1/tas)))
    #                     e = 0.611 * np.exp(5423 * ((1/273) - (1/tdps)))
    #                     ds['hurs'] = 100 * (e/es)
    #                 ds['hurs'].attrs['long_name'] = "relative_humidity" # Note some here will be avg. and some not, will have to align when we merge down.
    #                 ds['hurs'].attrs['standard_name'] = "relative_humidity" 
    #                 ds['hurs'].attrs['units'] = "percent" 
    #                 ds['hurs'].attrs['comment'] = "Calculated from air temperature and dew point temperature."

    #             if 'hurs_qc' in ds.keys():
    #                 ds['hurs_qc'].attrs['flag_values'] = "[1, 3, 9]"
    #                 ds['hurs_qc'].attrs['flag_meanings'] =  "passed_all_qc_checks failed_all_qc_checks missing"
                    
    #                 ds['hurs_flag'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]"
    #                 ds['hurs_flag'].attrs['flag_meanings'] =  "See QA/QC csv for network."

    #                 ds['hurs_temp_qc'].attrs['flag_values'] = "[1, 3, 9]"
    #                 ds['hurs_temp_qc'].attrs['flag_meanings'] =  "passed_all_qc_checks failed_all_qc_checks missing"
                    
    #                 ds['hurs_temp_flag'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]"
    #                 ds['hurs_temp_flag'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
    #             # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
    #             if 'rsds' in ds.keys():
    #                 #print(ds['rsds'].values[0])
    #                 if np.isnan(ds['rsds'].values[0]).all(): # If all values are NA
    #                     pass # Continue and don't add this variable.
    #                 else: # Otherwise, add variables.
    #                     ds['rsds'].attrs['long_name'] = "solar_radiation"
    #                     ds['rsds'].attrs['standard_name'] = "surface_downwelling_shortwave_flux_in_air"
    #                     ds['rsds'].attrs['units'] = "W m-2"
    #                     ds['rsds'].attrs['ancillary_variables'] = "rsds_duration rsds_qc rsds_flag" # List other variables associated with variable (QA/QC)
    #                     ds['rsds'].attrs['comment'] = "Waveband ranges from 0.4 - 2.3 micrometers."

    #                     ds['rsds_duration'].attrs['long_name'] = "time_period"
    #                     ds['rsds_duration'].attrs['units'] = "minutes"
    #                     ds['rsds_duration'].attrs['ancillary_variables'] = "rsds rsds_qc rsds_flag" # List other variables associated with variable (QA/QC)
    #                     ds['rsds_duration'].attrs['comment'] = "Time period over which solar radiation is integrated."

    #                 if 'rsds_qc' in ds.keys():
    #                     ds['rsds_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 9]"
    #                     ds['rsds_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                    
    #                     ds['rsds_flag'].attrs['flag_values'] = "00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10-93, 94-97, 98, 99"
    #                     ds['rsds_flag'].attrs['flag_meanings'] =  "See QA/QC csv for network."

                   
    #             # sfcWind : wind speed (m/s) (Method of calculation may vary, standardize during hourly merge or QA/QC process.)
    #             if "sfcWind" in ds.keys(): # No conversions needed, do not add raw column.
    #                 ds['sfcWind'].attrs['long_name'] = "wind_speed"
    #                 ds['sfcWind'].attrs['standard_name'] = "wind_speed"
    #                 ds['sfcWind'].attrs['units'] = "m s-1"
    #                 ds['sfcWind'].attrs['ancillary_variables'] = "sfcWind_qc sfcWind_method" # List other variables associated with variable (QA/QC)
    #                 ds['sfcWind'].attrs['comment'] = "Method of wind speed calculation varies, see sfcWind_method."

    #             if 'sfcWind_method' in ds.keys():
    #                 ds['sfcWind_method'].attrs['long_name'] = "wind_speed_calculation_method"
    #                 ds['sfcWind_method'].attrs['flag_values'] = "[A, B, C, H, N, R, Q, T, V, 9]"
    #                 ds['sfcWind_method'].attrs['flag_meanings'] = "abridged_beaufort beaufort calm 5-minute_average normal 60-minute_average squall 180-minute_average variable missing" 

    #             if 'sfcWind_qc' in ds.keys():
    #                 ds['sfcWind_method'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9]"
    #                 ds['sfcWind_method'].attrs['flag_meanings'] = "See QA/QC csv for network." 


    #             # sfcWind_dir: wind direction
    #             if "sfcWind_dir" in ds.keys(): # No conversions needed, do not make raw column.
    #                 ds['sfcWind_dir'].attrs['long_name'] = "wind_direction"
    #                 ds['sfcWind_dir'].attrs['standard_name'] = "wind_from_direction"
    #                 ds['sfcWind_dir'].attrs['units'] = "degrees_clockwise_from_north"
    #                 ds['sfcWind_dir'].attrs['ancillary_variables'] = "sfcWind_dir_qc" # List other variables associated with variable (QA/QC)
                    
    #             if 'sfcWind_dir_qc' in ds.keys():
    #                 ds['sfcWind_dir_qc'].attrs['flag_values'] = "[0, 1, 2, 3, 4, 5, 6, 7, 9]"
    #                 ds['sfcWind_dir_qc'].attrs['flag_meanings'] = "See QA/QC csv for network." 


    #             # Update attributes for any non-standard variables
                
    #             #QAQC process
    #             if "qaqc_process" in ds.keys():
    #                 ds['qaqc_process'].attrs['long_name'] = "qaqc_process_type"
    #                 ds['qaqc_process'].attrs['flag_values'] = "[V01, V02, V03]"
    #                 ds['qaqc_process'].attrs['flag_meanings'] = "no_qaqc automated_qaqc subjected_to_qaqc" 

    #             # Data source
    #             if "data_source" in ds.keys():
    #                 ds['data_source'].attrs['long_name'] = "source_of_data"
    #                 ds['data_source'].attrs['flag_values'] = "[6, 7]"
    #                 ds['data_source'].attrs['flag_meanings'] = "ASOSAWOS ASOSAWOS_merged_with_USAF_surface_hourly" 

    #             # # Reorder variables
    #             # In following order:
    #             desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
    #             desired_order = [i for i in desired_order if i in list(ds.keys())] # Only keep vars which are in ds.
    #             rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom.
    #             new_index = desired_order + rest_of_vars
    #             ds = ds[new_index]

    #             # Testing: Manually check values to see that they seem correctly scaled, no unexpected NAs.
    #             # for var in ds.variables:
    #             #     try:
    #             #         print([var, float(ds[var].min()), float(ds[var].max())]) 
    #             #     except:
    #             #         next
                
    #         except Exception as e: # If error in xarray reorganization
    #             print(file, e)
    #             errors['File'].append(station_id)
    #             errors['Time'].append(end_api)
    #             errors['Error'].append(e)
                    

    #         # Write station file to netcdf.
    #         if ds is None: # Should be caught be error handling above, but add in case.
    #             print("ds {} not saved.".format(file))
    #             continue
            
    #         else:
    #             try:
    #                 #print(ds) # For testing.
    #                 filename = station_id+".nc" # Make file name
    #                 filepath = homedir+"/"+savedir+filename # Write file path
    #                 #print(filepath)

    #                 ds.to_netcdf(path = filepath) # Save station file.
    #                 print("Saving {} with dims {}".format(filename, ds.dims))
    #                 ds.close() # Close dataframe.
    #             except Exception as e:
    #                 print(filename, e)
    #                 errors['File'].append(filename)
    #                 errors['Time'].append(end_api)
    #                 errors['Error'].append(e)
    #                 continue            

    # #Write errors to csv
    # filepath = homedir+"/"+savedir+"errors_asosawos_{}.csv".format(end_api) # Set path to save error file.
    # #print(filepath)
    # #print(errors)
    # with open(filepath, "w") as outfile:
    #     # pass the csv file to csv.writer function.
    #     writer = csv.writer(outfile)

    #     # pass the dictionary keys to writerow
    #     # function to frame the columns of the csv file
    #     writer.writerow(errors.keys())

    #     # make use of writerows function to append
    #     # the remaining values to the corresponding
    #     # columns using zip function.
    #     writer.writerows(zip(*errors.values()))

# Run function
clean_cwop(homedir, workdir, savedir)

# Testing:
## Import file.
# os.chdir(savedir)
# test = xr.open_dataset("ASOSAWOS_72061999999.nc")
# print(test) 

# # # # Test 1: multi-year merges work as expected.
# print(str(test['time'].min())) # Get start time
# print(str(test['time'].max())) # Get end time


# # ## Test 2: Inspect vars and attributes
# # ## 
# for var in test.variables: 
#     try:
#         print([var, float(test[var].min()), float(test[var].max())]) 
#     except:
#         continue

# # # Test 3: Get one month's data and test subsetting.
# print(test.sel(time = "2015-05"))

# # # Next file.
# test = xr.open_dataset("ASOSAWOS_72268999999.nc") 
# print(test)
# print(str(test['time'].min())) # Get start time
# print(str(test['time'].max())) # Get end time

# ## Inspect vars and attributes
# ## 
# for var in test.variables: 
#     try:
#         print([var, float(test[var].min()), float(test[var].max())]) 
#     except:
#         continue

# # # Get a few rows
# print(test.sel(time = "1989-06")) # Should only be 3 observations here. Very incomplete year.