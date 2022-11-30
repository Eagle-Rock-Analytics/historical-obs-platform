"""
This script performs data cleaning for CIMIS stations pulled through CADWR for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables and drop unnecessary variables
(2) Converts station metadata to standard format, with unique identifier
(3) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(4) Converts missing data to standard format
(5) Tracks existing qa/qc flag for review
(6) Merge files by station, and outputs cleaned variables as a single .nc file for each station in an individual network.
Inputs: Raw data for the network's stations, with each csv file representing data for all stations in a given year or month (current year's data)
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
Reference: https://www.ncei.noaa.gov/data/global-hourly/doc/isd-format-document.pdf

Note: QAQC flags formatted and uploaded manually. Last update Nov 22 2022.
Source: https://cimis.water.ca.gov/Content/PDF/CurrentFlags2.pdf (1995-present)
Source: https://cimis.water.ca.gov/Content/PDF/FormerFlags2.pdf (Pre 1995)
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
import zipfile
import openpyxl
from ftplib import FTP
from cleaning_helpers import var_to_unique_list

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


# Given a network name, return all relevant AWS filepaths for other functions.
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir

              
# ## FUNCTION: Clean CIMIS data.
# # Input: 
# # bucket_name: name of AWS bucket.
# # rawdir: path to where raw data is saved as .csv files, with each file representing a station's records from download start date to present (by default 01-01-1980).
# # cleandir: path to where cleaned files should be saved.
# # Output: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
# # Note: files are already filtered by bbox in the pull step, so additional geospatial filtering is not required here.

# # Function: take heads, read csv into pandas db and clean.
def clean_cimis(rawdir, cleandir):
    network = "CIMIS"

    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
    timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.
    
    # Column format changes in June 2014 for this dataset.
    newcols = ['Station ID', 'Date', 'Hour (PST)', 'Julian Date', 'Reference ETo (mm)', 'QC for Reference ETo',
                    'Precipitation (mm)', 'QC for Precipitation', 'Solar Radiation (W/m²)', 'QC for Solar Radiation',
                    'Vapor Pressure (kPa)', 'QC for Vapor Pressure', 'Air Temperature (°C)', 'QC for Air Temperature',
                    'Relative Humidity (%)', 'QC for Relative Humidity', 'Dew Point (°C)', 'QC for Dew Point', 
                    'Wind Speed (m/s)', 'QC for Wind Speed', 'Wind Direction (0-360)', 'QC for Wind Direction',
                    'Soil Temperature (°C)', 'QC for Soil Temperature']
    oldcols = ['Station ID', 'Date', 'Hour (PST)', 'Julian Date', 'QC for Reference ETo', 'Reference ETo (mm)',
                    'QC for Precipitation', 'Precipitation (mm)', 'QC for Solar Radiation', 'Solar Radiation (W/m²)',
                    'QC for Vapor Pressure', 'Vapor Pressure (kPa)', 'QC for Air Temperature', 'Air Temperature (°C)', 
                    'QC for Relative Humidity', 'Relative Humidity (%)', 'QC for Dew Point', 'Dew Point (°C)', 
                    'QC for Wind Speed', 'Wind Speed (m/s)', 'QC for Wind Direction', 'Wind Direction (0-360)', 
                    'QC for Soil Temperature', 'Soil Temperature (°C)']

    # Specify columns to remove
    removecols = ['Julian Date', 'QC for Soil Temperature', 'Soil Temperature (°C)', 'Reference ETo (mm)', 'QC for Reference ETo']

    try:
        # Get files
        files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir): 
            file = str(item.key)
            files += [file]

        # Get station file and read in metadata.
        station_file = [file for file in files if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_excel(BytesIO(obj['Body'].read()))    
        stations = station_file['Station Number'].dropna().astype(int)
        
        # Remove error, station files
        files = [file for file in files if '.zip' in file]
        files = [file for file in files if 'stationlist' not in file]
        files = [file for file in files if 'error' not in file]

    except Exception as e: # If unable to read files from rawdir, break function.
        print(e)
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    else: # If files read successfully, continue.
        for station in stations.sample(5): # subset for testing
            station_metadata = station_file.loc[station_file['Station Number']==float(station)]
            station_id = "CIMIS_"+str(station)
            
            dfs = []
            for file in files: # For each zip file (annual or monthly)
                try:
                    fileyear = file.split("/")[-1]
                    fileyear = fileyear.replace("hourlyStns", "")
                    # Set default columns to be columns for pre-June 2014 data.
                    if fileyear[0:4].isnumeric():
                        if int(fileyear[0:4])>=2014: # If data from 2014 on
                            columns = newcols
                        else:
                            columns = oldcols
                    else: # If data from current year
                        columns = newcols
                    
                    obj = s3.Bucket(bucket_name).Object(file)
                    with BytesIO(obj.get()["Body"].read()) as tf:
                        # rewind the file
                        tf.seek(0)
                        # Read the file as a zipfile and process the members
                        with zipfile.ZipFile(tf, mode='r') as zipf: # Unzip
                            for subfile in zipf.namelist(): # For each csv file in a zipfile
                                stationid = subfile.replace(".csv", '')
                                stationid = int(stationid[-3:])
                                if station == stationid:
                                    df= pd.read_csv(zipf.open(subfile), names = columns, skipinitialspace= True, na_values = ["*", "--"]) # Read into pandas
                                    # Fix non-standard NAs and whitespace issues while reading in.
                                    
                                    # Reorder columns into new column order
                                    df = df.reindex(columns = newcols)

                                    # Drop columns
                                    df = df.drop(columns = removecols, axis = 1)
                                    
                                    # Fix time format. Dataset reports 1-24h time, while datetime requires 0-23h.
                                    # Convert 24h to 0h, and assign to the day of the following date.

                                    df['Hour'] = df['Hour (PST)'].astype(str).str.zfill(4)
                                    df['Hour'] = [x[0:2]+":"+x[2:] for x in df['Hour']]
                                    df['Hour'] = np.where(df['Hour'] == "24:00", "00:00", df['Hour'])

                                    df['Date'] = pd.to_datetime(df['Date'])
                                    df['Date'] = np.where(df['Hour'] == "00:00", df['Date']+timedelta(days=1), df['Date'])

                                    df['time'] = pd.to_datetime(df['Date'].dt.strftime("%Y-%m-%d") + df['Hour'], format = "%Y-%m-%d%H:%M")
                                    # tz_localize expects PST and PDT time inputs, but we just have PST here so we take a timedelta to UTC.
                                    df['time'] = df['time']+timedelta(hours = 8)
                                    df['time'] = df['time'].dt.tz_localize('UTC')

                                    # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022.
                                    df = df.loc[(df['time']<'2022-09-01') & (df['time']>'1979-12-31')]
                                    
                                    # Remove unnecessary columns.
                                    df = df.drop(['Hour', 'Date', 'Hour (PST)'], axis = 1)
                                    
                                    dfs.append(df)
                                else:
                                    continue

                except Exception as e: # Handle exceptions thrown during individual file read in.
                    print(e)
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)
                    continue


            try:
                file_count = len(dfs)
                df_stat = pd.concat(dfs)
                
                # Drop any columns that only contain NAs.
                df_stat = df_stat.dropna(axis = 1, how = 'all')
                                                
                
                # Fix issue with "nan" and nan causing comparison errors
                df_stat = df_stat.replace("nan", np.nan)

                # Replace non-standard NAs
                df_stat = df_stat.replace(-9998, np.nan)
                df_stat = df_stat.replace(-9997, np.nan)
                df_stat = df_stat.replace(-6999, np.nan)
                df_stat = df_stat.replace(6999, np.nan)
                df_stat = df_stat.replace(-8484, np.nan)

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
                ds = ds.assign_attrs(history = "CIMIS_clean.py script run on {} UTC".format(timestamp))
                ds = ds.assign_attrs(comment = "Intermediate data product: may not have been subject to any cleaning or QA/QC processing")
                ds = ds.assign_attrs(license = "")
                ds = ds.assign_attrs(citation = "")
                ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")
                
                # Add station metadata
                # Station name
                ds = ds.assign_attrs(station_name = station_metadata['Name'].values[0]) 

                # Sensor heights
                ds = ds.assign_attrs(pyranometer_height_m = 2)
                ds = ds.assign_attrs(wind_vane_height_m = 2)
                ds = ds.assign_attrs(anemometer_height_m = 2)
                ds = ds.assign_attrs(air_temperature_height_m = 1.5)
                ds = ds.assign_attrs(humidity_height_m = 1.5)
                ds = ds.assign_attrs(rain_gauge_height_m = 1)

                ds = ds.assign_attrs(raw_files_merged = file_count) # Keep count of how many files merged per station.

                # Add dimensions: station ID and time.
                ds = ds.set_coords('time').swap_dims({'index': 'time'}) # Swap index with time.
                ds = ds.assign_coords(id = str(station_id))
                ds = ds.expand_dims("id") # Add station_id as index.
                ds = ds.drop_vars(("index")) # Drop station_id variable and index coordinate.
                ds = ds.rename({'id': 'station'}) # Rename id to station_id.
                ds = ds.drop('Station ID') # Drop station ID column

                # Update dimensions and coordinates
                                    
                # Add coordinates: latitude and longitude.
                lat = np.asarray([station_metadata['Latitude']]*len(ds['time']))
                lat.shape = (1, len(ds['time']))

                lon = np.asarray([station_metadata['Longitude']]*len(ds['time']))
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
                elev = np.asarray([station_metadata['ELEV']]*len(ds['time']))
                elev.shape = (1, len(ds['time']))
                ds['elevation'] = (['station', 'time'], elev)

                # Update dimension and coordinate attributes.
            
                # Update attributes.
                ds['time'] = pd.to_datetime(ds['time']) # Remove timezone data from string (to match other networks)
                ds['time'] = pd.to_datetime(ds['time'], unit = 'ns')
                ds['time'].attrs['long_name'] = "time"
                ds['time'].attrs['standard_name'] = "time"
                ds['time'].attrs['comment'] = "Converted from PST to UTC."
                
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
                # for var in list(ds.keys()):
                #     try:
                #         print([var, float(ds[var].min()), float(ds[var].max())])
                #     except:
                #         continue

                
                #tas: air surface temperature (K)
                if "Air Temperature (°C)" in ds.keys():
                    ds['tas'] = calc_clean._unit_degC_to_K(ds['Air Temperature (°C)'])
                    ds = ds.drop("Air Temperature (°C)")

                    ds['tas'].attrs['long_name'] = "air_temperature"
                    ds['tas'].attrs['standard_name'] = "air_temperature"
                    ds['tas'].attrs['units'] = "degree_Kelvin"
                    
                    if "QC for Air Temperature" in ds.keys():
                        # Flag values are listed in this column and separated with ; when more than one is used for a given observation.
                        ds = ds.rename({'QC for Air Temperature': 'tas_qc'})
                        ds['tas_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'tas_qc')
                        ds['tas_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."
                        ds['tas'].attrs['ancillary_variables'] = "tas_qc" # List other variables associated with variable (QA/QC)

                    ds['tas'].attrs['comment'] = "Converted from Celsius to Kelvin."

                # ps: surface air pressure (Pa)
                # No pressure sensors in this dataset.
                                                
                # tdps: dew point temperature (K)
                # This is calculated by CIMIS from vapor pressure and air temperature data, not a raw observation.
                if 'Dew Point (°C)' in ds.keys():
                    ds['tdps_derived'] = calc_clean._unit_degC_to_K(ds['Dew Point (°C)'])
                    ds = ds.drop("Dew Point (°C)")

                    # Set attributes for conversion.
                    ds['tdps_derived'].attrs['long_name'] = "dew_point_temperature"
                    ds['tdps_derived'].attrs['standard_name'] = "dew_point_temperature"
                    ds['tdps_derived'].attrs['units'] = "degree_Kelvin"

                    if 'QC for Dew Point' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'QC for Dew Point': 'tdps_derived_qc'})
                        ds['tdps_derived_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'tdps_derived_qc')
                        ds['tdps_derived_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."                  
                        ds['tdps_derived'].attrs['ancillary_variables'] = "tdps_derived_qc" # List other variables associated with variable (QA/QC)
                    
                    ds['tdps_derived'].attrs['comment'] = "Derived by CIMIS from vapor pressure and air temperature. Converted from Celsius to Kelvin."
                    
                # pr: precipitation
                if "Precipitation (mm)" in ds.keys():
                    ds = ds.rename({"Precipitation (mm)":"pr"})
                    ds['pr'].attrs['long_name'] = "precipitation_accumulation"
                    ds['pr'].attrs['units'] = "mm/hour"

                    if 'QC for Precipitation' in ds.keys(): # If QA/QC exists.
                        ds = ds.rename({'QC for Precipitation': 'pr_qc'})
                        ds['pr_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'pr_qc')
                        ds['pr_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."                  
                        ds['pr'].attrs['ancillary_variables'] = "pr_qc" # List other variables associated with variable (QA/QC)
                    
                    ds['pr'].attrs['comment'] = "Accumulated precipitation. Converted from inches to mm."
                
                                
                # hurs: relative humidity (%)
                if 'Relative Humidity (%)' in ds.keys(): # Already in %, no need to convert units.
                    ds = ds.rename({'Relative Humidity (%)': 'hurs'})
                    # Set attributes
                    ds['hurs'].attrs['long_name'] = "relative_humidity"
                    ds['hurs'].attrs['standard_name'] = "relative_humidity" 
                    ds['hurs'].attrs['units'] = "percent" 
                    
                    # If QA/QC column exists
                    if 'QC for Relative Humidity' in ds.keys():
                        ds = ds.rename({'QC for Relative Humidity': 'hurs_qc'})
                        ds['hurs_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'hurs_qc')
                        ds['hurs_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."                     
                        ds['hurs'].attrs['ancillary_variables'] = "hurs_qc" # List other variables associated with variable (QA/QC)
                    
                #rsds: surface_downwelling_shortwave_flux_in_air (solar radiation, w/m2)
                
                if 'Solar Radiation (W/m²)' in ds.keys(): # Already in w/m2, no need to convert units.
                    # If column exists, rename.
                    ds = ds.rename({'Solar Radiation (W/m²)': 'rsds'})
                    
                    # Set attributes
                    ds['rsds'].attrs['long_name'] = "solar_radiation"
                    ds['rsds'].attrs['standard_name'] = "surface_downwelling_shortwave_flux_in_air"
                    ds['rsds'].attrs['units'] = "W m-2"

                    #rsds: QA/QC flags
                    if 'QC for Solar Radiation' in ds.keys():
                        ds = ds.rename({'QC for Solar Radiation': 'rsds_qc'})
                        ds['rsds_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'rsds_qc')
                        ds['rsds_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."                     
                        ds['rsds'].attrs['ancillary_variables'] = "rsds_qc" # List other variables associated with variable (QA/QC)
                
                
                # sfcWind : wind speed (m/s)
                if "Wind Speed (m/s)" in ds.keys(): # Data originally in mph.
                    ds = ds.rename({'Wind Speed (m/s)':"sfcWind"})
                    ds['sfcWind'].attrs['long_name'] = "wind_speed"
                    ds['sfcWind'].attrs['standard_name'] = "wind_speed"
                    ds['sfcWind'].attrs['units'] = "m s-1"

                    #sfcWind: QA/QC flags
                    if 'QC for Wind Speed' in ds.keys():
                        ds = ds.rename({'QC for Wind Speed': 'sfcWind_qc'})
                        ds['sfcWind_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'sfcWind_qc')
                        ds['sfcWind_qc'].attrs['flag_meanings'] =  "See QA/QC csv for network."                 
                        ds['sfcWind'].attrs['ancillary_variables'] = "sfcWind_qc" # List other variables associated with variable (QA/QC)
                    
                    ds['sfcWind'].attrs['comment'] = "Converted from mph to m/s."

                
                # sfcWind_dir: wind direction
                if "Wind Direction (0-360)" in ds.keys(): # No conversions needed, do not make raw column.
                    ds = ds.rename({'Wind Direction (0-360)': 'sfcWind_dir'})
                    ds['sfcWind_dir'].attrs['long_name'] = "wind_direction"
                    ds['sfcWind_dir'].attrs['standard_name'] = "wind_from_direction"
                    ds['sfcWind_dir'].attrs['units'] = "degrees_clockwise_from_north"
                    
                    if 'QC for Wind Direction' in ds.keys():
                        ds = ds.rename({'QC for Wind Direction': 'sfcWind_dir_qc'})
                        ds['sfcWind_dir_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'sfcWind_dir_qc')
                        ds['sfcWind_dir_qc'].attrs['flag_meanings'] = "See QA/QC csv for network."  
                        ds['sfcWind_dir'].attrs['ancillary_variables'] = "sfcWind_dir_qc" # List other variables associated with variable (QA/QC)

                    ds['sfcWind_dir'].attrs['comment'] = "Wind direction is defined by the direction that the wind is coming from (i.e., a northerly wind originates in the north and blows towards the south)."

                    
                # Other variables: rename to match format

                # Partial vapor pressure (kPa -> Pa) ## No CMIP standard name for this var, just CF.
                # This is calculated by CIMIS from relative humidity and air temperature data.
                if 'Vapor Pressure (kPa)' in ds.keys():
                    ds['pvp_derived'] = calc_clean._unit_pres_kpa_to_pa(ds['Vapor Pressure (kPa)'])
                    ds = ds.drop("Vapor Pressure (kPa)")

                    ds['pvp_derived'].attrs['long_name'] = "partial_vapor_pressure"
                    ds['pvp_derived'].attrs['standard_name'] = "water_vapor_partial_pressure_in_air"
                    ds['pvp_derived'].attrs['units'] = "Pa"
                    
                    if "QC for Vapor Pressure" in ds.keys():
                        ds = ds.rename({'QC for Vapor Pressure': 'pvp_derived_qc'})
                        ds['pvp_derived_qc'].attrs['flag_values'] = var_to_unique_list(ds, 'pvp_derived_qc')
                        ds['pvp_derived_qc'].attrs['flag_meanings'] = "See QA/QC csv for network." 
                        ds['pvp_derived'].attrs['ancillary_variables'] = "pvp_qc" # List other variables associated with variable (QA/QC)

                    ds['pvp_derived'].attrs['comment'] = "Derived by CIMIS from relative humidity and air temperature measurements. Converted from kPa to Pa."
                    
                
                #Testing: Manually check values to see that they seem correctly scaled, no unexpected NAs.
                # for var in ds.variables:
                #     try:
                #         print([var, float(ds[var].min()), float(ds[var].max())]) 
                #     except:
                #         next

                # Quality control: if any variable is completely empty, drop it.
                for key in ds.keys():
                    try:
                        if np.isnan(ds[key].values).all():
                            print("Dropping {}".format(key))
                            ds = ds.drop(key)
                    except: # Add to handle errors for unsupported data types
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
                print(e)
                errors['File'].append(station) # Save ID of station.
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue

            #Write station file to netcdf.
            if ds is None: # Should be caught in error handling above, but add in case.
                print("{} not saved.".format(file))
                errors['File'].append(file) # If stat_files is none, this will default to saving ID of station.
                errors['Time'].append(end_api)
                errors['Error'].append("File has no data.")
                continue
            
            else:
                try:
                    filename = station_id+".nc" # Make file name
                    filepath = cleandir+filename # Write file path
                    print(filepath) # For testing

                    # Write locally
                    ds.to_netcdf(path = 'temp/temp.nc', engine = 'h5netcdf') # Save station file.

                    # Push file to AWS with correct file name.
                    s3.Bucket(bucket_name).upload_file('temp/temp.nc', filepath)

                    print("Saving {} with dims {}".format(filename, ds.dims))
                    ds.close() # Close dataframe.
                    
                except Exception as e:
                    print(e)
                    errors['File'].append(station)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)
                    continue  

        # # Save the list of removed variables to AWS
        removedvars = pd.DataFrame(removecols, columns = ["Variable"])
        csv_buffer = StringIO()
        removedvars.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"removedvars.csv")
        
   
    # Write errors.csv
    finally: # Always execute this.
        print(errors)
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"errors_{}_{}.csv".format(network, end_api))

   
# # Run functions
if __name__ == "__main__":
    rawdir, cleandir, qaqcdir = get_file_paths("CIMIS")
    print(rawdir, cleandir, qaqcdir)
    clean_cimis(rawdir, cleandir)
    
    # # Testing:
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
    
    