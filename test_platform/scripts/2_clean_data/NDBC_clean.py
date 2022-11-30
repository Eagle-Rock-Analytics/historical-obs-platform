"""
This script performs data cleaning for NDBC and MARITIME data pulled from NCEI for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables and drop unnecessary variables
(2) Converts station metadata to standard format, with unique identifier
(3) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(4) Converts missing data to standard format
(5) Tracks existing qa/qc flag for review
(6) Merge files by station, and outputs cleaned variables as a single .nc file for each station in an individual network.
Inputs: Raw data for the network's stations, with each .dat.gz file representing a single day's data for a specified year.
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime, date, timedelta
import re
import numpy as np
#import warnings
#warnings.filterwarnings(action = 'ignore', category = FutureWarning) # Optional: Silence pandas' future warnings about regex (not relevant here)
import pandas as pd
import boto3
from io import BytesIO, StringIO
import random
# import gzip
import zipfile
# To be able to open xarray files from S3, h5netcdf must also be installed, but doesn't need to be imported.


## Import cleaning stage calc functions for conversions
try:
    import calc_clean
except:
    print("Error importing calc_clean.py")
    pass

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

## Set up directory to save files temporarily, if it doesn't already exist.
try:
    os.mkdir('temp') # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass


## Given a network name, return all relevant AWS filepaths for other functions.
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir


## Clean MARITIME and NDBC data
# Input:
# bucket_name: name of AWS bucket
# rawdir: path to where the raw data is saved as .txt.gz files, with each file either representing a full year of data, or monthly for current year
# cleadir: path to where the cleaned data should be saved
# network: network to be cleaned
# Output: Cleaned data for an individual network, priority variables, all times. Organized by staiton as .nc file
def clean_buoys(rawdir, cleandir, network):
    # network = "NDBC"
    try:
        # Get files
        files = []
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir):
            file = str(item.key)
            files += [file]

        # Get station file and read in metadata
        station_file = [file for file in files if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj['Body'].read()))
        stations = station_file['STATION_ID'].dropna() # handful of stations have letters in place of numbers

        # Remove error, station files
        files = [file for file in files if '.txt.gz' in file]
        # files = [file for file in files if '.zip' in file] ## COME BACK TO ZIP FILES - THESE ARE THE CANADIAN BUOYS
        files = [file for file in files if 'stationlist' not in file]
        files = [file for file in files if 'error' not in file]

        # Set up error handling.
        errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
        end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv.
        timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S") # For attributes of netCDF file.

        # Set up list of variables to be removed
        removecols = ['GST', 'WVHT', 'DPD', 'APD', 'MWD', 'WTMP', 'VIS', 'TIDE']

    except Exception as e: # If unable to read files from rawdir, break function.
        print(e)
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    else: # If files read successfully, continue
        # for station in stations: # Full run
        for station in stations.sample(2): # SUBST FOR TESTING
            station_metadata = station_file.loc[station_file['STATION_ID']==float(station)] ## FLAGGING HERE - need to add capability for non-numeric station_ids for maritime
            station_id = network+"_"+str(station)
            print(station_id) # testing

            dfs = []
            for file in files: # For each zip file (annual or monthly)
                try:
                    fileyear = file.split(".txt.gz")[0]
                    fileyear = fileyear[-4:]

                    obj = s3.Bucket(bucket_name).Object(file)
                    with BytesIO(obj.get()["Body"].read()) as tf:
                        # TO DO: Canadian zip files will need extra layer here
                        # Read the file as a zipfile and process the members
                        # with zipfile.ZipFile(tf, mode='r') as zipf: # Unzip ## USE FOR CANADIAN ZIP FILES

                        allcols = ['#YY', 'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST', 'WVHT', 'DPD', 'MWD',
                                    'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']

                        unitstocheck = ['degT', 'm/s', 'm', 'sec', 'hPa', 'degC', 'mi', 'ft']

                        df = pd.read_csv(tf, names=allcols, low_memory=False)
                        df = df.drop(columns=removecols, axis=1) # Drop columns

                        # Fix time into datetime
                        df['time'] = pd.to_datetime(df[['#YY', 'MM', 'DD', 'hh', 'mm']]).dt.strftime('%Y-%m-%d%H:%M') # Assuming UTC time?

                        # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022
                        df = df.loc[(df['time']<'2022-09-01') & (df['time']>'1979-12-31')]

                        # Remove unnecssary time columns
                        df = df.drop(['#YY', 'MM', 'DD', 'hh', 'mm'], axis=1)

                        dfs.append(df)

                # Handle exceptions thrown during individaul file read in
                except Exception as e:
                    print(e)
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)


                # Next, start processing data columns
                # try:
                #     file_count = len(dfs)
                #     print(file_count) # TESTING
                #     df_stat = pd.concat(dfs)
                #
                #     # Replace non-standard NAs -- CHECK THAT THESE ARE NOT EXCLUDING GOOD DATA
                #     df_stat = df_stat.replace(999.0, np.nan) # Air temp, dewpoint temp
                #     df_stat = df_stat.replace(9999.0, np.nan) # Presure
                #     df_stat = df_stat.replace(999, np.nan) # Wind direction
                #     df_stat = df_stat.replace(99.0, np.nan) # Wind speed
                #
                #     # Duplicated/overlapping timesteps?
                #
                #     # Move df to xarray object
                #     ds = df_stat.to_xarray()
                #     del(df_stat)
                #
                #     # Update global attributes
                #     ds = ds.assign_attrs(title = network+" cleaned")
                #     ds = ds.assign_attrs(institution = 'Eagle Rock Analytics / Cal Adapt')
                #     ds = ds.assign_attrs(source = '')
                #     ds = ds.assign_attrs(history = 'NDBC_clean.py script run on {} UTC'.format(timestamp))
                #     ds = ds.assign_attrs(comment = 'Intermediate data product: may not have been subject to any cleaning or QA/QC processing.')
                #     ds = ds.assign_attrs(license = '')
                #     ds = ds.assign_attrs(citation = '')
                #     ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). \
                #     It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, \
                #     nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; \
                #     nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, \
                #     nor has the Commission passed upon the accuracy of the information in this document.")
                #     ds = ds.assign_attrs(station_name = station_name)
                #
                #     # Add dimensions and coordinates
                #     ds = ds.set_coords('time').swap_dims({'index': 'time'}) # Swap index with time
                #     ds = ds.assign_coords(id = str(station_id))
                #     ds = ds.expand_dims('id') # Add station_id as index
                #     ds = ds.drop_vars(('index')) # Drop station_id variable and index coordinate
                #     ds = ds.rename({'id': 'station'}) # Rename id to station_id
                #
                #     # Add coordinates: latitude and longitude
                #     ## TO DO
                #
                #     # Add variable: elevation
                #     ## TO DO
                #     ## Look at get_elev in MARITIME_clean.py
                #
                #     # Update dimension and coordinate attributes
                #     ds['time'] = pd.to_datetime(ds['time'].values, utc=True)
                #     ds['time'] = pd.to_datetime(ds['time'].values, unit='ns')
                #
                #     # Update attributes
                #     ds['time'].attrs['long_name'] = 'time'
                #     ds['time'].attrs['standard_name'] = 'time'
                #     ds['time'].attrs['comment'] = 'In UTC.'
                #
                #     # Station ID
                #     ds['station'].attrs['long_name'] = 'station_id'
                #     ds['station'].attrs['comment'] = 'Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network.'
                #
                #
                #     # # Latitude
                #     # ds['lat'].attrs['long_name'] = "latitude"
                #     # ds['lat'].attrs['standard_name'] = "latitude"
                #     # ds['lat'].attrs['units'] = "degrees_north"
                #     #
                #     # # Longitude
                #     # ds['lon'].attrs['long_name'] = "longitude"
                #     # ds['lon'].attrs['standard_name'] = "longitude"
                #     # ds['lon'].attrs['units'] = "degrees_east"
                #     #
                #     # # Elevation
                #     # ## TO DO
                #     # ds['elevation_raw'] = ds['elevation']
                #     # ds['elevation'] = calc_clean._unit_elev_ft_to_m(ds['elevation']) # Converts from feet to meters
                #     #
                #     # ds['elevation'].attrs['standard_name'] = 'height_above_mean_sea_level'
                #     # ds['elevation'].attrs['long_name'] = 'station_elevation'
                #     # ds['elevation'].attrs['units'] = 'meters'
                #     # ds['elevation'].attrs['positive'] = 'up' # Defines which direction is positive
                #     # ds['elevation'].attrs['comment'] = 'Converted from feet.'
                #     #
                #     # ds['elevation_raw'].attrs['standard_name'] = "height_above_mean_sea_level"
                #     # ds['elevation_raw'].attrs['long_name'] = "station_elevation"
                #     # ds['elevation_raw'].attrs['units'] = "feet"
                #     # ds['elevation_raw'].attrs['positive'] = "up" # Define which direction is positive
                #
                #
                #     # Next, update/add variable attributes and do unit conversions
                #     # tas: surface air temperature (K)
                #
                #     # ps: surface air pressure (Pa)
                #
                #     # tdps: dew point temperature (K)
                #
                #     # pr: precipitation
                #     # Note: not measured by NDBC or MARITIME
                #
                #     # hurs: relative humidity (%)
                #     # Note: not measured by NDBC or MARITIME
                #
                #     # rsds: solar radiation (W/m2)
                #     # Note: not measured by NDBC or MARITIME
                #
                #     # sfcWind: surface wind speed (m/s)
                #
                #     # sfcWind_dir: wind direction (degrees)
                #
                #
                #     # Reorder variables
                #     # desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
                #     # rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom
                #     # new_index = desired_order = rest_of_vars
                #     # ds = ds[new_index]
                #     #
                #     #
                #     # # Merge file to previous time period data for same station
                #     # if ds_stat is None:
                #     #     ds_stat = ds # For first file in station, set ds_stat to be the original ds
                #     #
                #     # else:
                #     #     ds_stat = xr.merge([ds_stat, ds], compat='no-conflicts') # Otherwise, merge new records to first ds
                #     #     print('Merging records for station {}'.fofrmat(i))
                #     #     file_count +=1
                #     #     ds_stat.attrs['raw_files_merged'] = file_count
                #
                #
                # except Exception as e:
                #     print(e)
                #     errors['File'].append(file)
                #     errors['Time'].append(end_api)
                #     errors['Error'].append(e)

                # # Write station file to netcdf format
                # if ds_stat is None:
                #     print("{} not saved.".format(file))
                #     continue
                # else:
                #     try:
                #         filename = station_id + ".nc" # Make file name
                #         filepath = cleandir + filename # Writes file path
                #
                #         ds_stat.to_netcdf(path = filepath)
                #         print('Saving {} with dims {}'.format(filename, ds_stat.dims))
                #         ds_stat.close() # Close dataframe
                #     except Exception as e:
                #         print(filename, e)
                #         errors['File'].append(filename)
                #         errors['Time'].append(end_api)
                #         errors['Error'].append(e)
                #         continue

# # Run functions
if __name__ == "__main__":
    network = "NDBC" # or "MARITIME"
    rawdir, cleandir, qaqcdir = get_file_paths(network)
    print(rawdir, cleandir, qaqcdir) # TESTING
    clean_buoys(rawdir, cleandir, network=network)


# To do upon completion:
# 1. Run through maritime_clean.py to ensure that missing pieces (elevation) are correctly incorporated
