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
# import warnings
# warnings.filterwarnings(action = 'ignore', category = FutureWarning) # Optional: Silence pandas' future warnings about regex (not relevant here)
import pandas as pd
import boto3
from io import BytesIO, StringIO
import gzip
import requests
from bs4 import BeautifulSoup
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


## FUNCTION: Given a network name, return all relevant AWS filepaths for other functions.
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir


## FUNCTION: Generate list of station and instrument elevations
# Input: url to tables on NDBC website
# Output: returns dataframe of stations and elevations. Called internally in clean_maritime() function
# This function generates a dataframe of station elevations and instrument elevations from the NDBC website, as this data is frequently missing from the data source
# This is joined by station to the xarray dataframes in the clean_maritime() function
def get_elevs(url):
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    tables = soup.find_all('pre')

    # Table 1 is smaller than table 2 and 3 by one column.
    # Start with table 1.
    tabletext = tables[0]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Barometer_Elevation"]
    table = tabletext.get_text().rsplit('ELEVATION',1)[1] # Remove headers.
    table = table.split() # Remove whitespace.
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+5] for x in range(0, len(table),5)] # Split into rows.
    df = pd.DataFrame(composite_list) # Turn to dataframe.
    df.columns = columns
    df['Tide_Reference'] = np.NAN # Add 6th column -- don't really  need this column, drop in update
#     df = df.reindex(columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference", "Barometer_Elevation"])

    # Table 2 has 6 columns
    tabletext = tables[1]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference", "Barometer_Elevation"]
    table = tabletext.get_text().rsplit('ELEVATION',1)[1] # Remove headers.
    table = table.split() # Remove whitespace.
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+6] for x in range(0, len(table),6)] # Split into rows.
    dftemp = pd.DataFrame(composite_list) # Turn to dataframe.
    dftemp.columns = columns
    df = pd.concat([df, dftemp])
    df=df.reset_index(drop=True)

#     Table 3 has 9 columns, but we only want the first 6.
    tabletext = tables[2]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference", "Barometer_Elevation"]
    table = tabletext.get_text().rsplit('CIRCLE',1)[1] # Remove headers.
    table = table.split() # Remove whitespace.
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+9] for x in range(0, len(table),9)] # Split into rows.
    dftemp = pd.DataFrame(composite_list) # Turn to dataframe.
    # Drop last three columns.
    dftemp = dftemp.iloc[:,0:6]
    dftemp.columns = columns
    df = pd.concat([df, dftemp])
    df=df.reset_index(drop=True)
    # print(df) # testing
    return df


## FUNCTION: Clean MARITIME and NDBC data
# Input:
# bucket_name: name of AWS bucket
# rawdir: path to where the raw data is saved as .txt.gz files, with each file either representing a full year of data, or monthly for current year
# cleadir: path to where the cleaned data should be saved
# network: network to be cleaned
# Output: Cleaned data for an individual network, priority variables, all times. Organized by staiton as .nc file
def clean_buoys(rawdir, cleandir, network):
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
        errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
        end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
        timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")
        print('successfully grabbed {} files'.format(len(files))) # testing

    except Exception as e: # If unable to read files from rawdir, break function
        print(e)
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    else: # If files read successfully, continue
        # for station in stations: # Full run
        for station in stations.sample(1): # SUBST FOR TESTING
            station_metadata = station_file.loc[station_file['STATION_ID']==float(station)] ## FLAGGING HERE - need to add capability for non-numeric station_ids for maritime
            station_id = network+"_"+str(station)
            print('Parsing: ', station_id) # testing

            for file in files:
                obj = s3.Object(bucket_name, file)
            try:
                with gzip.GzipFile(fileobj=obj.get()["Body"]) as gzipped_txt_file:
                    # Modified from: https://unidata.github.io/siphon/latest/_modules/siphon/simplewebservice/ndbc.html
                    col_names = ['year', 'month', 'day', 'hour', 'minute',
                                 'sfcWind_dir', 'sfcWind', 'sfcWind_gust',
                                 'wave_height', 'dominant_wave_period', 'average_wave_period', 'dominant_wave_dir',
                                 'ps', 'tas', 'water_temperature', 'tdps',
                                 'visibility', '3hr_pressure_tendency', 'water_level_above_mean']

                    col_units = {'sfcWind_dir': 'degrees',   # degT is "degrees true" for wind (might just be literally degrees, need to confirm)
                                 'sfcWind': 'm/s',
                                 'sfcWind_gust': 'm/s',
                                 'wave_height': 'm',
                                 'dominant_wave_period': 's',
                                 'average_wave_period': 's',
                                 'dominant_wave_direction': 'degrees',
                                 'ps': 'hPa',
                                 'tas': 'degC',
                                 'water_temperature': 'degC',
                                 'tdps': 'degC',
                                 'visibility': 'nautical_mile',
                                 '3hr_pressure_tendency': 'hPa',
                                 'water_level_above_mean': 'ft',
                                 'time': None}

                    df = pd.read_csv(gzipped_txt_file, comment='#', na_values='MM', names=col_names, sep='\s+')

                    # Convert time to datetime object
                    df['time'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']], utc=True)
                    df = df.drop(columns=['year', 'month', 'day', 'hour', 'minute'])
                    # df['units'] = col_units

                    # Standardize NA codes
                    try:
                        df = df.replace(999.0, np.nan) # sfcWind_dir, tas, tdps
                        df = df.replace(99.0, np.nan) # sfcWind
                        df = df.replace(9999.0, np.nan) # ps
                    except Exception as e:
                        print(e) # Check to see nan standardization worked as desired

                    # Drop unnecessary vars
                    cols_to_drop = ['sfcWind_gust','wave_height', 'dominant_wave_period', 'average_wave_period',
                                    'dominant_wave_dir', 'water_temperature', 'visibility', '3hr_pressure_tendency',
                                    'water_level_above_mean']
                    df = df.drop(columns = cols_to_drop)
                    print(df.head) # testing

                    # TO DO: Canadian zip files will need extra layer here
                    # Read the file as a zipfile and process the members
                    # with zipfile.ZipFile(tf, mode='r') as zipf: # Unzip ## USE FOR CANADIAN ZIP FILES
                    # https://stackoverflow.com/questions/18885175/read-a-zipped-file-as-a-pandas-dataframe

            # Handle exceptions thrown during individaul file read in
            except Exception as e:
                print(e)
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)

            # Next, start processing data columns
            if df.empty is True:
                print('df {} not saved'.format(file))

            if df.empty is False:
                try:
                    # TIME FILTER: Remove any rows before Jan 01 1980 and after August 30 2022
                    df = df.loc[(df['time']<'2022-09-01') & (df['time']>'1979-12-31')]

                    # Move df to xarray object
                    ds = df.to_xarray()

                    # Update global attributes
                    ds = ds.assign_attrs(title = network+" cleaned")
                    ds = ds.assign_attrs(institution = 'Eagle Rock Analytics / Cal Adapt')
                    ds = ds.assign_attrs(source = '')
                    ds = ds.assign_attrs(history = 'MARITIME_clean.py script run on {} UTC'.format(timestamp)) # note that the script is currently called NDBC_clean, but will rename once complete to match pull stage
                    ds = ds.assign_attrs(comment = 'Intermediate data product: may not have been subject to any cleaning or QA/QC processing.')
                    ds = ds.assign_attrs(license = '')
                    ds = ds.assign_attrs(citation = '')
                    ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). \
                    It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, \
                    nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; \
                    nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, \
                    nor has the Commission passed upon the accuracy of the information in this document.")
                    ds = ds.assign_attrs(station_name = station_id)
                    # ds = ds.assign_attrs(raw_files_merged = file_count) # Keep count of how many files merged per station.

                    # Add dimensions and coordinates
                    ds = ds.set_coords('time').swap_dims({'index': 'time'}) # Swap index with time
                    ds = ds.assign_coords(id = str(station_id))
                    ds = ds.expand_dims('id') # Add station_id as index
                    ds = ds.drop_vars(('index')) # Drop station_id variable and index coordinate
                    ds = ds.rename({'id': 'station'}) # Rename id to station_id

                    # Add coordinates: latitude and longitude
                    # Note: Datafiles do not have lat/lon information, grab from stationlist file
                    try:
                        if station_file[station_file['STATION_ID'].str.contains(station)].index.values.size > 0:
                            print('station is in list for lat/lon') # testing, can delete
                            idx = station_file[station_file['STATION_ID'].str.contains(station)].index.values
                            ds['latitude'] = float(station_file['LATITUDE'].iloc[idx])
                            ds['longitude'] = float(station_file['LONGITUDE'].iloc[idx])
                        else:
                            print('This station is not in the station list -- please check') # Can delete
                    except Exception as e:
                        print(e)
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append(e)

                    # Add variable: elevation
                    # Note: Datafriles do not have elevation information, grab from NDBC
                    url = 'https://www.ndbc.noaa.gov/bmanht.shtml'
                    elevs_df = get_elevs(url)

                    try:
                        if elevs_df[elevs_df['Station_ID'].str.contains(station)].index.values.size > 0:
                            print('station is in list for elevation') # testing
                            idx = elevs_df[elevs_df['Station_ID'].str.contains(station)].index.values

                            if (elevs_df['Site_Elevation'].iloc[idx].values) != 'NA':
                                print('elevation is not NA') # testing
                                ds['elevation'] = float(elevs_df['Site_Elevation'].iloc[idx])
                            else:
                                print('elevation is NA') # testing
                                ds['elevation'] = np.nan # Some stations have NA in for elevation
                        else:
                            print('This station is not in the elevation list -- setting to NaN') # Can delete
                            ds['elevation'] = np.nan # Handles if station is not in the list above

                    except Exception as e:
                        print(e)
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append(e)


                    # Update dimension and coordinate attributes
                    ds['time'] = pd.to_datetime(ds['time'].values, utc=True)
                    ds['time'] = pd.to_datetime(ds['time'].values, unit='ns')

                    # Update attributes
                    ds['time'].attrs['long_name'] = 'time'
                    ds['time'].attrs['standard_name'] = 'time'
                    ds['time'].attrs['comment'] = 'In UTC'

                    # Station ID
                    ds['station'].attrs['long_name'] = 'station_id'
                    ds['station'].attrs['comment'] = 'Unique ID created by Eagle Rock Analytics. Includes network name appended to original unique station ID provided by network.'

                    # Latitude
                    ds['lat'].attrs['long_name'] = "latitude"
                    ds['lat'].attrs['standard_name'] = "latitude"
                    ds['lat'].attrs['units'] = "degrees_north"

                    # Longitude
                    ds['lon'].attrs['long_name'] = "longitude"
                    ds['lon'].attrs['standard_name'] = "longitude"
                    ds['lon'].attrs['units'] = "degrees_east"

                    # Elevation
                    ds['elevation'].attrs['standard_name'] = 'height_above_mean_sea_level'
                    ds['elevation'].attrs['long_name'] = 'station_elevation'
                    ds['elevation'].attrs['units'] = 'meters'
                    ds['elevation'].attrs['positive'] = 'up' # Defines which direction is positive


                    # Next, update/add variable attributes and do unit conversions
                    # tas: surface air temperature (K)
                    if 'tas' in ds.keys():
                        ds['tas'] = calc_clean._unit_degC_to_K(ds['tas'])
                        ds['tas'].attrs['long_name'] = 'air_temperature'
                        ds['tas'].attrs['standard_name'] = 'air_temperature'
                        ds['tas'].attrs['units'] = 'degree_Kelvin'
                        ds['tas'].attrs['comment'] = 'Converted from degC to K'

                    # ps: surface air pressure (Pa)
                    if 'ps' in ds.keys():
                        ds['ps'] = calc_clean._unit_pres_hpa_to_pa(ds['ps'])
                        ds['ps'].attrs['long_name'] = 'station_air_pressure'
                        ds['ps'].attrs['standard_name'] = 'air_pressure'
                        ds['ps'].attrs['units'] = 'Pa'
                        ds['ps'].attrs['comment'] = 'Converted from hPa to Pa'

                    # tdps: dew point temperature (K)
                    if 'tdps' in ds.keys():
                        ds['tdps'] = calc_clean._unit_degC_to_K(ds['tdps'])
                        ds['tdps'].attrs['long_name'] = 'dew_point_temperature'
                        ds['tdps'].attrs['standard_name'] = 'dew_point_temperature'
                        ds['tdps'].attrs['units'] = 'degree_Kelvin'
                        ds['tdps'].attrs['comment'] = 'Converted from degC to K'

                    # pr: precipitation
                    # Note: not measured by NDBC or MARITIME

                    # hurs: relative humidity (%)
                    # Note: not measured by NDBC or MARITIME
                    # Note: will be calcualted with tas and tdps

                    # rsds: solar radiation (W/m2)
                    # Note: not measured by NDBC or MARITIME

                    # sfcWind: surface wind speed (m/s)
                    if 'sfcWind' in ds.keys():
                        ds['sfcWind'].attrs['long_name'] = 'wind_speed'
                        ds['sfcWind'].attrs['standard_name'] = 'wind_speed'
                        ds['sfcWind'].attrs['units'] = 'm s-1'

                    # sfcWind_dir: wind direction (degrees)
                    if 'sfcWind_dir' in ds.keys():
                        ds['sfcWind_dir'].attrs['long_name'] = 'wind_direction'
                        ds['sfcWind_dir'].attrs['standard_name'] = 'wind_from_direction'
                        ds['sfcWind_dir'].attrs['units'] = 'degrees_clockwise_from_north'


                    # Need to consider this - many buoys will not have any of the required data
                    # Quality control: if any variable is completely empty, drop it.
                    # for key in ds.keys():
                    #     try:
                    #         if np.isnan(ds[key].values).all():
                    #             print("Dropping {}".format(key))
                    #             ds = ds.drop(key)
                    #     except: # Add to handle errors for unsupported data types
                    #         next

                    # Need to consider what to do if the entire dataset will be empty because no data is observed

                    # Reorder variables
                    desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
                    rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom
                    new_index = desired_order = rest_of_vars
                    ds = ds[new_index]

                    print(ds) # Testing

                except Exception as e:
                    print(e)
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)

    #             # Write station file to netcdf format
    #             if ds is None:
    #                 print("ds {} not saved.".format(file))
    #                 continue
    #             else:
    #                 try:
    #                     filename = station_id + ".nc" # Make file name
    #                     filepath = cleandir + filename # Writes file path
    #                     print(filepath) # for testing
    #
    #                     ds.to_netcdf(path = filepath)
    #
    #                     s3.Bucket(bucket_name).upload_file(filepath)
    #
    #                     print('Saving {} with dims {}'.format(filename, ds.dims))
    #                     ds.close() # Close dataframe
    #
    #                 except Exception as e:
    #                     print(filename, e)
    #                     errors['File'].append(filename)
    #                     errors['Time'].append(end_api)
    #                     errors['Error'].append(e)
    #                     continue
    #
    # # Write errors to csv
    # finally:
    #     print(errors) # Testing
    #     errors = pd.DataFrame(errors)
    #     csv_buffer = StringIO()
    #     errors.to_csv(csv_buffer)
    #     content = csv_buffer.getvalue()
    #     s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"errors_{}_{}.csv".format(network, end_api)) # Make sure error files save to correct directory


# Run functions
if __name__ == "__main__":
    network = "NDBC" # or "MARITIME"
    rawdir, cleandir, qaqcdir = get_file_paths(network)
    print(rawdir, cleandir, qaqcdir) # TESTING
    clean_buoys(rawdir, cleandir, network=network)


# To do upon completion:
# 1. Run through maritime_clean.py to ensure that missing pieces (elevation) are correctly incorporated and delete when complete
# 2. Rename to maritime_clean.py so to be consistent with maritime_pull.py
# 3. go through imports to remove ones not used
# 4. drop tide_reference column in get_elevs
