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
from datetime import datetime
import numpy as np
import pandas as pd
import boto3
from io import BytesIO, StringIO
import gzip
import requests
from bs4 import BeautifulSoup
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


## FUNCTION: Given a network name, return all relevant AWS filepaths for other functions.
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir


## FUNCTION: Generate list of station and instrument elevations
# Input: url to tables on NDBC website
# Output: returns dataframe of stations and elevations. Called internally in clean_buoys() function
# This function generates a dataframe of station elevations and instrument elevations from the NDBC website, as this data is frequently missing from the data source
# This is joined by station to the xarray dataframes in the clean_buoys() function
def get_elevs(url):
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    tables = soup.find_all('pre')

    # Table 1 is smaller than table 2 and 3 by one column
    # Start with table 1
    tabletext = tables[0]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Barometer_Elevation"]
    table = tabletext.get_text().rsplit('ELEVATION',1)[1] # Remove headers
    table = table.split() # Remove whitespace
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+5] for x in range(0, len(table),5)] # Split into rows
    df = pd.DataFrame(composite_list)
    df.columns = columns
    df['Tide_Reference'] = np.NAN # Add 6th column -- don't really  need this column, drop in update
#     df = df.reindex(columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference", "Barometer_Elevation"])

    # Table 2 has 6 columns
    tabletext = tables[1]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference", "Barometer_Elevation"]
    table = tabletext.get_text().rsplit('ELEVATION',1)[1] # Remove headers
    table = table.split() # Remove whitespace
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+6] for x in range(0, len(table),6)] # Split into rows
    dftemp = pd.DataFrame(composite_list)
    dftemp.columns = columns
    df = pd.concat([df, dftemp])
    df=df.reset_index(drop=True)

#   Table 3 has 9 columns, but we only want the first 6
    tabletext = tables[2]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference", "Barometer_Elevation"]
    table = tabletext.get_text().rsplit('CIRCLE',1)[1] # Remove headers
    table = table.split() # Remove whitespace
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+9] for x in range(0, len(table),9)] # Split into rows
    dftemp = pd.DataFrame(composite_list)
    dftemp = dftemp.iloc[:,0:6] # Drop last three columns
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
    emptybuoy = [] # intialize list of buoys that do not have wx data for various reasons
    try:
        files = [] # Get files
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = rawdir):
            file = str(item.key)
            files += [file]

        # Get station file and read in metadata
        station_file = [file for file in files if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj['Body'].read()))
        stations = station_file['STATION_ID'].dropna()

        # Remove error, station files
        # files = [file for file in files if '.txt.gz' in file] # US-owned stations
        # files = [file for file in files if '.zip' in file] # Canada-owned stations
        files = [file for file in files if 'stationlist' not in file]
        files = [file for file in files if 'error' not in file]

        # Set up error handling.
        errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
        end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
        timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    except Exception as e: # If unable to read files from rawdir, break function
        print(e)
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    else: # If files read successfully, continue
        # for station in stations: # Full run
        # for station in stations.sample(4): # SUBSET FOR TESTING
        for station in ['46138']:
            station_id = network+"_"+str(station)
            print('Parsing: ', station_id)

            df_stat = None # Initialize merged df
            try:
                stat_files = [k for k in files if station in k] # Gets list of files from the same station

                if not stat_files: # If station has no files downloaded
                    emptybuoy += [station_id]
                    print('{} does not have any raw data downloaded to AWS at present. Several buoys do not record any meteorological observations.'.format(station_id))
                    errors['File'].append(station_id)
                    errors['Time'].append(end_api)
                    errors['Error'].append('No raw data found for this station.')
                    othercols = None # setting to none here so script doesn't fail if the first station grabbed is one with no data
                    continue # Skip this station

                for file in stat_files:
                    try:
                        # handling for different file extensions (.txt.gz = US-owner, .zip = Canada-owner)
                        if file.endswith('.txt.gz'): # US-based owner
                            obj = s3.Object(bucket_name, file)
                            with gzip.GzipFile(fileobj=obj.get()["Body"]) as file_to_parse:
                                # modified from: https://unidata.github.io/siphon/latest/_modules/siphon/simplewebservice/ndbc.html
                                df = pd.read_csv(file_to_parse, sep='\s+', low_memory=False)

                                # older files are missing the minute column
                                if {'mm'}.issubset(df.columns) == False:
                                    df.insert(loc=4, column='mm', value='00') # manually setting to top of hour, our process will collapse all other obs to top of hour at next stage

                                # fix year label mismatch
                                yr_raw = str(df.iloc[1][0]).split('.')[0]
                                if len(yr_raw) != 4: # older files have a two-digit year
                                    if df.iloc[1][0] >= 80: # 1980-1999
                                        df['YYYY'] = df['YY'].apply(lambda x: "{}{}".format('19', x))
                                        df = df.iloc[:, 1:]

                                    elif df.iloc[1][0][0] <= 23: # 2000-present
                                        df['YYYY'] = df['YY'].apply(lambda x: "{}{}".format('20', x))
                                        df = df.iloc[:, 1:]

                                if (df.columns[0][0].isdigit()) == False: # newer files have a 2-line header with comments
                                    df.drop([0], inplace=True)
                                    df.rename(columns={'#YY':'YYYY'}, inplace=True)

                                # fix variable label mismatch
                                if {'WD', 'BAR'}.issubset(df.columns) == True: # older files have different var names
                                    df.rename(columns={'WD':'WDIR', 'BAR':'PRES'}, inplace=True)

                                df['WDIR'] = pd.to_numeric(df['WDIR'])
                                df['WSPD'] = pd.to_numeric(df['WSPD'])
                                df['PRES'] = pd.to_numeric(df['PRES'])
                                df['ATMP'] = pd.to_numeric(df['ATMP'])
                                df['DEWP'] = pd.to_numeric(df['DEWP'])

                                # standardize NA codes -- CHECK ON THIS
                                try:
                                    df.replace(999, np.nan, inplace=True) # sfcWind_dir
                                    df.replace(999.0, np.nan, inplace=True) # sfcWind_dir, tas, tdps
                                    df.replace(99.0, np.nan, inplace=True) # sfcWind
                                    df.replace(9999.0, np.nan, inplace=True) # ps

                                    # shouldn't have to do this twice, but not working without it
                                    df.replace('999', np.nan, inplace=True) # sfcWind_dir
                                    df.replace('999.0', np.nan, inplace=True) # sfcWind_dir, tas, tdps
                                    df.replace('99.0', np.nan, inplace=True) # sfcWind
                                    df.replace('9999.0', np.nan, inplace=True) # ps

                                except Exception as e:
                                    print(e)

                                # convert date to datetime
                                df.rename(columns={'YYYY':'year', 'MM':'month', 'DD':'day', 'hh':'hour', 'mm':'minute'}, inplace=True)
                                df['time'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']], utc=True)
                                df = df.drop(columns=['year', 'month', 'day', 'hour', 'minute'])

                                # time filter: remove any rows before Jan 01 1980 and after August 30 2022
                                df = df.loc[(df['time']<'2022-09-01') & (df['time']>'1979-12-31')]

                                # drop variables if not desired variable
                                cols_to_keep = ['WDIR', 'WSPD', 'PRES', 'ATMP', 'DEWP', 'time']
                                othercols = [col for col in df.columns if col not in cols_to_keep]
                                df = df[df.columns.intersection(cols_to_keep)] # drop all columns not in cols_to_keep list

                                df.rename(columns={'WDIR':'sfcWind_dir',
                                                  'WSPD':'sfcWind',
                                                  'PRES':'ps',
                                                  'ATMP':'tas',
                                                  'DEWP':'tdps'}, inplace=True)

                                # if more than one file per station, merge files together
                                if df_stat is None:
                                    df_stat = df
                                    del(df) # deleting for memory
                                else:
                                    if len(stat_files) > 1: # if there is more than one file per station
                                        df_stat = pd.concat([df_stat, df], axis=0, ignore_index=True)

                        elif file.endswith('.zip'): # canadian buoy
                            obj = s3.Bucket(bucket_name).Object(file)

                            with BytesIO(obj.get()["Body"].read()) as tf:
                                tf.seek(0) # rewind the file

                                with zipfile.ZipFile(tf, mode='r') as zipf:

                                    for subfile in zipf.namelist():
                                        filestn = subfile.replace(".csv", "")
                                        filestn = int(filestn[1:]) # drops proceeding 'c'

                                        if int(station) == filestn:
                                            df = pd.read_csv(zipf.open(subfile), low_memory=False)
                                            # https://www.meds-sdmm.dfo-mpo.gc.ca/isdm-gdsi/waves-vagues/formats-eng.html

                                            # convert date to datetime
                                            df['time'] =  pd.to_datetime(df['DATE'], format='%m/%d/%Y %H:%M', utc=True)

                                            # time filter: remove any rows before Jan 01 1980 and after August 30 2022
                                            df = df.loc[(df['time']<'2022-09-01') & (df['time']>'1979-12-31')]

                                            # drop variables if not desired variable
                                            # cols_to_keep = ['time', 'LATITUDE', 'LONGITUDE', 'WDIR', 'WDIR.1', 'WSPD.1', 'WSPD', 'ATMS', 'ATMS.1', 'DRYT'] # waiting on response about ".1" - thinking these are secondary sensors, but with different nan codes
                                            cols_to_keep = ['time', 'Q_FLAG', 'LATITUDE', 'LONGITUDE', 'WDIR', 'WSPD', 'ATMS', 'DRYT'] # 'SLEV' is sea level height, could also be useful
                                            othercols = [col for col in df.columns if col not in cols_to_keep]
                                            df = df[df.columns.intersection(cols_to_keep)] # drop all columns not in cols_to_keep list

                                            df.rename(columns={'LATITUDE':'lat',
                                                                'LONGITUDE':'lon',
                                                                'WDIR':'sfcWind_dir',
                                                                'WSPD':'sfcWind',
                                                                'ATMS':'ps',
                                                                'DRYT':'tas',
                                                                'Q_FLAG':'qc_flag',
                                                                'time':'time'}, inplace=True)

                                            # missing data flags - mainly as a catchall in case qc check did not catch
                                            df.replace(999, np.nan, inplace=True) # sfcWind_dir
                                            df.replace(999.0, np.nan, inplace=True) # sfcWind_dir, tas, tdps
                                            df.replace(99.0, np.nan, inplace=True) # sfcWind
                                            df.replace(9999.0, np.nan, inplace=True) # ps

                                            # if more than one file per station, merge files together
                                            if df_stat is None:
                                                df_stat = df
                                                del(df) # deleting for memory
                                            else:
                                                if len(stat_files) > 1: # if there is more than one file per station
                                                    df_stat = pd.concat([df_stat, df], axis=0, ignore_index=True)

                    except Exception as e:
                        print(e)
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append(e)
                        continue

                # Format joined station file
                file_count = len(stat_files)

                # Sort by time and remove any overlapping timesteps
                df_stat = df_stat.sort_values(by='time')
                df_stat = df_stat.drop_duplicates()

                # Move df to xarray object
                ds = df_stat.to_xarray()
                del(df_stat)

                # Update global attributes
                ds = ds.assign_attrs(title = network+" cleaned")
                ds = ds.assign_attrs(institution = 'Eagle Rock Analytics / Cal Adapt')
                ds = ds.assign_attrs(source = '')
                ds = ds.assign_attrs(history = 'MARITIME_clean.py script run on {} UTC'.format(timestamp))
                ds = ds.assign_attrs(comment = 'Intermediate data product: may not have been subject to any cleaning or QA/QC processing.')
                ds = ds.assign_attrs(license = '')
                ds = ds.assign_attrs(citation = '')
                ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")
                # ds = ds.assign_attrs(station_name = station_file['NAME'].values[0]) ## FIX
                ds = ds.assign_attrs(raw_files_merged = file_count) # Keep count of how many files merged per station

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
                        # print('station is in list for lat/lon') # testing, can delete
                        idx = station_file[station_file['STATION_ID'].str.contains(station)].index.values
                        lat = np.asarray([float(station_file['LATITUDE'].iloc[idx])] * len(ds['time']))
                        lat.shape = (1, len(ds['time']))
                        lon = np.asarray([float(station_file['LONGITUDE'].iloc[idx])] * len(ds['time']))
                        lon.shape = (1, len(ds['time']))
                        ds = ds.assign_coords(lat = (["station","time"],lat), lon = (["station","time"], lon))
                    else:
                        print('This station is not in the station list -- please check') # shouldn't ever be the case, can delete, but good for an initial flag
                except Exception as e:
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)

                # Add variable: elevation
                # Note: Datafiles do not have elevation information, grab from NDBC
                # Note: We are keeping the valid NA elevations for stations at this stage. Will infill from DEM in Stage 3: QA/QC
                url = 'https://www.ndbc.noaa.gov/bmanht.shtml'
                elevs_df = get_elevs(url)
                try:
                    if elevs_df[elevs_df['Station_ID'].str.contains(station)].index.values.size > 0: # station is in NDBC elevation list
                        idx = elevs_df[elevs_df['Station_ID'].str.contains(station)].index.values
                        if (elevs_df['Site_Elevation'].iloc[idx].values) != 'NA': # elevation has a non-NA valid value
                            elev = np.asarray([float(elevs_df['Site_Elevation'].iloc[idx])] * len(ds['time']))
                            elev.shape = (1, len(ds['time']))
                        else:
                            elev = np.asarray([np.nan] * len(ds['time'])) # some stations have NA in for elevation
                            elev.shape = (1, len(ds['time']))
                    else:
                        print('This station is not in the elevation list -- setting to NaN') # station in our stationlist but not in the NDBC elevation list
                        elev = np.asarray([np.nan] * len(ds['time']))
                        elev.shape = (1, len(ds['time']))
                    ds['elevation'] = (['station', 'time'], elev)
                except Exception as e:
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)

                # Add sensor heights, grab from NDBC
                # defaulting to nan, as we are assuming most stations will only have some of these sensor heights available
                ds = ds.assign_attrs(thermometer_height_m = np.nan)
                ds = ds.assign_attrs(barometer_elevation_m = np.nan)
                ds = ds.assign_attrs(anemometer_height_m = np.nan)
                try:
                    if elevs_df[elevs_df['Station_ID'].str.contains(station)].index.values.size > 0: # station exists in NDBC list
                        idx = elevs_df[elevs_df['Station_ID'].str.contains(station)].index.values

                        # air temperature
                        if (elevs_df['Air_Temp_Elevation'].iloc[idx].values) != 'NA': # air temp sensor has value
                            ds.attrs['thermometer_height_m'] = float(elevs_df['Air_Temp_Elevation'].iloc[idx])

                        # air pressure
                        if (elevs_df['Barometer_Elevation'].iloc[idx].values) != 'NA': # air pressure sensor has value
                            ds.attrs['barometer_elevation_m'] = float(elevs_df['Barometer_Elevation'].iloc[idx])

                        # wind
                        if (elevs_df['Anemometer_Elevation'].iloc[idx].values) != 'NA': # wind sensor has value
                            ds.attrs['anemometer_height_m'] = float(elevs_df['Anemometer_Elevation'].iloc[idx])
                    else:
                        # station does not exist in NDBC sensor height list but is in our station_list
                        print('This station is not in the sensor height list -- setting to NaN')
                except Exception as e:
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
                # Note: will be calculated with tas and tdps

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

                # qc_flag: quality control flag
                # note this flag appears that it applies to every observation at that time stamp - waiting on confirmation
                if "qc_flag" in ds.keys():
                    ds['qc_flag'].attrs['flag_values'] = "0 1 3 4 5 6 7 8 9"
                    ds['qc_flag'].attrs['flag_meanings'] =  "https://www.meds-sdmm.dfo-mpo.gc.ca/isdm-gdsi/waves-vagues/formats-eng.html"

                # drop any column that does not have any valid (non-nan data)
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
                    except: # Add to handle errors for unsupported data types
                        next

                # removes elevation and qc_flag (canadian buoy only) if the only remaining variables (occurs at least once)
                for key in ds.keys():
                    if 'qc_flag' in list(ds.keys()):
                        if len(ds.keys())==2:
                            print("Dropping empty var: {}".format('elevation'))
                            print("Dropping empty var: {}".format('qc_flag'))
                            ds = ds.drop('qc_flag')
                            ds = ds.drop('elevation')
                            continue

                # Reorder variables, in the following order
                desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
                desired_order = [i for i in desired_order if i in list(ds.keys())] # only keep vars that are in ds
                rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom (elevation)
                new_index = desired_order + rest_of_vars
                ds = ds[new_index]

            except Exception as e:
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue

            # Write station file to netcdf format
            if len(ds.keys())==0:   # skip station if the entire dataset will be empty because no data is observed (as in only ocean obs are recorded, but not needed)
                emptybuoy += [station_id]
                print("{} has no data for all meteorological variables of interest throughout its current reporting; station not cleaned.".format(station_id))
                errors['File'].append(station_id)
                errors['Time'].append(end_api)
                errors['Error'].append('Station has no data.')
                continue
            else:
                try:
                    filename = station_id + ".nc" # Make file name
                    filepath = cleandir + filename # Writes file path

                    # Write locally
                    ds.to_netcdf(path = 'temp/temp.nc', engine = 'h5netcdf') # Save station file.

                    # Push file to AWS with correct file name
                    s3.Bucket(bucket_name).upload_file('temp/temp.nc', filepath)

                    print('Saving {} with dims {}'.format(filename, ds.dims))
                    print(ds) # Testing, but nice to have
                    ds.close() # Close dataframe

                except Exception as e:
                    print(filename, e)
                    errors['File'].append(filename)
                    errors['Time'].append(end_api)
                    errors['Error'].append(e)
                    continue

        # Save list of removed variables to AWS
        removedvars = pd.DataFrame(othercols, columns = ['Variable'])
        csv_buffer = StringIO()
        removedvars.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"removedvars.csv")

    # Write errors to csv
    finally:
        print(errors) # Testing
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"errors_{}_{}.csv".format(network, end_api)) # Make sure error files save to correct directory

    # Save list of all buoy stations that do not report any wx data, but have raw data files on AWS
    print('Empty buoys to remove: ', emptybuoy)
    no_wx_buoys = pd.DataFrame(emptybuoy, columns = ['Station_ID'])
    csv_buffer = StringIO()
    no_wx_buoys.to_csv(csv_buffer)
    content = csv_buffer.getvalue()
    s3_cl.put_object(Bucket=bucket_name, Body=content, Key=cleandir+"emptybuoys.csv")

# Run functions
if __name__ == "__main__":
    network = "NDBC" # or "MARITIME"
    rawdir, cleandir, qaqcdir = get_file_paths(network)
    print(rawdir, cleandir, qaqcdir) # TESTING
    clean_buoys(rawdir, cleandir, network=network)
