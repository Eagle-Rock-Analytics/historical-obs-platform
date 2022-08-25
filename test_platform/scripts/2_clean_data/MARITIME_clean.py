"""
This script is a template structure for data cleaning for the MARITIME network for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables, and calculates derived priority variables if not observed
(2) Drops unnecessary variables
(3) Converts station metadata to standard format, with unique identifier
(4) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(5) Converts missing data to standard format
(6) Tracks existing qa/qc flag for review
(7) Merge files by station, and outputs cleaned variables as a single .nc file for an individual network.
Inputs: Raw data for MARITIME stations, with each file representing a month of a year.
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
"""
# Small errors still coming up in error file, to be fixed with more time/capacity.
# Resolve HCF errors
# Station 46089: throws error, can't convert "NA" string into float. Some data (likely lat/lon) is tripping this.

# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime
import re
import geopandas as gp
import numpy as np
from calc import _calc_relhumid, get_wecc_poly, _calc_ps
import csv
import itertools
import pandas as pd
import requests
from bs4 import BeautifulSoup

# Set envr variables and calculate any needed variables
homedir = os.getcwd()
workdir = "test_platform/data/1_raw_wx/MARITIME/"
savedir = "test_platform/data/2_clean_wx/MARITIME/"
wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 

# Set up directory to save files, if it doesn't already exist.
try:
    os.mkdir(savedir) # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass

# Step 0: Get list of all variables and generate list of variables to remove from dataset.

## FUNCTION: Generate list of variables, drop variables from collection of files.
# Input: workdir and savedir. Output: 2 csvs: 1) of removed variables and 2) error file. Also returns list of removed variables.
# This function natively generates a list of all variables (including from groups), compares them against a list of variables to keep
# and prints a list of removed variables to a csv in the save directory. This should be run occasionally but not daily.

def get_vars(homedir, workdir, savedir, **options):
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download
    os.chdir(homedir)
    os.chdir(workdir) # Change directory to where raw files saved.
    files = os.listdir() # Gets list of files in directory to work with
    files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names
    files_old = [i for i in files if any(i for j in list(range(1980, 2011)) if str(j) in i)] # Get all file names which include any year from 1980-2010.
    files_new = [i for i in files if i not in files_old] # Get all other file names.
    keys = []
    coords = [] # Coordinates need to be dropped seperately. Get list of coordinates.
    dropvars = []
    badvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.

    for file in files:
        if file in files_new: # If file from 2011 or later, process to remove group structure.
            try:
                # Get global attributes from main file
                ds = xr.open_dataset(file, drop_variables=badvars)
                
                # Get list of sensors from global attributes
                # Some sensors have a 'sensor_suite' attribute, which we use to get the exact paths of sensors.
                if 'sensor_suite' in ds.attrs:
                    sensor_list = ds.attrs['sensor_suite'].split(',')
                    sensor_list = [k for k in sensor_list if k.endswith("1")] # Get data from primary sensor - may change depending on email response.
                    sensor_list = sorted(sensor_list)
                    sensors = set([i[11:] for i in sensor_list]) # Get unique list of sensors found in dataset.
                    paths = []
                    # Decision: Only keep first instance of a sensor if reported across multiple payloads. [ revisit pending email response ]
                    for sensor in sensors:
                        index = [idx for idx, s in enumerate(sensor_list) if str(sensor) in s][0]
                        path = sensor_list[index]
                        paths.append(path)
                    

                else: # Otherwise, try all possible sensors manually. To do this, combine payload 1-5(?) with sensor 1.
                    payloads = ['/payload_1/', '/payload_2/', '/payload_3/', '/payload_4/', '/payload_5/']
                    sensors = ['anemometer_1', 'humidity_sensor_1', 'ocean_temperature_sensor_1', 'air_temperature_sensor_1', 'barometer_1', 'gps_1', 'wave_sensor_1']
                    if options.get("sensors") == "secondary": # If function specified to download backup sensors.
                        sensors = ['anemometer_2', 'humidity_sensor_2', 'ocean_temperature_sensor_2', 'air_temperature_sensor_2', 'barometer_2', 'gps_1', 'wave_sensor_2']
                    paths = [''.join(chunks) for chunks in itertools.product(payloads, sensors)] # Generate all possible combinations of paths.

                for path in paths: # For each sensor
                    try:
                        ds_temp = xr.open_dataset(file, drop_variables=badvars, group = path) # Open sensor data
                    except Exception as e:
                        #print(file, e) # To revisit: not reporting error here because some groups don't have sensor list (and so we have to try the full suite.)
                        continue # Continue to next branch
                    ds = xr.merge([ds, ds_temp], compat = 'override') # Move sensor data to main ds. Where sensors have the same name, only override if one sensor has missing data.
                    del(ds_temp)
            except Exception as e:
                #print(file, e) # For testing. Remove in final run.
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue
        elif file in files_old:
            try:
                ds = xr.open_dataset(file, drop_variables=badvars)
            except Exception as e:
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue
        try:
            key = list(ds.keys())
            coord = list(ds.coords)
            keys+=key
            coords+=coord
        except Exception as e:
            errors['File'].append(file)
            errors['Time'].append(end_api)
            errors['Error'].append(e)
            continue

    # Get list of all variables
    mar_vars = list(set(keys)) # Get list of unique vars 
    mar_vars = mar_vars + badvars # add troublesome vars back in.

    #print(sorted(mar_vars)) # Use to manually scan list and identify variables to keep

    # Identify all variables to keep (this is a very conservative list, could prob cut more.)
    vars_to_keep = ['latitude',
                    'longitude',
                    'solar_radiation_36', # Solar radiation
                    'wind_speed', # Wind speed
                    'speed_averaging_method', # Info on variable calculation method (vector or scalar)
                    'air_temperature', # Air temp
                    'wind_sampling_duration', # Duration of averaging for wind speed.
                    'air_pressure_at_sea_level', # Air pressure 
                    'air_pressure', ## QUESTION - air pressure at sea level or air pressure??
                    #'time10_bnds', 'timem_bnds', 'time_bnds', # Time bounds - not saved for now as they cause errors.
                    'dew_point_temperature',  # Dew point
                    'precipitation', # Precipitation
                    'wind_direction', # Wind direction
                    'relative_humidity',
                    'anemometer_height',
                    ] 

    vars_to_keep = vars_to_keep + [item + '_qc' for item in vars_to_keep] + [item + '_detail_qc' for item in vars_to_keep] # Add qc flags.

    # Remove these variables from the list of all variables to get drop list.
    dropvars = np.setdiff1d(mar_vars, vars_to_keep)

    # Write the list of removed variables to csv for future reference. Keep these in one centralized "removedvars.csv" file that gets appended to.
    vars = []
    try:
        os.chdir(homedir+'/'+savedir)
        # Do the reading
        if os.path.isfile('removedvars.csv'): # If removedvars.csv already exists, read in data to vars. Otherwise, vars is empty.
            with open('removedvars.csv','r') as f: 
                rows = list(csv.DictReader(f, delimiter=','))
                for row in rows:
                    vars.append(row['Variable'])
        
        for i in dropvars:
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
        
    # Write errors to a csv.
    filepath = "errors_maritime_vars_{}.csv".format(end_api) # Set path to save error file.
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

    # Return to home directory
    os.chdir(homedir)

    return dropvars

## FUNCTION: Generate list of station and instrument elevations.
# Input: url to tables on NDBC website. Output: returns dataframe of stations and elevations. Called internally in clean_maritime() function.
# This function generates a dataframe of station elevations and instrument elevations from the NDBC website, as this data is frequently missing from the data source.
# This is joined by station to the xarray dataframes in the clean_maritime() function.

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
    df.columns= columns
    df['Tide_Reference'] = np.NAN # Add 6th column.
    df = df.reindex(columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference","Barometer_Elevation"])
    
    # Table 2 has 6 columns
    tabletext = tables[1]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference","Barometer_Elevation"]
    table = tabletext.get_text().rsplit('ELEVATION',1)[1] # Remove headers.
    table = table.split() # Remove whitespace.
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+6] for x in range(0, len(table),6)] # Split into rows.
    dftemp = pd.DataFrame(composite_list) # Turn to dataframe.
    dftemp.columns= columns
    df = pd.concat([df, dftemp])
    
    # Table 3 has 9 columns, but we only want the first 3.
    tabletext = tables[2]
    columns = ["Station_ID", "Site_Elevation", "Air_Temp_Elevation", "Anemometer_Elevation", "Tide_Reference","Barometer_Elevation"]
    table = tabletext.get_text().rsplit('CIRCLE',1)[1] # Remove headers.
    table = table.split() # Remove whitespace.
    # Should be 5 for table 0 and 6 for table 1+2
    composite_list = [table[x:x+9] for x in range(0, len(table),9)] # Split into rows.
    dftemp = pd.DataFrame(composite_list) # Turn to dataframe.
    # Drop last three columns.
    dftemp = dftemp.iloc[:,0:6]
    dftemp.columns= columns
    df = pd.concat([df, dftemp])
    return df
    
# FUNCTION: Get list of ASOS/AWOS IDs in ISD dataset. 
# Pulls in ISD station list and ASOSAWOS station list (two separate csvs), joins by ICAO and returns list of station IDs.
def get_ids(homedir):
    os.chdir("test_platform/scripts/2_clean_data/")
    isd = pd.read_csv("isd-history.csv")
    isd = isd[isd['CTRY']=="US"]

    asosawos = pd.read_csv('asosawos_stations.csv')
    asosawos = asosawos.loc[(asosawos['A']=="A") | (asosawos['A']=="W")]

    # Use ICAO as matching column because neither datasource has any NA values here.
    # Fix column types.
    asosawos['ICAO'] = asosawos['ICAO'].astype(str)
    isd['ICAO'] = isd['ICAO'].astype(str)
    # Merge on ICAO.
    join = pd.merge(asosawos, isd, on = 'ICAO', how = 'inner')
    # Get station IDs from USAF and WBAN id columns
    join['station-id'] = join['USAF'].astype(str)+"-"+join['WBAN'].astype(str)
    os.chdir(homedir) # Return to homedirectory
    return join['station-id']

test = get_ids(homedir)
print(test)
# missing['LAT_clean'] = missing['LAT_clean'].round(3)
# missing['LON_clean'] = missing['LON_clean'].round(3)
# print(missing)

# test2 = pd.merge(missing, isd, left_on=['LAT_clean', 'LON_clean'], right_on=["LAT", "LON"], how = "inner") # About 9 more stations.
# print(test2)
# Try match by name.
#test2 = pd.merge(missing, isd, left_on = 'STATION', right_on = "STATION NAME", how = 'inner') # 10 more match.
#print(test2)
# # Group by station and get last ending date.
# finaldate = test.loc[test.groupby('ICAO').END.idxmax()]
# print(finaldate['END'].min())
# print(len(finaldate[(finaldate['END'].astype('str').str.contains("2022")==False)]))

# finaldate['END'] = pd.to_datetime(finaldate["END"])
# finaldate.hist(column = 'END')
#print(asosawos)

## FUNCTION: Clean MARITIME data.
# Input: workdir, savedir.
# Note that this function searches for a csv called removevars.csv in the savedir folder. If file doesn't exist, script will call get_vars function.
def clean_maritime(homedir, workdir, savedir, **options):
    
    # Set directory to parent git directory (historical-obs-platform folder)
    os.chdir(homedir)

    # Set up bounding box to filter data.
    try:
        wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp' # Harded-coded because these should not move.
        wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 
        t, m, bbox = get_wecc_poly(wecc_terr, wecc_mar)
        lonmin, lonmax = float(bbox['minx']), float(bbox['maxx']) 
        latmin, latmax = float(bbox['miny']), float(bbox['maxy']) 
    except: # If geospatial call fails, hardcode.
        lonmin, lonmax = -139.047795, -102.03721
        latmin, latmax = 30.142739, 60.003861

    ## Set up csv to record any files not cleaned and reason.
    errors = {'File':[], 'Time':[], 'Error':[]}
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    # Read in sensor heights from NDBC website.
    # From https://www.ndbc.noaa.gov/bmanht.shtml
    try:
        url = "https://www.ndbc.noaa.gov/bmanht.shtml"
        elevs = get_elevs(url)
    except Exception as e: # May want to update this error handling, this will kill function if it does not work.
        print(e)
        return

    ## Step 1: Read in files and first pass of clean and organize
    timestamp = datetime.now()

    # Set target dimension and variable order
    # For example, air_temperature would have (time x lat x lon x elev x data) dimensions
    clean_dims = ['time', 'latitude', 'longitude', 'elevation']
    clean_vars = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']

    # Read in files.
    os.chdir(workdir) # Change directory to where raw files saved.
    files = os.listdir() # Gets list of files in directory to work with
    files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names

    ### Note: Useful for code where stations are together, but not necessary for this script. Leaving here commented in case needed.
    # # Get list of station IDs from filename and clean.
    ids = list()
    for file in files:
        id = re.sub('_[0-9]{6,6}.*.nc', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
        id = re.sub('NDBC_', "", id) # Remove leading NDBC
        id = re.sub("_", "", id) # Remove any remaining underscores
        if id not in ids:
            ids.append(id)

    # ## Step 2: Read in datafile and drop variables that aren't of interest
    # # # Read in list of column variables to drop. If no list, call get_vars().
    # # # This method is robust to variables getting added over time.
    
    os.chdir(homedir)
    dropvars = []
    path = savedir+"/removedvars.csv" 
    if os.path.exists(path): # If path to csv exists, attempt to read it. 
        with open(path,'r') as f: 
                rows = list(csv.DictReader(f, delimiter=','))
                for row in rows:
                    dropvars.append(row['Variable'])
    else: # Otherwise, call get_vars function.
        dropvars = get_vars(homedir, workdir, savedir)
    
    # Identify list of coordinates to drop
    coords_to_remove = ['depth', 'timem', 'time_wpm_20', 'time_wpm_40', 'time_bnds', 'time10', 'wave_wpm'] 
    # Hardcoded for now, could get moved. These are a list of coords that throw errors consistently.
    
    # # Split files into those written before 2015 and those written after.
    files_old = [i for i in files if any(i for j in list(range(1980, 2011)) if str(j) in i)] # Get all file names which include any year from 1980-2010.
    files_new = [i for i in files if i not in files_old] # Get all other file names.
    
    # Read in data files by station ID, clean, and merge.
    os.chdir(workdir)
    for i in ids: # For each station
    #for i in ['46023', '46125', '46027']: # For testing
        ds_stat = None # Initialize ds_stat
        file_count = 0
        stat_files = [k for k in files if i in k] # Get list of files with station ID in them.
        for file in stat_files: # For each file
            try:
                if file in files_new: # If file newer than 2015, process to remove group structure.
                    
                    # Get global attributes from main file
                    ds = xr.open_dataset(file, drop_variables=dropvars)

                    ## Before we clean, require lat/lon coordinates and filter by location in WECC.
                    if 'geospatial_lat_min' and 'geospatial_lat_max' and 'geospatial_lon_min' and 'geospatial_lat_max' in ds.attrs:
                        if float(ds.attrs['geospatial_lat_min']) < latmin or float(ds.attrs['geospatial_lat_max']) > latmax or float(ds.attrs['geospatial_lon_min']) < lonmin or float(ds.attrs['geospatial_lon_max']) > lonmax:
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(float(ds.attrs['geospatial_lat_min']), float(ds.attrs['geospatial_lon_min'])))
                            continue
                        else:
                            pass 
                            #print("File {} in WECC!".format(file)) # Just for testing, remove for final.
                    else:
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("No geospatial coordinates found in file")
                        continue

                    # Get list of sensors from global attributes
                    # Some sensors have a 'sensor_suite' attribute, which we use to get the exact paths of sensors.
                    if 'sensor_suite' in ds.attrs:
                        full_sensor_list = ds.attrs['sensor_suite'].split(',')
                        sensor_list = [k for k in full_sensor_list if k.endswith("1")] # Get data from primary sensors
                        if options.get("sensors") == "secondary":
                            sensor_list = [k for k in full_sensor_list if k.endswith("2")] # Get data from secondary sensors
                            sensor_list.append("gps_1") # Need lat and lon for script to work, there is no gps 2 typically.
                        sensor_list = sorted(sensor_list)
                        sensors = set([i[11:] for i in sensor_list]) # Get unique list of sensors found in dataset.
                        paths = []
                        # Decision: Only keep first instance of a sensor if reported across multiple payloads.
                        # This decision is robust to the E & M sensor types discussed in email communications with Dawn, 
                        # where E is primary and M is secondary, because E sensors are always on payloads preceding M sensors.
                        for sensor in sensors:
                            index = [idx for idx, s in enumerate(sensor_list) if str(sensor) in s][0]
                            path = sensor_list[index]
                            paths.append(path)
                        
                    else: 
                        payloads = ['/payload_1/', '/payload_2/', '/payload_3/', '/payload_4/', '/payload_5/']
                        sensors = ['anemometer_1', 'humidity_sensor_1', 'ocean_temperature_sensor_1', 'air_temperature_sensor_1', 'barometer_1', 'gps_1', 'wave_sensor_1']
                        if options.get("sensors") == "secondary":
                            sensors = ['anemometer_2', 'humidity_sensor_2', 'ocean_temperature_sensor_2', 'air_temperature_sensor_2', 'barometer_2', 'gps_1', 'wave_sensor_2'] # GPS as 1 because we need lat and lon always.
                        
                        # Option 1: Simple, easy, works with "override" feature of xarray.merge to keep the first instance of a variable 
                        # (by default, the one from the lowest number payload.)
                        # Ignores the "M" feature altogether (I think impact here is minimal, see notes in option 2)
                        paths = [''.join(chunks) for chunks in itertools.product(payloads, sensors)] # Generate all possible combinations of paths.


                    if paths is None: # If no sensors, go to next file.
                        continue

                    for path in paths: # For each sensor
                        try:
                            ds_temp = xr.open_dataset(file, drop_variables=dropvars, group = path) # Open sensor data and drop variables from dropvar list.
                            
                        except Exception as e:
                            #print(file, e) # To revisit: not reporting error here because some groups don't have sensor list (and so we have to try the full suite.)
                            continue # Continue to next branch
                        
                        ds = xr.merge([ds, ds_temp], compat = 'override') # Move sensor data to main ds
                        
                        # Add sensor attributes to main dataset if they exist.
                        if 'manufacturer' and 'part_number' in ds_temp.attrs:
                            sensor_name = str(ds_temp.attrs['manufacturer'])+" "+str(ds_temp.attrs['part_number'])# Get sensor data as string from dataframe.
                            sensor_type = re.sub('[^a-zA-Z]+', '', path.rsplit('/', 1)[-1]) # Remove all spaces and numbers
                        elif 'manufacturer' in ds_temp.attrs:
                            sensor_name = str(ds_temp.attrs['manufacturer'])# Get sensor data as string from dataframe.
                            sensor_type = re.sub('[^a-zA-Z]+', '', path.rsplit('/', 1)[-1])
                        else:
                            sensor_type = re.sub('[^a-zA-Z]+', '', path.rsplit('/', 1)[-1])
                            sensor_name = "Unknown"
                        
                        sensor_name = [sensor_name]*len(ds['time'])
                        ds[sensor_type] = ('time', sensor_name) # Assign sensor manufacturing info as attribute.
                        
                        # Install date not included in dataset but available in original metadata.
                        # If needed, will be inferred later from sensor column.

                        if 'height_of_instrument' in ds_temp.attrs: # Make height of instrument an attribute (when taken directly from metadata)
                            # In QA/QC, compare this value to the one taken from the function get_elev().
                            if 'height_of_instrument' in ds_temp.attrs:
                                sensor_height = sensor_type+"_height"  
                                ds.attrs[sensor_height] = ds_temp.attrs['height_of_instrument'] # Assign sensor height as attribute where provided.
                        del(ds_temp)

                    # If after merging, ds is empty, record error and skip file.
                    if not list(ds.data_vars):
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("No variables recorded in file")
                        continue
                    
                # For older files, simply read in and drop vars from dropvar list.
                if file in files_old:
                    
                    ds = xr.open_dataset(file, drop_variables=dropvars)
                    
                    #if options.get("sensors") == "secondary": # If function specified to download backup sensors.
                    #    pass # TO DO - Are there secondary sensors in pre 2010 files? Waiting on email response.
                        #print(ds.attrs)
                        #print(ds.keys())
                        #continue # Older files do not have secondary sensors????, so simply skip them.

                    # Quality control
                    ## Before we clean, skip any files that don't have data in them (and save to error list).
                    if not list(ds.data_vars):
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("No variables recorded in file")
                        continue

                    ## Before we clean, require lat/lon coordinates and filter by location in WECC.
                    if 'lat' or 'lon' in ds.keys():
                        #Filter to only keep lat and lon in WECC.
                        if float(max(ds['lat'])) > latmax or float(min(ds['lat'])) < latmin or float(max(ds['lon'])) > lonmax or float(min(ds['lon'])) < lonmin:
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(float(max(ds['lat'])), float(min(ds['lat']))))
                            continue
                        elif float(ds.attrs['geospatial_lat_min']) < latmin or float(ds.attrs['geospatial_lat_max']) > latmax or float(ds.attrs['geospatial_lon_min']) < lonmin or float(ds.attrs['geospatial_lon_max']) > lonmax:
                            errors['File'].append(file)
                            errors['Time'].append(end_api)
                            errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(float(ds.attrs['geospatial_lat_min']), float(ds.attrs['geospatial_lon_min'])))
                            continue
                        #else:
                        #    print("File {} in WECC!".format(file)) # Just for testing, remove for final.
                    else:
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("No geospatial coordinates found in file")
                        continue

                # Now, both old and new files are in ds format.
                # Remove unwanted coordinates.
                ctr = [x for x in coords_to_remove if x in ds.keys()] # Filter coords to remove list by coordinates found in dataset.
                ds = ds.drop_dims(ctr) # Drop coordinates.
                
                ## Step 3: Convert station metadata to standard format
                
                # 3.1  Organize netCDF dimensions and coordinates.
                ## Two dims - station ID and time.
                # lat and lon as coordinates dependent on station and time.
            
                # Generate unique ID across entire dataset by combining network name and ID.
                id = "MARITIME_"+ds.attrs["id"]
                ds['original_id'] = ds.attrs['id'] # Keep original ID as var.
                    
                if file in files_old: # Works to reorder old files.
                    # Convert lat to array to make it a function of time
                    lat = [float(ds['lat'])]
                    lat = np.asarray(lat*len(ds['time']))
                    lat.shape = (1, len(ds['time']))
                    ds = ds.reset_coords('lat', drop = True)

                    # Convert lon to array to make it a function of time
                    lon = [float(ds['lon'])]
                    lon = np.asarray(lon*len(ds['time']))
                    lon.shape = (1, len(ds['time']))
                    ds = ds.reset_coords('lon', drop = True)

                    # Reconfigure station ID such that ID, rather than a number, is the index.
                    ds = ds.assign_coords(station_id = str(id))
                    ds = ds.expand_dims("station_id")
                    ds = ds.squeeze("station")
                    ds = ds.reset_coords("station", drop = True)
                    ds = ds.rename({"station_id": "station"})

                    # reassign lat and lon as coordinates
                    ds = ds.assign_coords(lat = (["station","time"],lat), lon = (["station","time"], lon))

                if file in files_new:
                    try:
                        ds = ds.rename({'latitude': 'lat', 'longitude': 'lon'})
                    except:
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("No geospatial coordinates found in file")
                        continue
                    ds = ds.assign_coords(station = str(id))
                    ds = ds.expand_dims("station")
                    ds = ds.assign_coords(lat = (["station","time"],ds['lat'].data), lon = (["station","time"], ds['lon'].data))
                    

                # 3.3 Set elevation for station and sensors
                stat_elevs = elevs[elevs["Station_ID"] == ds["original_id"]]
                
                if stat_elevs.empty == False:
                    
                    #Fix NAs
                    stat_elevs = stat_elevs.replace("NA", np.nan)
                    
                    # reassign station and sensor elevations as coordinates
                    elev = np.asarray([float(stat_elevs['Site_Elevation'])]*len(ds['time']))
                    tas_elev = np.asarray([float(stat_elevs['Air_Temp_Elevation'])]*len(ds['time']))
                    anem_elev = np.asarray([float(stat_elevs['Anemometer_Elevation'])]*len(ds['time']))
                    bar_elev = np.asarray([float(stat_elevs['Barometer_Elevation'])]*len(ds['time']))
                
                    
                    
                elif stat_elevs.empty: # A few stations are missing in the table. If no value returned, insert NAs.
                    #print("Station elevations not found") # Helpful for ID'ing stations with missing data, comment out for final run.
                    elev = np.asarray([np.nan]*len(ds['time']))
                    tas_elev = np.asarray([np.nan]*len(ds['time']))
                    anem_elev = np.asarray([np.nan]*len(ds['time']))
                    bar_elev = np.asarray([np.nan]*len(ds['time']))
                    
                # Reshape arrays
                elev.shape = (1, len(ds['time']))
                tas_elev.shape = (1, len(ds['time']))
                anem_elev.shape = (1, len(ds['time']))
                bar_elev.shape = (1, len(ds['time']))
                    
                # Add variables
                ds['elevation'] = (['station', 'time'], elev)
                ds['tas_elev'] = (['station', 'time'], tas_elev)
                ds['anemometer_elev'] = (['station', 'time'], anem_elev)
                ds['barometer_elev'] = (['station', 'time'], bar_elev)
                
                # 3.4 Standardize sensor metadata
                ### Note that anemometer height from dataset (as att) and elev (as var) from table are not equivalent.
                #  Not used currently, but flag for future ref.

                ## Step 4: Convert dataset metadata in standard format -- CF compliance - to be finalized, overwrite existing metadata.
                ds = ds.assign_attrs(title = "MARITIME cleaned")
                ds = ds.assign_attrs(institution = "Eagle Rock Analytics / Cal Adapt")
                ds = ds.assign_attrs(source = "")
                ds = ds.assign_attrs(history = "MARITIME_clean.py")
                ds = ds.assign_attrs(comment = "Intermediate data product: may not have been subject to any cleaning or QA/QC processing")
                ds = ds.assign_attrs(license = "")
                ds = ds.assign_attrs(citation = "")
                ds = ds.assign_attrs(disclaimer = "This document was prepared as a result of work sponsored by the California Energy Commission (PIR-19-006). It does not necessarily represent the views of the Energy Commission, its employees, or the State of California. Neither the Commission, the State of California, nor the Commission's employees, contractors, or subcontractors makes any warranty, express or implied, or assumes any legal liability for the information in this document; nor does any party represent that the use of this information will not infringe upon privately owned rights. This document has not been approved or disapproved by the Commission, nor has the Commission passed upon the accuracy of the information in this document.")
                
                ## Step 5: Convert missing data to common format -- CF compliance
                for var in ds.variables: 
                    try:
                        ds = ds.where(ds[var] != 9.969209968386869e+36) # Replace any identified NA values by using mask.
                        #print([var, float(ds[var].min()), float(ds[var].max())]) 
                        # ^ Code to run to print min and max to visually ID non-NAN NA values. Commented out but useful for future scripts.
                    except:
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("Error removing masking NA values in variables..")
                        continue  
                
                # Where lat or lon is = this value, remove observation.
                try:
                    ds = ds.where(ds['lat']!= 9.969209968386869e+36, drop = True)
                    ds = ds.where(ds['lon']!= 9.969209968386869e+36, drop = True)
                except:
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Error removing observations missing lat/lon coordinates.")
                    continue
                
                # Code to test this works.
                # for var in ["lat", "lon"]: 
                #     try:
                #         print([var, float(ds[var].min()), float(ds[var].max())]) 
                #         # ^ Code to run to print min and max to visually ID non-NAN NA values. Commented out but useful for future scripts.
                #     except Exception as e:
                #         print(e)
                #         break  
                
                
                # Step 6: Standardize variables and variable metadata.
                # If not observed, calculate derived primary variables if possible.
                #** indicates primary approach
                   
                # tas : air surface temperature
                if "air_temperature" in ds.keys():
                    ds['tas'] = ds["air_temperature"]+273.15 # Convert from C to K
                    ds['tas'].attrs['units'] = "degree_Kelvin" # Update units
                    ds = ds.rename({'air_temperature': 'tas_raw'})

                # ps: surface air pressure
                ## Two variables provided for some datasets, and they should be identical (sensors at sea level).
                ## Take "air pressure" as primary, and use converted "air pressure at sea level" if not.

                if "air_pressure" in ds.keys():
                    ds['ps'] = ds["air_pressure"]*100 # Convert from mb to Pa
                    ds['ps'].attrs['units'] = "Pascal" # Update units
                    ds = ds.rename({'air_pressure': 'ps_raw'}) # Convert from mb to Pa
                    if "air_pressure_at_sea_level" in ds.keys():
                        ds = ds.drop("air_pressure_at_sea_level") # Remove second air pressure var if it exists.

                elif "air_pressure_at_sea_level" in ds.keys():
                    # Note there are a few stations which don't have barometer elevation from chart or sensor data.
                    # These cannot be converted, but we leave the raw data here.
                    if 'barometer_elev' in ds.keys():
                        if stat_elevs['Barometer_Elevation'] is not np.nan:
                            # Convert to air pressure at station level.
                            # Inputs: psl, elev, temp
                            ps = _calc_ps(ds['air_pressure_at_sea_level'], ds['barometer_elev'], ds['tas'])
                            if np.isnan(ps).all(): # If all values NA, don't save as variable, just keep raw data.
                                ds = ds.rename({'air_pressure_at_sea_level': 'psl_raw'})
                            else:
                                ds['ps'] = ps*100 # Convert from mb to Pa
                                ds['ps'].attrs['units'] = "Pascal" # Update units
                                ds = ds.rename({'air_pressure_at_sea_level': 'psl_raw'}) # Rename original obs.

                    elif 'barometer_height' in ds.attrs:
                        # Haven't seen an instance where barometer height gleaned from metadata but station not in table.
                        # Print statement here should catch this if so, remove after testing if not triggered.
                        print(ds.attrs['barometer_height'])
                        #ps = _calc_ps(ds['air_pressure_at_sea_level'], ds.attrs['barometer_height'], ds['tas'])
                        break  # For testing.
                                
                # tdps: dew point temperature
                # dew point temperature calculation (necessary input vars: requires at least 2 of three - air temp + relative humidity + vapor pressure)
                # Only more recent stations have dewpoint temp (dew_point_temperature)
                # No stations have vapor pressure or RH, so no calculations possible.
                if "dew_point_temperature" in ds.keys(): # If variable already exists, rename.
                    ds['tdps'] = ds["dew_point_temperature"]+273.15 # Convert from C to K
                    ds['tdps'].attrs['units'] = "degree_Kelvin"
                    ds = ds.rename({'dew_point_temperature': 'tdps_raw'})
                
                # # pr: precipitation 
                # DO CONVERSIONS FROM RAINFALL TO RAINFALL RATE (Hourly accum) DURING NEXT STAGE, as we will aggregate over the hour.
                if "precipitation" in ds.keys():
                    ds = ds.rename({'precipitation': 'p_raw'})

                # hurs: relative humidity
                # relative humidity calculation (necessary input vars: air temp + dew point**, air temp + vapor pressure, air pressure + vapor pressure)
                ## Vapor pressure not a variable measured by this network.
                if 'relative_humidity' in ds.keys():
                    ds = ds.rename({'relative_humidity': 'hurs'})
                    ds['hurs_raw'] = ds['hurs'] # Add raw column if data originally provided.
                elif 'tas' and 'tdps' in ds.keys():
                    ds['hurs'] = _calc_relhumid(ds['tas'], ds['tdps'])

                # # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation)
                # Note that this occurs so infrequently that I still need to check units are correct. Thus print statement left in for now.
                if 'solar_radiation_36' in ds.keys():
                    ds = ds.rename({'solar_radiation_36' : 'rsds'}) # To do.
                    print(ds['rsds'].attrs)
                    
                # # sfcWind : wind speed
                ### In maritime we have:
                ## wind_speed: average wind speed (with duration of averaging specified in wind_sampling_duration)
                ## speed_averaging_method: contains QAQC flags.
                if "wind_speed" in ds.keys(): # If variable already exists, rename. Units already in m s-1.
                    ds = ds.rename({'wind_speed': 'sfcWind'})
                
                # # sfcWind_dir: wind direction
                if "wind_direction" in ds.keys(): # If variable already exists, rename.
                    ds = ds.rename({'wind_direction': 'sfcWind_dir'}) # Already in degrees.
                
                # ## Step 6: Tracks existing QA/QC flags to standard format
                # Files older than 2015 have no built in QA/QC flags (as far as I see).
                ## Keep both flags (detail and qc) for now.

                # ps: surface air pressure
                if 'air_pressure_qc' in ds.keys():
                    ds = ds.rename({'air_pressure_qc': 'ps_qc',
                                    'air_pressure_detail_qc': 'ps_detail_qc'})
                    if "air_pressure_at_sea_level_qc" in ds.keys():
                        ds = ds.drop("air_pressure_at_sea_level_qc")
                        ds = ds.drop("air_pressure_at_sea_level_detail_qc")

                elif 'air_pressure_at_sea_level_qc' in ds.keys():
                    ds = ds.rename({'air_pressure_at_sea_level_qc': 'ps_qc',
                                    'air_pressure_at_sea_level_detail_qc': 'ps_detail_qc'})
                
                # tas : air surface temperature
                if 'air_temperature_qc' in ds.keys():
                    ds = ds.rename({'air_temperature_qc': 'tas_qc',
                                    'air_temperature_detail_qc': 'tas_detail_qc'})
                
                # tdps: dew point temperature
                if 'dew_point_temperature_qc' in ds.keys():
                    ds = ds.rename({'dew_point_temperature_qc': 'tdps_qc',
                                    'dew_point_temperature_detail_qc': 'tdps_detail_qc'})
                
                # # pr: precipitation 
                if 'precipitation_qc' in ds.keys():
                    ds = ds.rename({'precipitation_qc': 'pr_qc',
                                    'dew_point_temperature_detail_qc': 'pr_detail_qc'})

                # hurs: relative humidity
                if 'relative_humidity_qc' in ds.keys():
                    ds = ds.rename({'relative_humidity_qc': 'hurs_qc',
                                    'relative_humidity_detail_qc': 'hurs_detail_qc'})
                
                # # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation)
                if 'solar_radiation_36_qc' in ds.keys():
                    ds = ds.rename({'solar_radiation_36_qc': 'rsds_qc',
                                    'solar_radiation_36_detail_qc': 'rsds_detail_qc'})
                
                # # sfcWind : wind speed
                # Add wind speed calculation method to flags here.
                if 'wind_speed_qc' in ds.keys():
                    ds = ds.rename({'wind_speed_qc': 'sfcWind_qc',
                                    'wind_speed_detail_qc': 'sfcWind_detail_qc'})

                # # sfcWind_dir: wind direction
                if 'wind_direction_qc' in ds.keys():
                    ds = ds.rename({'wind_direction_qc': 'sfcWind_dir_qc',
                                    'wind_direction_detail_qc': 'sfcWind_dir_detail_qc'})
            
                # latitude
                if 'latitude_qc' in ds.keys():
                    ds = ds.rename({'latitude_qc': 'lat_qc',
                                    'latitude_detail_qc': 'lat_detail_qc'})
            
                # longitude
                if 'longitude_qc' in ds.keys():
                    ds = ds.rename({'longitude_qc': 'lon_qc',
                                    'longitude_detail_qc': 'lon_detail_qc'})

            
                # Final cleaning steps.

                # Anemometer height: make either into attribute or into var across all dfs.
                # Speed averaging method and wind sampling duration - rename to conform or make attrs if the same across all values?
                #print(np.unique(ds['wind_sampling_duration']))
                #print(np.unique(ds['speed_averaging_method']))

                #print(list(ds.variables))
                # REMOVE ANY EXTRA METADATA - to do.
                # Remove any pesky coordinates.
                if 'time_wpm_40' in ds.coords:
                    ds = ds.drop_dims('time_wpm_40')
                if 'time_bnds' in ds.coords:
                    ds = ds.drop_dims('time_bnds')
                    print(ds)

                # # Reorder variables
                # In following order:
                desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
                desired_order = [i for i in desired_order if i in list(ds.keys())] # Only keep vars which are in ds.
                rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom.
                new_index = desired_order + rest_of_vars
                ds = ds[new_index]

                ## Step 7: Merge files by station

                # Merge file to previous time period data for same station.
                if ds_stat is None:
                    ds_stat = ds # For first file in station, set ds_stat to be original ds.
                else:
                    ds_stat = xr.merge([ds_stat,ds]) # Otherwise, merge new records to first ds.
                    if ds_stat.attrs["time_coverage_start"]>ds.attrs["time_coverage_start"]: # If dataset merged begins before ds_stat
                        ds_stat.attrs["time_coverage_start"]=ds.attrs["time_coverage_start"] # Update time coverage to match.
                    elif ds_stat.attrs["time_coverage_end"]<ds.attrs["time_coverage_end"]: # Otherwise if dataset merged begins after ds_stat
                        ds_stat.attrs["time_coverage_end"]=ds.attrs["time_coverage_end"] # Update time coverage to match. 
                    file_count +=1
                    ds_stat.attrs['raw_files_merged'] = file_count # Keep count of how many files merged per station.

                # Testing only: delete later.
                # time = ds.attrs["time_coverage_start"][:7]
                # #print(timestamp)
                # filename = ds_stat.attrs["station"]+"_"+time+".nc" # Make file name
                # if options.get("sensors") == "secondary": # If function specified to download backup sensors.
                #     filename = ds_stat.attrs["station"]+"_"+time+"_secondary.nc"
                # filepath = homedir+"/"+savedir+filename # Write file path
                # #print(filepath)
                # ds.to_netcdf(path = filepath) # Save station file.
                # print("Saving {}".format(filename))
                
            except Exception as e:
                print(file, e) # For testing only.
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue

        # For each station, save one file.
        # Write files to netCDF
        if ds_stat is None: # If file is skipped, ds_stat will be equivalent to none.
            print("File {} not saved.".format(file))
            continue
        else:
            try:
                filename = id+".nc" # Make file name
                if options.get("sensors") == "secondary": # If function specified to download backup sensors.
                    filename = id+"_secondary.nc"
                filepath = homedir+"/"+savedir+filename # Write file path
                ds_stat.to_netcdf(path = filepath) # Save station file.
                print("Saving {} with dims {}".format(filename, ds_stat.dims))
                ds.close() # Close dataframe.
            except Exception as e:
                print(filename, e)
                errors['File'].append(filename)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue

    #Write errors to csv
    filepath = homedir+"/"+savedir+"errors_maritime_{}.csv".format(end_api) # Set path to save error file.
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

#clean_maritime(homedir, workdir, savedir)   
#clean_maritime(homedir, workdir, savedir, sensors = "secondary")   # For backup data.


## notes
# 
# Option 2: Incorporating "M" data. (Scrapped for now, could get used in future iterations) - for sensor/payload merging.
                        # Future work could also incorporate release variables.

                        # I wrote this method
                        # paths = []
                        # for sensor in sensors:                
                        #     for payload in payloads:
                                
                        #         try:
                        #             # For each sensor, open payloads in order.
                        #             ds_temp = xr.open_dataset(file, drop_variables=dropvars, group = payload) 
                        #             payload_id = ds_temp.attrs["payload_id"]
                        #             group = payload+sensor
                        #             ds_temp = xr.open_dataset(file, drop_variables=dropvars, group = group) # If file able to open
                        #             paths.append(group) # Append to path list
                        #             break # and move to next sensor
                        #         except:
                        #             # If there are no "sensor_2" for any given sensor in a payload, make a last check to see if any payloads are designated as "M" (secondary) in type.
                        #             if "M" in payload_id: # If payload designated as payload of secondary sensors, take sensor even if labelled as 1.
                        #                 try:
                        #                     print(payload_id)
                        #                     sensor = sensor.replace(sensor[len(sensor) - 1:], "1")
                        #                     print(payload_id, sensor) # Haven't yet seen this actually trigger so this works theoretically but needs to be tested still.
                        #                     group = payload+sensor
                        #                     ds_temp = xr.open_dataset(file, drop_variables=dropvars, group = group) # If file able to open
                        #                     paths.append(group)
                        #                 except:
                        #                     continue # Onto the next sensor.
                        #             continue

# TESTING.


# A station that should not have 'ps' bc barometer elevation doesn't exist from either source.
#os.chdir(savedir)
#test = xr.open_dataset("MARITIME_46023.nc")

# Open a 'normal' file and run through var values.
# os.chdir(savedir)
# test = xr.open_dataset("MARITIME_46040.nc")
# print(test)
# for var in test.variables: 
#     try:
#         print([var, float(test[var].min()), float(test[var].max())]) 
#     except:
#         continue

### NOTES on sensor structure.
# E sensors always precede M in payload.
# E sensors may not be the same as M sensors (e.g., the list of sensors might not be overlapping at all.)
# Strategic decision: always pull barometer 1 from payload 1, then barometer 1 from payload 2 if payload 1 doesn't have it.

# os.chdir(homedir)

# # Read in files.
# os.chdir(workdir) # Change directory to where raw files saved.
# files = os.listdir() # Gets list of files in directory to work with
# files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names

# # # Split files into those written before 2015 and those written after.
# files_old = [i for i in files if any(i for j in list(range(1980, 2011)) if str(j) in i)] # Get all file names which include any year from 1980-2014.
# files_new = [i for i in files if i not in files_old] # Get all other file names.

# badvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.
# payloads = ['/payload_1/', '/payload_2/', '/payload_3/', '/payload_4/', '/payload_5/']
# paths = []
# # Read in data files by station ID, clean, and merge. 
# for file in files_new: # If file newer than 2015, process to remove group structure.
#     sensors = ['anemometer_1', 'humidity_sensor_1', 'ocean_temperature_sensor_1', 'air_temperature_sensor_1', 'barometer_1', 'gps_1', 'wave_sensor_1']
#     for sensor in sensors:                
#         for payload in payloads:
#             try:
#                 # Get global attributes from main file
#                 ds = xr.open_dataset(file, drop_variables=badvars, group = payload)
#                 group = payload+sensor
#                 ds = xr.open_dataset(file, drop_variables=badvars, group = group)
#                 print(file, group)
#                 paths.append(group)
#                 break               
#             except:
#                 continue
# #print(paths)



        #         print(file, payload, ds.attrs['payload_id'])
        #         for j in sensors:
        #             group = payload+j
        #             try:
        #                 ds = xr.open_dataset(file, drop_variables=badvars, group = group)
        #                 print(j, group)
        #                 next
        #             except:
        #                 continue
        #     if "M" in ds.attrs["payload_id"]:
        #         print(file, payload, ds.attrs['payload_id'])
        #         sensors = ['anemometer_1', 'humidity_sensor_1', 'ocean_temperature_sensor_1', 'air_temperature_sensor_1', 'barometer_1', 'gps_1', 'wave_sensor_1']
        #         for j in sensors:
        #             group = payload+j
        #             try:
        #                 ds = xr.open_dataset(file, drop_variables=badvars, group = group)
        #                 print(j, group)
        #                 next
        #             except:
        #                 continue
        #         #group = payload+"air_temperature_sensor_1"
        #         #ds = xr.open_dataset(file, drop_variables=badvars, group = group)
        #         #print(ds)
        # except Exception as e:
        #     #print(e)
        #     continue

# ### TESTING.
# # Test vars
# #%%
# import matplotlib.pyplot as plt
# import xarray as xr
# savedir = "data/2_clean_wx/MARITIME/"
# os.chdir(savedir)
# ds = xr.open_dataset("MARITIME_46023.nc")
# #ds['ps'].plot() # Looks good.
# ds['tdps'].plot()
# test = ds['tdps'].where(ds.tdps<9.96921e+36, drop = True)
# test.plot()
# # # Test merge.
# os.chdir(savedir)
# ds = xr.open_dataset("MARITIME_46023.nc")
# ds1 = xr.open_dataset("MARITIME_46023_198901.nc")
# ds2 = xr.open_dataset("MARITIME_46023_198902.nc")
# ds3 = xr.open_dataset("MARITIME_46023_198903.nc")
# ds4 = xr.open_dataset("MARITIME_46023_198904.nc")
# ds5 = xr.open_dataset("MARITIME_46023_198905.nc")
# ds6 = xr.open_dataset("MARITIME_46023_198906.nc")

#print(ds)
#ds_test = xr.merge([ds1,ds2,ds3,ds4,ds5,ds6])

# Check they are the same - all values should return true.
# for key in ds.keys():
#     print ('checking %s ' % key)
#     print ('-- identical in ds and ds_test : %s' % \
#            np.allclose(ds[key], ds_test[key], equal_nan=True))



### NOTES


        # for group in dsslim.groups:
    #     if group in list(filter(lambda k: '1' in k, list(dsslim.groups))):
    #         for variable in dsslim[group].variables:
    #             var_name = str(variable)+" "+group
    #             var_name = dsslim[group][variable]
    #             print(testds)
    # #dsone = dsslim[(list(filter(lambda k: '1' in k, list(dsslim.groups))))] # Filter one df to include all "1" sensors - doesn't work.
    # print(dsslim)
    # print(list(filter(lambda k: '2' in k, list(dsslim.groups)))) # Filter one df to include all "2" sensors
    
    
    #dsslim = nc4.delncattr('payload_1')
    #print(dsslim)
    #print(ds['payload_1/barometer_2'])
    #print(ds['payload_1/barometer_1/air_pressure_at_sea_level'][:]-ds['payload_1/barometer_2/air_pressure_at_sea_level'][:]) # Ranges from 40-100PA difference.
    #print(ds['payload_1/barometer_1/air_pressure_at_sea_level'][:]-ds['payload_2/barometer_1/air_pressure_at_sea_level'][:]) # Returns all zeros or nas
##
    #netcdf_flattener.flatten(ds, test)
    #print(test)    
    #break
    #dataset = xr.open_dataset(xr.backends.NetCDF4DataStore(ds))


# # read in the data
# ds1 = nc.open_data("file1.nc")
# ds2 = nc.open_data("file2.nc")

# # create a mask, where 1 is non-missing data, missing is missing
# # change var to whatever the variable is named
# ds2.assign(mask = lambda x: x.var == x.var, drop = True)
# # multiply the first file by the masked file
# ds1.multiply(ds2)
# # save the output file
# ds1.to_nc("file1_fixed.nc")

# exit()

# NetCDF files fom 2015-present use a nested group structure. As of 08/22 the netcdf4 package, rather than xarray is used to access groups and variables.

## From: https://stackoverflow.com/questions/31931483/programmatically-list-all-variables-of-a-netcdf-file-using-netcdf4-and-python
#os.chdir(workdir)

#dropvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.
#ds = netCDF4.Dataset('NDBC_46001_201501_D2_v00.nc')
#print(ds['/payload_1'].groups['barometer_2'])

#exit()

# def expand_var_list(var_list, group):
#     for var_key, _ in group.variables.items():
#         var_list.append(var_key)
    
#     for _, sub_group in group.groups.items():
#         expand_var_list(var_list, sub_group)

# all_vars = []
# with netCDF4.Dataset("NDBC_46001_201501_D2_v00.nc", "r") as root_group:
#     expand_var_list(all_vars, root_group)
# print(all_vars)

# exit()

# ## PLAY:
# os.chdir(workdir)
# dropvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.
# ds = Dataset('NDBC_46001_201501_D2_v00.nc')
# print(ds.groups)
# #ds = xr.open_dataset('NDBC_46001_201501_D2_v00.nc', drop_variables=dropvars) #, group = 'payload_1/anemometer_1'
# print(ds['/payload_1/anemometer_1/'])    
# exit()



#### NOTES

# For files after 2011
# # Each file has multiple payloads for the same data.
# e.g.
# print(ds['payload_1/anemometer_1/wind_speed'][:]-ds['payload_1/anemometer_2/wind_speed'][:]) # Range from -0.9 to 0.2
# print(ds['payload_1/anemometer_1/wind_speed'][:]-ds['payload_2/anemometer_1/wind_speed'][:]) # Returns all zeros.
# print(ds['payload_1/barometer_1/air_pressure_at_sea_level'][:]-ds['payload_1/barometer_2/air_pressure_at_sea_level'][:]) # Ranges from 40-100PA difference.
# print(ds['payload_1/barometer_1/air_pressure_at_sea_level'][:]-ds['payload_2/barometer_1/air_pressure_at_sea_level'][:]) # Returns all zeros or nas
# QC flags are also the same.
### Payload 1 appears to be rounded. Payload two appears to be the same data values but unrounded.
# Decision: take payload 1 and instrument 1 unless NA. If na, replace value with that of anemometer 2 / barometer 2, etc.

## one option: could do this in nco using -G: http://nco.sourceforge.net/nco.html#Group-Path-Editing
## but all the var names are the same within groups so this would not


#### TESTING NOTES
        # Quick and dirty list of all variables and their attributes
        #print(list(ds.attrs)) # Global attributes
        #print(list(ds.keys())) # Get all variables
        #for i in list(ds.keys()):
        #    print(i, ds[i].attrs) # Attributes for each variable
            
        


# For each station, join files together (by year or by all time?).
# Naming conventions for files change over time.
# MARITIME does not provide a regularly updated spreadsheet / .txt. file of station names (that I could find past 2020.)
# Thus, we manually strip station IDs from file names.

# File naming format changes over time as follows.
# Reference: https://www.ncei.noaa.gov/data/oceans/ndbc/cmanwx/
# 1970 - 06-1972
# xerb#_YYYYMM.nc

# 10-1972 - 1977
# eb##__YYYYMM.nc

# 1978 - 2010
# ?????_YYYYMM.nc

# 2011 - 
# NDBC_?????_YYYYMM_D*.nc

# Get list of file names
# files = os.listdir(workdir)
# files = list(filter(lambda f: f.endswith('.nc'), files))

# # Get list of IDs from file names
# ids = list()
# for file in files:
#     id = re.sub('_[0-9]{6,6}.*.nc', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
#     id = re.sub('NDBC_', "", id) # Remove leading NDBC
#     id = re.sub("_", "", id) # Remove any remaining underscores
#     if id not in ids:
#         ids.append(id)

# # Use ID to grab all files linked to station.
# for id in ids:
#     subfiles = list(filter(lambda f: id in f, files)) # Get all filenames containing unique ID.
#     # Next, open all these files and merge.
# #print(filestest)

# #filename = "32302_198901.nc"
# #file = os.path.join(workdir, filename) 

# #for stationname in stations:
# filename = "41001_200401.nc"
# file = os.path.join(workdir, filename) 
# def preprocess(ds):
#     return ds.drop_vars(dropvars)

# #ds = xr.open_mfdataset(file, combine = 'nested', compat="override", drop_variables = dropvars) # Does not work!.

# # Read in file and drop variables that aren't of interest.
# ds = xr.open_dataset(file, drop_variables = dropvars)
# print(ds)


# # Convert date-time to CF-compliant format (datetime object, UTC).
# ## Data is CF1.6 compliant, time in datetime format.

# # Convert station metadata to standard format.
# ## One column for lat (y), one column for lon (x).

# # Add dataset metadata in standard format.

# # Convert missing data to common format.
# ### Use NaN.


# Convert existing QA/QC flags to standard format.

# Deduplicate / merge files?
#print(ds.attrs)

# Step two: 
# In the long run use this code to open and merge multiple files (by station? by year?)

#file = os.path.join(workdir, "*.nc") 
# def preprocess(ds):
    # add any preprocessing steps in here.
    #return preprocessed ds
#ds = xr.open_mfdataset(file, preprocess = preprocess OR drop_variables = dropvars) # But how to keep dropped variables here?

#print(ds)
#vardrop = ['']
#ds = ds.drop_vars('sea_surface_temperature')
#print(ds)


# Step 2: join data for buoy by year.

#def clean0_maritime(workdir):

# %%
