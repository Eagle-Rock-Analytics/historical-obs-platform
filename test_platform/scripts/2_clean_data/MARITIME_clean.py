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

# Step 0: Environment set-up
# Import libraries
import os
from regex import D
import xarray as xr
from datetime import datetime
import re
import geopandas as gp
import numpy as np
import xarray as xr
from calc import _calc_relhumid, get_wecc_poly
import csv
import netCDF4 as nc4

### TO DO:
# Keep anemometer height
# Keep anemometer manufacturer + part number

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

# Testing
os.chdir(workdir) # Change directory to where raw files saved.
files = os.listdir() # Gets list of files in directory to work with
files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names
files_old = [i for i in files if any(i for j in list(range(1980, 2015)) if str(j) in i)] # Get all file names which include any year from 1980-2014.
files_new = [i for i in files if i not in files_old] # Get all other file names.
dropvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.

for file in files:
    if file in files_new: # If file newer than 2015, process to remove group structure.
        # Get global attributes from main file
        ds = xr.open_dataset(file, drop_variables=dropvars)
        print(ds.attrs['geospatial_lat_max'], ds.attrs['geospatial_lat_min'])

# Step 0: Get list of all variables and generate list of variables to remove from dataset.
## FUNCTION: Generate list of variables, drop variables from collection of files.
# Input: workdir and savedir. Output: 2 csvs: 1) of removed variables and 2) error file. Also returns list of removed variables.
# This function natively generates a list of all variables (including from groups), compares them against a list of variables to keep
# and prints a list of removed variables to a csv in the save directory. This should be run occasionally but not daily.

def get_vars(workdir, savedir):
    homedir = os.getcwd()
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling.
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download

    os.chdir(workdir) # Change directory to where raw files saved.
    files = os.listdir() # Gets list of files in directory to work with
    files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names
    files_old = [i for i in files if any(i for j in list(range(1980, 2015)) if str(j) in i)] # Get all file names which include any year from 1980-2014.
    files_new = [i for i in files if i not in files_old] # Get all other file names.
    keys = []
    coords = [] # Coordinates need to be dropped seperately. Get list of coordinates.
    dropvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.

    for file in files:
        if file in files_new: # If file newer than 2015, process to remove group structure.
            try:
                # Get global attributes from main file
                ds = xr.open_dataset(file, drop_variables=dropvars)
                # Get list of sensors from global attributes
                if 'sensor_suite' in ds.attrs:
                    sensor_list = ds.attrs['sensor_suite'].split(',')
                    prim_sensors = [k for k in sensor_list if 'payload_1' in k] # Get data from primary payload
                    prim_sensors = [k for k in prim_sensors if k.endswith("1")] # Get data from primary sensor - may change depending on email response.
                else: # Otherwise, try all possible sensors manually.
                    prim_sensors = ['/payload_1/anemometer_1', '/payload_1/humidity_sensor_1', '/payload_1/ocean_temperature_sensor_1', '/payload_1/air_temperature_sensor_1', '/payload_1/barometer_1', '/payload_1/gps_1', '/payload_1/wave_sensor_1']
                for path in prim_sensors: # For each sensor
                    try:
                        ds_temp = xr.open_dataset(file, drop_variables=dropvars, group = path) # Open sensor data
                    except Exception as e:
                        #print(file, e) # To revisit: not reporting error here because some groups don't have sensor list (and so we have to try the full suite.)
                        continue # Continue to next branch
                    ds = xr.merge([ds, ds_temp]) # Move sensor data to main ds
                    del(ds_temp)
            except Exception as e:
                #print(file, e) # For testing. Remove in final run.
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue
        elif file in files_old:
            try:
                ds = xr.open_dataset(file, drop_variables=dropvars)
            except Exception as e:
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)
                continue
            # Flag to return to: is there external data on sensor types we want to incorporate here?
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
    mar_vars = mar_vars + dropvars # add troublesome vars back in.

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
                    'time10_bnds', 'timem_bnds', 'time_bnds', # Time bounds for later (could prob only keep one)
                    'dew_point_temperature',  # Dew point
                    'precipitation', # Precipitation
                    'wind_direction', # Wind direction
                    'relative_humidity',
                    'anemometer_height', # For QA/QC (TO DO: var or attribute?)
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
    
    return dropvars

#test = get_vars(workdir, savedir)
#print(test)

## FUNCTION: Clean MARITIME data.
# Input: workdir, savedir.
# Note that this function searches for a csv called removevars.csv in the savedir folder. If file doesn't exist, script will call get_vars function.
def clean_maritime(workdir, savedir):
    
    homedir = os.getcwd() # Store home directory.
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
    # ids = list()
    # for file in files:
    #     id = re.sub('_[0-9]{6,6}.*.nc', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
    #     id = re.sub('NDBC_', "", id) # Remove leading NDBC
    #     id = re.sub("_", "", id) # Remove any remaining underscores
    #     if id not in ids:
    #         ids.append(id)

    # ## Step 2: Read in datafile and drop variables that aren't of interest
    # # # Read in list of column variables to drop. If no list, call get_vars().
    # # # This method is robust to variables getting added over time.
    
    dropvars = []
    path = savedir+"/removedvars.csv" 
    if os.path.exists(path): # If path to csv exists, attempt to read it. 
        with open(path,'r') as f: 
                rows = list(csv.DictReader(f, delimiter=','))
                for row in rows:
                    dropvars.append(row['Variable'])
    else: # Otherwise, call get_vars function.
        dropvars = get_vars(workdir, savedir)
    
    # Identify list of coordinates to drop
    coords_to_remove = ['depth', 'timem', 'time_wpm_20', 'time10', 'wave_wpm'] # Hardcoded for now, could get moved. These are a list of coords that throw errors consistently.
    
    # # Split files into those written before 2015 and those written after.
    files_old = [i for i in files if any(i for j in list(range(1980, 2015)) if str(j) in i)] # Get all file names which include any year from 1980-2014.
    files_new = [i for i in files if i not in files_old] # Get all other file names.
    
    # Read in data files.
    for file in files:
        try:
            if file in files_new: # If file newer than 2015, process to remove group structure.
                # Get global attributes from main file
                ds = xr.open_dataset(file, drop_variables=dropvars)
                
                ## Before we clean, require lat/lon coordinates and filter by location in WECC.
                ### STOPPED HERE ON 08/09: lat/lon recorded differently between new and old file types.
                if 'geospatial_lat_min' and 'geospatial_lat_max' and 'geospatial_lon_min' and 'geospatial_lat_max' in ds.attrs:
                    if float(ds.attrs['geospatial_lat_min']) < latmin or float(ds.attrs['geospatial_lat_max']) > latmax or float(ds.attrs['geospatial_lon_min']) < lonmin or float(ds.attrs['geospatial_lon_max']) > lonmax:
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(float(ds.attrs['geospatial_lat_min']), float(ds.attrs['geospatial_lon_min'])))
                        continue
                    else:
                        print("File {} in WECC!".format(file)) # Just for testing, remove for final.
                else:
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append("No geospatial coordinates found in file")
                    continue

                # Get list of sensors from global attributes
                if 'sensor_suite' in ds.attrs:
                    sensor_list = ds.attrs['sensor_suite'].split(',')
                    prim_sensors = [k for k in sensor_list if 'payload_1' in k] # Get data from primary payload
                    prim_sensors = [k for k in prim_sensors if k.endswith("1")] # Get data from primary sensor - may change depending on email response.
                else: # Otherwise, try all possible sensors manually.
                    prim_sensors = ['/payload_1/anemometer_1', '/payload_1/humidity_sensor_1', '/payload_1/ocean_temperature_sensor_1', '/payload_1/air_temperature_sensor_1', '/payload_1/barometer_1', '/payload_1/gps_1', '/payload_1/wave_sensor_1']
                for path in prim_sensors: # For each sensor
                    try:
                        ds_temp = xr.open_dataset(file, drop_variables=dropvars, group = path) # Open sensor data
                    except Exception as e:
                        #print(file, e) # To revisit: not reporting error here because some groups don't have sensor list (and so we have to try the full suite.)
                        continue # Continue to next branch
                    ds = xr.merge([ds, ds_temp]) # Move sensor data to main ds

                    # Add sensor attributes to main ds.
                    if 'manufacturer' and 'part_number' in ds_temp.attrs:
                        sensor_name = str(ds_temp.attrs['manufacturer'])+" "+str(ds_temp.attrs['part_number'])# Get sensor data as string from dataframe.
                        sensor_type = re.sub('[^a-zA-Z]+', '', path.rsplit('/', 1)[-1]) # Remove all spaces and numbers
                    elif 'manufacturer' in ds_temp.attrs:
                        sensor_name = str(ds_temp.attrs['manufacturer'])# Get sensor data as string from dataframe.
                        sensor_type = re.sub('[^a-zA-Z]+', '', path.rsplit('/', 1)[-1])
                    else:
                        sensor_type = re.sub('[^a-zA-Z]+', '', path.rsplit('/', 1)[-1])
                        sensor_name = "Unknown"
                    ds.attrs[sensor_type] = sensor_name # Assign sensor manufacturing info as attribute.

                    install_name = sensor_type+"_install_date"  
                    if 'install date' in ds_temp.attrs:
                        ds.attrs[install_name] = ds_temp.attrs['install_date'] # Assign install date as attribute
                    else:
                        ds.attrs[install_name] = "Unknown"

                    if 'height_of_instrument' in ds_temp.attrs: # TO RESOLVE: MAKE ANEMOMETER HEIGHT A VAR (as it is in pre 2015 data) OR AN ATTR?
                        sensor_height = sensor_type+"_height"  
                        ds.attrs[sensor_height] = ds_temp.attrs['height_of_instrument'] # Assign install date as attribute
                    del(ds_temp)
                
                # If after merging, ds is empty, record error and skip file.
                if not list(ds.data_vars):
                    errors['File'].append(file)
                    errors['Time'].append(end_api)
                    errors['Error'].append("No variables recorded in file")
                    continue
            
            # For older files, simply read in.
            if file in files_old:
                ds = xr.open_dataset(file, drop_variables=dropvars)
        
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

### STOPPED HERE ON 08/09. Above needs testing as well.

            ## Step 3: Convert station metadata to standard format

            # 3.1 Generate unique ID across entire dataset by combining network name and ID.
            ds = ds.assign_attrs(station_id = "MARITIME_"+ds.attrs["id"])
            # Rename original ID column
            ds.attrs['original_id'] = ds.attrs.pop('id')
            
            # 3.2 Check lat / lon - original dataset is CF-compliant, data should be in standard format.
            # 3.3 Check elev. - original dataset is CF-compliant, data should be in standard format.

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
            #for var in ds.variables: # Get range of variable values to check for non-NAN NA values. How??
            ### THIS IS EXTREMELY UGLY. FIND BETTER METHOD!
            # nas = ['-999', '-9999','n/a', 'NA', 'none', "NaN", 'NAN', 'nan', "1989-05-01T00:00:00.000000000"]
            # stamp = "1989-05-01T00:00:00.000000000"
            # na_list = []
            # #print(ds.sel(time = stamp))
            # for item in nas:
            #     try:
            #         #print(item)
            #         print(ds.sel(time = item))
            #         na_list += item
            #         print(ds.sel(wind_speed = item))
            #         na_list += item
            #         print(ds.sel(air_temperature = item))
            #         na_list += item
            #         # If any of these return errors, then the value has not been found for that variable.
            #     except:
            #         continue
            # print(na_list)
            
            # ## Use NaN
            # # ds_masked = ds.where(ds['var'] != -9999.) # Replace -9999 here with any identified NA values. 
            
            # Step 6: Standardize variables and variable metadata.
            # If not observed, calculate derived primary variables if possible.
            #** indicates primary approach

            # ps: surface air pressure
            if "air_pressure_at_sea_level" in ds.keys():
                ds['ps'] = ds["air_pressure_at_sea_level"]*100 # Convert from mb to Pa
                ds['ps'].attrs['units'] = "Pascal" # Update units
                ds = ds.rename({'air_pressure_at_sea_level': 'ps_raw'}) # Convert from mb to Pa
                
            # tas : air surface temperature
            if "air_temperature" in ds.keys():
                ds['tas'] = ds["air_temperature"]+273.15 # Convert from C to K
                ds['tas'].attrs['units'] = "degree_Kelvin" # Update units
                ds = ds.rename({'air_temperature': 'tas_raw'})

            # tdps: dew point temperature
            # dew point temperature calculation (necessary input vars: requires at least 2 of three - air temp + relative humidity + vapor pressure)
            # Only more recent stations have dewpoint temp (dew_point_temperature)
            # No stations have vapor pressure or RH, so no calculations possible.
            if "dew_point_temperature" in ds.keys(): # If variable already exists, rename.
                ds['tdps'] = ds["dew_point_temperature"]+273.15 # Convert from C to K
                ds['tdps'].attrs['units'] = "degree_Kelvin"
                ds = ds.rename({'dew_point_temperature': 'tdps_raw'})
            
            # # pr: precipitation - DO CONVERSIONS FROM RAINFALL TO RAINFALL RATE DURING NEXT STAGE.
            if "precipitation" in ds.keys():
                ds = ds.rename({'precipitation': 'p_raw'})
                #ds['p'] = ds['precipitation'] # Convert mm to kg m-2 s-1 - do in next stage.

            # hurs: relative humidity
            # relative humidity calculation (necessary input vars: air temp + dew point**, air temp + vapor pressure, air pressure + vapor pressure)
            ## Vapor pressure not a variable measured by this network.
            if 'tas' and 'tdps' in ds.keys():
                ds['hurs'] = _calc_relhumid(ds['tas'], ds['tdps'])

            # # rsds: surface_downwelling_shortwave_flux_in_air (solar radiation)
            #### Still to do.
            if 'solar_radiation_36' in ds.keys():
                try:
                    print(ds['solar_radiation_36'].attrs['long_name']) #solar_radiation_wavelength_less_than_3_6_um
                    print(ds['solar_radiation_36'].attrs)
                except:
                    continue

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

            # # Reorder variables
            # In following order:
            desired_order = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
            desired_order = [i for i in desired_order if i in list(ds.keys())] # Only keep vars which are in ds.
            rest_of_vars = [i for i in list(ds.keys()) if i not in desired_order] # Retain rest of variables at the bottom.
            new_index = desired_order + rest_of_vars
            ds = ds[new_index]


        except Exception as e:
            print(e) # For testing only.
            errors['File'].append(file)
            errors['Time'].append(end_api)
            errors['Error'].append(e)
            continue
        
# Write cleaned files to netcdf - in progress.
os.chdir(homedir)
os.chdir(savedir) # Change directory to where files saved.

# Write errors to csv
filepath = "errors_cwop_{}.csv".format(end_api) # Set path to save error file.
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

    ## Step 7: Open datafile and merge files by station

    # For each station ID, read in all files.
    
    ## Save each cleaned file in the new working directory. 
    # Then, as last line of script, scan to see if there is more than one file per station id
    # and read in and merge cleaned files together if so.
    
    ## Change WD

    # For each ID, get list of files with ID in name.
    # for id in ids:
    #     file_sub = list(filter(lambda k: id in k, files)) # Add, "MARITIME" in k.
    #     print(id, file_sub)
    #     id_comb = xr.open_mfdataset(file_sub,combine = 'by_coords', concat_dim="time") # Theoretically works.
    #     id.to_netcdf('MARITIME_{}.nc'.format(id)) # Write to one netcdf file
    #     # Clean up file_sub files if successful.

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

# For files after 2015
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
