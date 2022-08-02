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
import xarray as xr
from datetime import datetime
import re
import geopandas as gp
import numpy as np
import xarray as xr
from calc import _calc_relhumid, get_wecc_poly
import csv

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

try:
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
### note: add functionality here to only pull ~new/uncleaned~ files
os.chdir(workdir) # Change directory to where raw files saved.
files = os.listdir() # Gets list of files in directory to work with
files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names

# Get list of station IDs from filename and clean.
ids = list()
for file in files:
    id = re.sub('_[0-9]{6,6}.*.nc', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
    id = re.sub('NDBC_', "", id) # Remove leading NDBC
    id = re.sub("_", "", id) # Remove any remaining underscores
    if id not in ids:
        ids.append(id)


## Step 2: Read in datafile and drop variables that aren't of interest
# # Define list of column variables to drop.
# # This method is robust to variables getting added over time.
    
# Identify variables of interest and variable names.
keys = []
coords = [] # Coordinates need to be dropped seperately. Get list of coordinates.
dropvars = ['wave_wpm_bnds', 'wave_f40_bnds', 'wave_f40', 'wave_wa', 'wave_wa_bnds'] # Get rid of troublesome var.
    

for file in files:
    try:
        ds = xr.open_dataset(file, drop_variables=dropvars)
        key = list(ds.keys())
        coord = list(ds.coords)
        keys+=key
        coords+=coord
    except Exception as e:
        errors['File'].append(file)
        errors['Time'].append(end_api)
        errors['Error'].append(e)
        next

# Get list of all variables
mar_vars = list(set(keys)) # Get list of unique vars 
mar_vars = mar_vars + dropvars # add troublesome vars back in.

#print(mar_vars)

# Identify all variables to keep (this is a very conservative list, could prob cut more.)
vars_to_keep = ['solar_radiation_36', # Solar radiation
                'wind_speed', # Wind speed
                'speed_averaging_method', # Info on QAQC flags for wind.
                'air_temperature', # Air temp
                'wind_sampling_duration', # Duration of averaging for wind speed.
                'air_pressure_at_sea_level', # Air pressure
                'time10_bnds', 'timem_bnds', 'time_bnds', # Time bounds for later (could prob only keep one)
                'dew_point_temperature',  # Dew point
                'precipitation', # Precipitation
                'wind_direction', # Wind direction
                'anemometer_height'] # For QA/QC


# Identify list of coordinates to drop
mar_coords = list(set(coords))
coords_to_remove = ['depth', 'timem', 'time_wpm_20', 'time10', 'wave_wpm']

## Vars I removed we may want:
## continuous_wind_speed: with duration specified in comment (generally 10min). optional, while wind_speed always reported.
## wind gust    

# Remove these variables from the list of all variables to get drop list.
dropvars = np.setdiff1d(mar_vars, vars_to_keep)

#print(dropvars)

for file in files:
    try:
        ## Step 2: Read in datafile and drop variables that aren't of interest
        ds = xr.open_dataset(file, drop_variables=dropvars)

        # Quality control
        ## Before we clean, skip any files that don't have data in them. (save to error list?)
        if not list(ds.data_vars):
            errors['File'].append(file)
            errors['Time'].append(end_api)
            errors['Error'].append("No variables recorded in file")
            next

        ## Before we clean, require lat/lon coordinates and filter by location in WECC.
        
        if 'lat' or 'lon' in ds.keys():
            #Filter to only keep lat and lon in WECC.
            if float(max(ds['lat'])) > latmax or float(min(ds['lat'])) < latmin or float(max(ds['lon'])) > lonmax or float(min(ds['lon'])) < lonmin:
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(float(max(ds['lat'])), float(min(ds['lat']))))
                next 
            elif float(ds.attrs['geospatial_lat_min']) < latmin or float(ds.attrs['geospatial_lat_max']) > latmax or float(ds.attrs['geospatial_lon_min']) < lonmin or float(ds.attrs['geospatial_lon_max']) > lonmax:
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append("File not in WECC. Lat: {} Lon: {}".format(float(ds.attrs['geospatial_lat_min']), float(ds.attrs['geospatial_lon_min'])))
                next
            else:
                print("File {} in WECC!".format(file))
        else:
            errors['File'].append(file)
            errors['Time'].append(end_api)
            errors['Error'].append("No geospatial coordinates found in file")
            next
        
        # Remove unwanted coordinates.
        ctr = [x for x in coords_to_remove if x in ds.keys()] # Filter coords to remove list by coordinates found in dataset.
        ds = ds.drop_dims(ctr) # Drop coordinates.


        ## Step 3: Convert station metadata to standard format -- CF compliance

        # 3.1 Generate unique ID across entire dataset by combining network name and ID.
        ds = ds.assign_attrs(station_id = "MARITIME_"+ds.attrs["id"])
        # Rename original ID column
        ds.attrs['original_id'] = ds.attrs.pop('id')
        
        # 3.2 Check lat / lon - TO DO 
        # 3.3 Check elev. - TO DO

        ## Step 4: Convert dataset metadata in standard format -- CF compliance
        ### SKIPPING FOR NOW.
        
        ## Step 5: Convert missing data to common format -- CF compliance
        ## Use NaN
        #for var in ds.variables: # Get range of variable values to check for non-NAN NA values. How??
        
        # Hacky way to get count of NAs in a column.
        #print(ds['wind_speed'].sum(skipna=False) - ds['wind_speed'].sum(skipna=True))

        #break
        
        # Step 6: Standardize variables and variable metadata.
        # If not observed, calculate derived primary variables
        # note!!!! not all stations have all columns. Add try / exceptions here or other code to help manage breaks!!
        #** indicates primary approach

        # In following order:
        # clean_vars = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']

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
        # Still to do.
        if 'solar_radiation_36' in ds.keys():
            try:
                print(ds['solar_radiation_36'].attrs['long_name']) #solar_radiation_wavelength_less_than_3_6_um
                print(ds['solar_radiation_36'].attrs)
            except:
                next

        # # sfcWind : wind speed
        ### In maritime we have:
        ## wind_speed: average wind speed (with duration of averaging specified in wind_sampling_duration)
        ## speed_averaging_method: contains QAQC flags.
        if "wind_speed" in ds.keys(): # If variable already exists, rename. Units already in m s-1.
            ds = ds.rename({'wind_speed': 'sfcWind'})
        
        # # sfcWind_dir: wind direction
        if "wind_direction" in ds.keys(): # If variable already exists, rename.
            ds = ds.rename({'wind_direction': 'sfcWind_dir'}) # Already in degrees.
        
        # # Reorder variables
        #print(ds.keys())
        # ds = ds['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
        #print(ds)
        # ## Step 6: Tracks existing QA/QC flags to standard format
        ### again, maritime should have flags but i don't see them currently....
        # old_flag = []
        # histobs_flag = []

        # Reorder variables

    except Exception as e:
        errors['File'].append(file)
        errors['Time'].append(end_api)
        errors['Error'].append(e)
        next
        
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


#### NOTES

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
