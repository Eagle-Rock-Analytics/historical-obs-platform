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
#import calc # import calc.py
import numpy as np
import xarray as xr

# Set envr variables
workdir = "/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/Maritime/"
#workdir = "test_platform/data/1_raw_wx/MARITIME/"
wecc_terr = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = 'test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 

## Step 1: Read in files and first pass of clean and organize
timestamp = datetime.now()

# Set target dimension and variable order
# For example, air_temperature would have (time x lat x lon x elev x data) dimensions
clean_dims = ['time', 'latitude', 'longitude', 'elevation']
clean_vars = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']

# Read in files.
### note: add functionality here to only pull ~new/uncleaned~ files
os.chdir(workdir) # Change directory to where files saved.
files = os.listdir(workdir) # Gets list of files in directory to work with
files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names

# Get list of station IDs from filename and clean.
ids = list()
for file in files:
    id = re.sub('_[0-9]{6,6}.*.nc', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
    id = re.sub('NDBC_', "", id) # Remove leading NDBC
    id = re.sub("_", "", id) # Remove any remaining underscores
    if id not in ids:
        ids.append(id)

#print(ids)
# for file in files:
#     try:
#         ## Step 2: Read in datafile and drop variables that aren't of interest
#         # Define list of column variables to drop.
#         dropvars = ['wave_wpm', 'wave_wpm_bnds', 'bottom_depth', 'magnetic_variation','sampling_rates_waves', 'sampling_duration_waves', 'total_intervals_waves', 'significant_wave_height', 'average_wave_period', 'dominant_wave_period','spectral_density_c',
#                         'end_of_wave_data_acquisition_k', 'c11_k']

#         ds = xr.open_dataset(file, drop_variables=dropvars)
        
#         ## Step 3: Convert station metadata to standard format -- CF compliance
#         # Generate unique ID across entire dataset by combining network name and ID.
#         ds = ds.assign_attrs(station_id = "MARITIME_"+ds.attrs["id"])
#         # Rename original ID column
#         ds.attrs['original_id'] = ds.attrs.pop('id')

#         # Check lat / lon - TO DO 
#         # Check elev. - TO DO

#         #### TESTING NOTES
#         # Quick and dirty list of all variables and their attributes
#         print(list(ds.attrs)) # Global attributes
#         print(list(ds.keys())) # Variables
#         for i in list(ds.keys()):
#             print(i, ds[i].attrs) # Attributes for each variable
            
#         ## Step 4: Convert dataset metadata in standard format -- CF compliance
#         ## Unit check? Conversions needed?
#         ### SKIPPING FOR NOW.
        
#         ## Step 5: Convert missing data to common format -- CF compliance
#         ## Use NaN
#         #for var in ds.variables: # Get range of variable values to check for non-NAN NA values. How??
        
#         # Hacky way to get count of NAs in a column.
#         print(ds['wind_speed'].sum(skipna=False) - ds['wind_speed'].sum(skipna=True))

#         break
        
#         ## Step 6: If not observed, calculate derived primary variables
#         # ** indicates primary approach

#         # #dew point temperature calculation (necessary input vars: requires at least 2 of three - air temp + relative humidity + vapor pressure)
#         # ## MARITIME network should have dewpoint but does not appear to????
#         # tdps = calc._calc_dewpointtemp(tas, hurs, e)

#         # # relative humidity calculation (necessary input vars: air temp + dew point**, air temp + vapor pressure, air pressure + vapor pressure)
#         # hurs = calc._calc_relhumid(tas, tdps)

#         # # wind speed (necessary input vars: u and v components)
#         # ## Maritime already has wind speed calculated.
#         # sfcWind = calc._calc_windmag(u10, v10)

#         # # wind direction (necessary input vars: u and v components)
#         # ## Maritime already has wind speed calculated.
#         # sfcWind_dir = calc._calc_winddir(u10, v10)

#         # ## Step 6: Tracks existing QA/QC flags to standard format
#         ### again, maritime should have flags but i don't see them currently....
#         # old_flag = []
#         # histobs_flag = []

#     except Exception as e:
#         print("Error processing file {}: {}".format(file, e)) # To do: keep track of these.
    
    
    ## Step 7: Open datafile and merge files by station

    # For each station ID, read in all files.
    
    ## Save each cleaned file in the new working directory. 
    # Then, as last line of script, scan to see if there is more than one file per station id
    # and read in and merge cleaned files together if so.
    
    ## Change WD

    # For each ID, get list of files with ID in name.
    for id in ids:
        file_sub = list(filter(lambda k: id in k, files)) # Add, "MARITIME" in k.
        print(id, file_sub)
        id_comb = xr.open_mfdataset(file_sub,combine = 'by_coords', concat_dim="time") # Theoretically works.
        id.to_netcdf('MARITIME_{}.nc'.format(id)) # Write to one netcdf file
        # Clean up file_sub files if successful.


#### NOTES
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
