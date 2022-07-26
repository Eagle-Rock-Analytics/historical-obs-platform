"""
This script is a template structure for data cleaning for a variety of data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Read through variables, and calculates derived priority variables if not observed
(2) Drops unnecessary variables
(3) Converts station metadata to standard format, with unique identifier
(4) Converts data metadata to standard format, and converts units into standard units if not provided in standard units.
(5) Converts missing data to standard format
(6) Tracks existing qa/qc flag for review
(7) Merge files by station, and outputs cleaned variables as a single .nc file for an individual network.
Inputs: Raw data for an individual network
Outputs: Cleaned data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime
import re
import pandas
import numpy as np

# Set envr variables
workdir = "/path/to/working/directory/where/data/is/for/testing/"


## Step 1: Read in files and first pass of clean and organize
timestamp = datetime.now()

## -----------------------------------------------------------------------------------------------------------------
## This may happen via various methods depending on how the data is stored
files = os.listdir(workdir) # Gets list of files in directory to work with

## netCDF
files = list(filter(lambda f: f.endswith(".nc"), files))
# to be continued

## CSV
files = list(filter(lambda f: f.endswith(".csv"), files))
# to be continued

## Text files
files = list(filter(lambda f: f.endswith(".txt"), files))
# to be continued


# -------------------------------------------------------------------------------------------------------------------

# dimension/variable order
# For example, air_temperature would have (time x lat x lon x elev x data) dimensions
clean_dims = ['time', 'latidue', 'longitude', 'elevation']
clean_vars = ['ps', 'tas', 'tdps', 'pr', 'hurs', 'rsds', 'sfcWind', 'sfcWind_dir']
## match these back to their cf compliant names

# If not observed, calculate derived primary variables
# ** indicates primary approach

# dew point temperature calculation (necessary input vars: requires at least 2 of three - air temp + relative humidity + vapor pressure)
def _calc_dewpointtemp(tas, hurs, e):
    es = 0.611 * exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure
    e = (es * hurs)/100                        # calculates vapor pressure, IF NOT ALREADY OBSERVED -- will need ifelse statement
    tdps = ((1/273) - 0.0001844 * ln(e/0.611))^-1   # calculates dew point temperature, units = K

return tdps

# relative humidity calculation (necessary input vars: air temp + dew point**, air temp + vapor pressure, air pressure + vapor pressure)
def _calc_relhumid(tas, tdps):
    es = 0.611 * exp(5423 * ((1/273) - (1/tas)))   # calculates saturation vapor pressure using air temp
    e = 0.611 * exp(5423 * ((1/273) - (1/tdps)))   # calculates vapor pressure using dew point temp
    hurs = 100 * (e/es)

return hurs

# wind speed (necessary input vars: u and v components)
def _calc_windmag(u10, v10):
    sfcWind = np.sqrt((u10)^2  + (v10)^2)   # calculates wind magnitude, units = ms-1

return sfcWind

# wind direction (necessary input vars: u and v components)
def _calc_winddir(u10, v10):
    pass        # this is a complicated calculation -- looking for options

return sfcWind_dir


## Step 2: Read in datafile and drop variables that aren't of interest
dropvars = ['var1', 'var2', 'var3'] # will depend on each datasource and what it provides

### QUESTION FOR GRACE/OWEN: Do we keep the moisture variables if they were needed for RH/TD calculations or drop?

# ----------------------------------------------------------------------------------------------------------
# This is based on if files are listed by station
ids = list()
for file in files:
    id = re.sub(".*.nc", "", file) # removes all YYYYMM (and trailing metadata/file extension, if .nc)
    id = re.sub('XXXX_', "", id) # removes leading ndbc (example)
    id = re.sub("_", "", id) # removes remaining underscores
    if id not in ids:
        ids.append(id)

# This is based on if files are listed by date (CWOP_SR is an example)
# TBD


# ----------------------------------------------------------------------------------------------------------

file = os.path.join(workdir, filename)
def preprocess(ds):
    return ds.drop_vars(dropvars) # Drops variables that aren't of interest

ds = xr.open_dataset(file, drop_variables = dropvars)
print(ds)


# TO DO (potential?)
# consistent date-time format?
# a flag/grouping for hourly/sub-hourly

### QUESTION FOR GRACE/OWEN: Do we keep the 'original' and 'converted' variables? Only 'converted' variables would be read by DATA_qaqc.py

## Step 3: Convert station metadata to standard format -- CF compliance
# Generate unique ID across entire dataset
# new_ID = network_name + network_id ?


## Step 4: Convert dataset metadata in standard format -- CF compliance
## Unit check? Conversions needed?


## Step 5: Convert missing data to common format -- CF compliance
## Use NaN


## Step 6: Tracks existing QA/QC flags to standard format
old_flag = []
histobs_flag = []


## Step 7: Open datafile and merge files by station

# GOAL: outputs cleaned .nc file with dimensions + variables in correct order with standardized metadata
