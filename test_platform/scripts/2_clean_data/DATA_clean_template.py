## Historical Observations Platform
## Template Script for Stage: CLEAN DATA

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries
import os
import xarray as xr
from datetime import datetime
import re
import pandas

# Set envr variables
workdir = "/path/to/working/directory/where/data/is/for/testing/"


# Step 1: First pass of clean
timestamp = datetime.now()

# Create list of variables to remove
dropvars = ['var1', 'var2', 'var3'] # will depend on each datasource and what it provides


# Step 2: Read in datafile and drop variables that aren't of interest
files = os.listdir(workdir) # Gets list of files in directory to work with
files = list(filter(lambda f: f.endswith(".nc"), files)) # may need to change extension as needed

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

# Step 3: Convert station metadata to standard format


# Step 4: Convert dataset metadata in standard format
## Unit check? Conversions needed?


# Step 5: Convert missing data to common format
## Use NaN


# Step 6: Convert existing QA/QC flags to standard format
old_flag = []
histobs_flag = []


# Step 7: Open datafile and merge files by station -- in this script or the next?


# What needs to be returned here?
## data should return cleaned (and organized?) dataframe per source
## in terms of organization, should determine variable order?
