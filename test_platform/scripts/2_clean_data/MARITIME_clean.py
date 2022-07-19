# Import library
import os
import xarray as xr
from datetime import datetime
import re
import geopandas as gp

# Read in downloaded data and do first clean.
## Code written first by file, then rewrite to include merging mult. files (by time or by station?)

# Get current time
timestamp = datetime.now()

# Set envr vars.
workdir = "/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/Maritime/"
wecc_terr = '/home/ella/Desktop/Eagle Rock/Historical Data Platform /historical-obs-platform/test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp'
wecc_mar = '/home/ella/Desktop/Eagle Rock/Historical Data Platform /historical-obs-platform/test_platform/data/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp' 

# Create list of variables to remove
dropvars = ['wave_wpm', 'wave_wpm_bnds', 'bottom_depth', 'magnetic_variation','sampling_rates_waves', 'sampling_duration_waves', 'total_intervals_waves', 'significant_wave_height', 'average_wave_period', 'dominant_wave_period','spectral_density_c',
                'end_of_wave_data_acquisition_k', 'c11_k']

# TEST: Only keep stations in WECC
def get_wecc_poly(terrpath, marpath):
    ## get bbox of WECC to use to filter stations against
    ## Read in terrestrial WECC shapefile.
    t = gp.read_file(terrpath)
    ## Read in marine WECC shapefile.
    m = gp.read_file(marpath)
    ## Combine polygons and get bounding box of union.
    bbox = t.union(m).bounds
    return t,m, bbox

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
files = os.listdir(workdir)
files = list(filter(lambda f: f.endswith('.nc'), files))

# Get list of IDs from file names
ids = list()
for file in files:
    id = re.sub('_[0-9]{6,6}.*.nc', "",file) # Remove all YYYYMM (and trailing metadata/file extension)
    id = re.sub('NDBC_', "", id) # Remove leading NDBC
    id = re.sub("_", "", id) # Remove any remaining underscores
    if id not in ids:
        ids.append(id)

# Use ID to grab all files linked to station.
for id in ids:
    subfiles = list(filter(lambda f: id in f, files)) # Get all filenames containing unique ID.
    # Next, open all these files and merge.
#print(filestest)

#filename = "32302_198901.nc"
#file = os.path.join(workdir, filename) 

#for stationname in stations:
filename = "41001_200401.nc"
file = os.path.join(workdir, filename) 
def preprocess(ds):
    return ds.drop_vars(dropvars)

#ds = xr.open_mfdataset(file, combine = 'nested', compat="override", drop_variables = dropvars) # Does not work!.

# Read in file and drop variables that aren't of interest.
ds = xr.open_dataset(file, drop_variables = dropvars)
print(ds)


# Convert date-time to CF-compliant format (datetime object, UTC).
## Data is CF1.6 compliant, time in datetime format.

# Convert station metadata to standard format.
## One column for lat (y), one column for lon (x).

# Add dataset metadata in standard format.

# Convert missing data to common format.
### Use NaN.


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
