# Import library
import os
import xarray as xr
from datetime import datetime

# Read in downloaded data and do first clean.
## Code written first by file, then rewrite to include merging mult. files (by time or by station?)

# Get current time
timestamp = datetime.now()

# Set envr vars.
workdir = "/home/ella/Desktop/Eagle-Rock/Historical-Data-Platform/Maritime/"
file = os.path.join(workdir, filename) 

# Create list of variables to remove
dropvars = ['wave_wpm', 'wave_wpm_bnds', 'bottom_depth', 'magnetic_variation','sampling_rates_waves', 'sampling_duration_waves', 'total_intervals_waves', 'significant_wave_height', 'average_wave_period', 'dominant_wave_period','spectral_density_c']

# Only keep stations in WECC

# For each station, join files together (by year or by all time?).
filename = "32302_198901.nc"
files = os.listdir(workdir)
files = list(filter(lambda f: f.endswith('.nc'), files))
stations = 
#for stationname in stations:
#    ds = xr.open_mfdataset('{}_*.nc'.format(stationname), combine = 'by_coord', concat_dim = 'time')

# Read in file and drop variables that aren't of interest.
ds = xr.open_dataset(file, drop_variables = dropvars)

# Convert date-time to CF-compliant format (datetime object, UTC).
## Data is CF1.6 compliant, time in datetime format.

# Convert station metadata to standard format.
## One column for lat (y), one column for lon (x).

# Add dataset metadata in standard format.

# Convert missing data to common format.
### Use NaN.


# Convert existing QA/QC flags to standard format.

# Deduplicate / merge files?
print(ds.attrs)

# Step two: 
# In the long run use this code to open and merge multiple files (by station? by year?)

#file = os.path.join(workdir, "*.nc") 
# def preprocess(ds):
    # add any preprocessing steps in here.
    #return preprocessed ds
#ds = xr.open_mfdataset(file, preprocess = preprocess OR drop_variables = dropvars) # But how to keep dropped variables here?

print(ds)
#vardrop = ['']
#ds = ds.drop_vars('sea_surface_temperature')
#print(ds)


# Step 2: join data for buoy by year.

#def clean0_maritime(workdir):
