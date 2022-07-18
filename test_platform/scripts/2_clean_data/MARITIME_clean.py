# Read in downloaded data and do first clean.
## Code written first by file, then rewrite to include merging mult. files (by time or by station?)

# Step one: read in, and drop variables that are not of interest.
filename = "32302_198901.nc"
file = os.path.join(workdir, filename) 

# Create list of variables to remove
dropvars = ['wave_wpm', 'wave_wpm_bnds', 'bottom_depth', 'magnetic_variation','sampling_rates_waves', 'sampling_duration_waves', 'total_intervals_waves', 'significant_wave_height', 'average_wave_period', 'dominant_wave_period','spectral_density_c']

# Read in file and drop variables that aren't of interest.
ds = xr.open_dataset(file, drop_variables = dropvars)
print(ds)

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
