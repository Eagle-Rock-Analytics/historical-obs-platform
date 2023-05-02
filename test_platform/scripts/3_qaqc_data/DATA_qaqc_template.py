"""
This script is a template structure for data qa/qc protocols for a variety of data sources for
ingestion into the Historical Observations Platform.
Approach:
(1) Remove duplicate stations
(2) Handle variables that report at different intervals and/or change frequency over time (convert to hourly?)
(3) QA/QC testing, including consistency checks, gaps, checks against climatological distributions, and cross variable checks.
(4) Case study analysis for accuracy
Inputs: Cleaned data for an individual network
Outputs: QA/QC-processed data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries
import os
import datetime
import pandas as pd
import xarray as xr
import boto3
import s3fs
from io import BytesIO, StringIO


## Import qaqc stage calc functions
try:
    from calc_qaqc import *
except:
    print("Error importing calc_qaqc.py")

## Set up directory to save files temporarily, if it doesn't already exist.
try:
    os.mkdir('temp') # Make the directory to save data in. Except used to pass through code if folder already exists.
except:
    pass

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"



## Step 0: Development standards
## Calculate false positive rate (our goal is <2%), evaluate data flagged by original network

# Suggest we use ASOS/AWOS; CWOP; MARITIME as test networks for development
# From each network station list, randomly sample 100 stations

# Write function to do a check

# For distributional checks or timeseries checks (repeated values, gaps):
    # When a value is flagged save out a figure (PDF) of the flagged value +/- 24 hours on either side

# Evaluate function thresholds (where applicable, start with HadISD values and tune iteratively based on below results):
    # Review flagged values for each network, calculate false positive rate (things flagged that shouldn't be/total flagged vals)
    # Write out csv of values in station that were flagged by network QA/QC but not by ours
        # Summarize # values by that station's flags
        
        # PR for a QA/QC function ready for merging to main should include summary of above results
        # History of threshold tuning to get to target FP rate

## Step 1: whole station checks
    # Remove station if NA in lat-lon
    # Remove station if NA elevation or elevation out of range (must be >= -282 ft below sea level - Death Valley)
    # Remove station if outside of WECC

## Function: Conducts whole station qa/qc checks (lat-lon, within WECC, elevation)
def whole_station_qaqc(network, cleandir, qaqcdir):

    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
    end_api = datetime.datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    try:
        files = [] # Get files
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
            file = str(item.key)
            files += [file]

        # Get cleaned station file and read in metadata
        station_file = [file for file in files if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj['Body'].read()))
        stations = station_file['ERA-ID'].dropna()

        files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of cleaned file names

    except Exception as e:
        print(e) # testing
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append("Error in whole network: {}".format(e))

    else: # if files successfully read in
        # for station in stations: # full run
        for station in stations.sample(5): # TESTING SUBSET
            file_name = cleandir+station+".nc"

            if file_name not in files: # dont run qa/qc on a station that isn't cleaned
                print("{} was not cleaned - skipping qa/qc".format(station))
                errors['File'].append(station)
                errors['Time'].append(end_api)
                errors['Error'].append("No cleaned data for this station, does not proceed to qa/qc: see cleaned station list for reason")
                continue
            else:
                print('Running QA/QC on: ', station) # testing

                try:
                    fs = s3fs.S3FileSystem()
                    aws_url = "s3://wecc-historical-wx/"+file_name

                    with fs.open(aws_url) as fileObj:
                        ds = xr.open_dataset(fileObj, engine='h5netcdf') # CHECK THE ENGINE HERE

                        ## Add qc_flag variable for all variables, including elevation; defaulting to nan for fill value that will be replaced with qc flag
                        exclude_qaqc = ["time", "station", "lat", "lon"] # lat and lon have a different qc check
                        raw_qc_vars = [] # qc_variable for each data variable, will vary station to station

                        for var in ds.variables:
                             if '_qc' in var:
                                  raw_qc_vars.append(var) # raw qc variables, need to keep for comparison, then drop

                        for var in ds.variables:
                            if var not in exclude_qaqc and var not in raw_qc_vars:
                                qc_var = var + "_qc" # variable/column label
                                # qc_vars.append(qa_var)
                                ds[qc_var] = xr.full_like(ds[var], np.nan) # adds new variable in shape of original variable with designated nan fill value
                                                
                        stn_to_qaqc = ds.to_dataframe()

                        ## QAQC Part 1 (Whole Network) Functions
                        ## Lat-lon -- does not proceed through qaqc if failure
                        stn_to_qaqc = qaqc_missing_latlon(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{0} has a missing lat-lon, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_missing_latlon')
                            continue # skipping station
                        print('pass qaqc_missing_latlon') #testing

                        ## Within WECC -- does not proceed through qaqc if failure
                        stn_to_qaqc = qaqc_within_wecc(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{0} lat-lon is out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_within_wecc')
                            continue # skipping station
                        print('pass qaqc_within_wecc') #testing

                        ## Elevation -- if DEM in-filling fails, does not proceed through qaqc
                        stn_to_qaqc = qaqc_elev_range(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{} elevation out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Failure on qaqc_elev_range')
                        print('pass qaqc_elev_range') # testing

                        stn_to_qaqc = qaqc_elev_infill(stn_to_qaqc) # nan infilling must be before range check
                        if len(stn_to_qaqc.index) == 0:
                            print('DEM in-filling for {} failed, may not mean station does not pass qa/qc -- check'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('DEM in-filling error, may not mean station does not pass qa/qc -- check')
                            continue # skipping station

                        print(stn_to_qaqc) # testing
                
                except Exception as e:
                    # print(e) # testing
                    errors['File'].append(station)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Cannot read files in from AWS: {0}".format(e))










## -------------------------------------------------------------------------------------------------------------------------------------------

## Step 2: QA/QC (these could potentially be individual steps too)
## Testing with random subsample to ensure acceptably low level of false positive rates for specified thresholds and tests for each variables

# How to handle existing QA/QC flags -- use tracked qa/qc flags from DATA_clean.py

# Sensor height checks
## Wind must be 10m, temperature at 2m (flag if outside tolerance of +/- 0.3048 m)

# Internal variable consistency checks
    # logical checks for climate variables
        # Precipitation should not be negative
        # Precip accum in shorter time interval less than accum in longer time interval
        # Relative humidity values between 0 and 100 (or 0 and 1)
        # Flag values outside of world record for N AM
    # Gaps for jumps in a value
        # Use delta distributions to test
    # Climatological checks against distributions
        # WMO: 0.3 - 99.7 as percentiles
        # Generate distribution to compare to by using monthly temps of all time at that station
    # Repeated values
    # Values occurring more frequently than expected (histograms with >0.5 in a given bin)

# Cross-variable checks
    # Logic checks
        # Wind speed and direction
        # If speed is 0, direction should be 0
        
        # Dew-point not greater than air temp

# Drop original data variables to only keep qa/qc processed, converted variables

