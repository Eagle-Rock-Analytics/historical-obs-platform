"""
This script performs qa/qc protocols for cleaned NDBC/MARITIME data for ingestion into the Historical Observations Platform.
Approach:
(1) Remove duplicate stations
(2) Handle variables that report at different intervals and/or change frequency over time (convert to hourly?)
(3) QA/QC testing, including consistency checks, gaps, checks against climatological distributions, and cross variable checks.
(4) Case study analysis for accuracy -- SHOULD THIS BE A SEPARATE SCRIPT/PROCESS?
Inputs: Cleaned data for an individual network
Outputs: QA/QC-processed data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

## TO DO LIST
## Any notes critical for further development, e.g.: AWS implementation

# Step 0: Environment set-up
# Import libraries
from datetime import datetime
import numpy as np
import pandas as pd
import xarray as xr
import boto3
from random import sample
import s3fs
import geopandas as gp

# To be able to open xarray files from S3, h5netcdf must also be installed, but doesn't need to be imported.

## Import qaqc stage calc functions
try:
    from calc_qaqc import *
except:
    print("Error importing calc_qaqc.py")


## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## Read in data, and subset random sample for testing
def get_cleaned_stations(cleandir, qaqcdir, terrpath, marpath):
    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    try:
        files = [] # Get files
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
            file = str(item.key)
            files += [file]

        files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of file names
        # files = [file for file in files if "error" not in file] # Remove error handling files
        print('Number of cleaned files for testing: ', len(files)) # testing

    except Exception as e:
        print(e) # testing
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    else: # if files successfully read in
        files = sample(files,5) # subset for testing
        for file in files:
            print('Running QA/QC on: ', file) # testing

            try:
                fs = s3fs.S3FileSystem()
                aws_url = "s3://wecc-historical-wx/"+file

                with fs.open(aws_url) as fileObj:
                    ds = xr.open_dataset(fileObj, engine='h5netcdf')
                    file_to_qaqc = ds.to_dataframe()
                    print(file_to_qaqc.head()) # testing

                    file_to_qaqc = qaqc_missing_latlon(file_to_qaqc)
                    if file_to_qaqc.empty == True:
                        print('{} has a missing lat-lon, skipping'.format(file)) # testing
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append('Missing lat or lon, skipping qa/qc.')
                        continue # skipping station
                    # else:
                    #     print("pass: lat lon present") # testing, will delete

                    file_to_qaqc = qaqc_within_wecc(file_to_qaqc)
                    if file_to_qaqc.empty == True:
                        print('{} lat-lon is out of range for WECC, skipping'.format(file)) # testing
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append('Latitude or Longitude out of range for WECC, skipping qa/qc')
                        continue # skipping station
                    # else:
                    #     print("pass: within wecc") # testiing, will delete

                    file_to_qaqc = qaqc_elev_check(file_to_qaqc)
                    if file_to_qaqc.empty == True:
                        print('{} elevation out of range for WECC, skipping'.format(file)) # testing
                        errors['File'].append(file)
                        errors['Time'].append(end_api)
                        errors['Error'].append('Elevation out of range for WECC, skippinig qa/qc')
                        continue # skipping station
                    # else:
                    #     print("pass: elevation in range, but no infilling") # testing, will delete


                    print("{} passes qa/qc round 1".format(file)) # testing
                    ds.close()


            except Exception as e:
                print(e) # testing
                errors['File'].append(file)
                errors['Time'].append(end_api)
                errors['Error'].append(e)



    # # Write errors to csv
    # finally:
    #     print(errors) # Testing
    #     errors = pd.DataFrame(errors)
    #     csv_buffer = StringIO()
    #     errors.to_csv(csv_buffer)
    #     content = csv_buffer.getvalue()
    #     s3_cl.put_object(Bucket=bucket_name, Body=content, Key=qaqcdir+"errors_{}_{}.csv".format(network, end_api)) # Make sure error files save to correct directory



# Run function
if __name__ == "__main__":
    rawdir, cleandir, qaqcdir, mergedir = get_file_paths("NDBC")
    print(cleandir, qaqcdir) # testing
    get_cleaned_stations(cleandir, qaqcdir, wecc_terr, wecc_mar)
