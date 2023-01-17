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
def get_cleaned_stations(cleandir, qaqcdir):

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
        # files = [file for file in files if "error" not in file] # Remove error handling files.
        # files = [file for file in files if "station" not in file] # Remove error handling files.
        print('Number of cleaned files for testing: ', len(files)) # testing

    except Exception as e:
        print(e) # testing
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append(e)

    else: # if files successfully read in
        # for file in files: # full run
        files = sample(files,1)
        for file in files: # subset for testing
            print('Running QA/QC on: ', file) # testing

            try:
                fs = s3fs.S3FileSystem()
                aws_urls = ["s3://wecc-historical-wx/"+file for file in files]

                with fs.open(aws_urls[0]) as fileObj:
                    ds = xr.open_dataset(fileObj, engine='h5netcdf')
                    file_to_qaqc = ds.to_dataframe()
                    print(file_to_qaqc.head()) # testing

                    # qaqc_missing_latlon
                    qaqc_latlon = qaqc_missing_latlon(file_to_qaqc) # eventually in calc_qaqc

                    # if file_to_qaqc['lat'].isnull().values.any() == True:
                    #     print('{} has a missing latitude, skipping'.format(file))
                    #     errors['File'].append(file)
                    #     errors['Time'].append(end_api)
                    #     errors['Error'].append('Missing latitude, skipping qa/qc.')
                    #     continue
                    #
                    # if file_to_qaqc['lon'].isnull().values.any() == True:
                    #     print('{} has a missing longitude, skipping'.format(file))
                    #     errors['File'].append(file)
                    #     errors['Time'].append(end_api)
                    #     errors['Error'].append('Missing longitude, skipping qa/qc.')
                    #     continue

                    # qaqc_within_wecc -- IN PROGRESS
                    # qaqc_wecc = qaqc_within_wecc(file_to_qaqc) # eventually in calc_qaqc

                    # qaqc_elev_check -- IN PROGRESS
                    qaqc_elev = qaqc_elev_check(file_to_qaqc)

                    ds.close()
            except Exception as e:
                print(e) # testing
                errors['File'].append("Whole network")
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
    get_cleaned_stations(cleandir, qaqcdir)
