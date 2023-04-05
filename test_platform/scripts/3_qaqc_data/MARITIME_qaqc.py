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

# Step 0: Environment set-up
# Import libraries
import os
from datetime import datetime
# import numpy as np
import pandas as pd
import xarray as xr
import boto3
# from random import sample
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
# wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
# wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## Function: Conducts whole station qa/qc checks (lat-lon, within WECC, elevation)
def whole_station_qaqc(network, cleandir, qaqcdir):

    # Set up error handling.
    errors = {'File':[], 'Time':[], 'Error':[]} # Set up error handling
    end_api = datetime.now().strftime('%Y%m%d%H%M') # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    try:
        files = [] # Get files
        for item in s3.Bucket(bucket_name).objects.filter(Prefix = cleandir):
            file = str(item.key)
            files += [file]

        # Get cleaned station file and read in metadata
        station_file = [file for file in files if 'stationlist_' in file]
        obj = s3_cl.get_object(Bucket=bucket_name, Key=station_file[0])
        station_file = pd.read_csv(BytesIO(obj['Body'].read()))
        stations = station_file['STATION_ID'].dropna()

        files = list(filter(lambda f: f.endswith(".nc"), files)) # Get list of cleaned file names
        print('Number of cleaned files for testing: ', len(files)) # testing

    except Exception as e:
        print(e) # testing
        errors['File'].append("Whole network")
        errors['Time'].append(end_api)
        errors['Error'].append("Error in whole network: {}".format(e))

    else: # if files successfully read in
        # for station in stations: # full run
        for station in stations.sample(5): # TESTING SUBSET
            station = network + "_" + station
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
                        ds = xr.open_dataset(fileObj, engine='h5netcdf') # FLAG CHECK THE ENGINE HERE
                        stn_to_qaqc = ds.to_dataframe()
                        print(stn_to_qaqc.head()) # testing

                        ## Lat-lon
                        stn_to_qaqc = qaqc_missing_latlon(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{} has a missing lat-lon, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Missing lat or lon, skipping qa/qc.')
                            continue # skipping station

                        ## Within WECC
                        stn_to_qaqc = qaqc_within_wecc(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{} lat-lon is out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Latitude or Longitude out of range for WECC, skipping qa/qc')
                            continue # skipping station

                        ## Elevation
                        stn_to_qaqc = qaqc_elev_demfill(stn_to_qaqc) # nan infilling must be before range check
                        if len(stn_to_qaqc.index) == 0:
                            print('This station reports a NaN for elevation, infilling from DEM')
                            continue

                        stn_to_qaqc = qaqc_elev_check(stn_to_qaqc)
                        if len(stn_to_qaqc.index) == 0:
                            print('{} elevation out of range for WECC, skipping'.format(station)) # testing
                            errors['File'].append(station)
                            errors['Time'].append(end_api)
                            errors['Error'].append('Elevation out of range for WECC, skippinig qa/qc')
                            continue # skipping station

                except Exception as e:
                    print(e) # testing
                    errors['File'].append(station)
                    errors['Time'].append(end_api)
                    errors['Error'].append("Cannot read files in from AWS: {}".format(e))
                        
                print("{} passes qa/qc round 1".format(station)) # testing
                # Update global attributes
                ds = ds.assign_attrs(title = network+" quality controlled")
                ds = ds.assign_attrs(history = 'MARITIME_qaqc.py script run on {} UTC'.format(timestamp))
                ds = ds.assign_attrs(comment = 'Intermediate data product: may not have been subject to any cleaning or QA/QC processing.') # do we modify this now?


                # Write station file to netcdf format
                try:
                    filename = station + ".nc" # Make file name
                    filepath = qaqcdir + filename # Writes file path

                    # Write locally
                    ds.to_netcdf(path = 'temp/temp.nc', engine = 'netcdf4') # Save station file.

                    # Push file to AWS with correct file name
                    s3.Bucket(bucket_name).upload_file('temp/temp.nc', filepath)

                    print('Saving {} with dims {}'.format(filename, ds.dims))
                    ds.close() # Close dataframe

                except Exception as e:
                    print(filename, e)
                    errors['File'].append(filename)
                    errors['Time'].append(end_api)
                    errors['Error'].append('Error saving ds as .nc file to AWS bucket: {}'.format(e))
                    continue

    # Write errors to csv
    finally:
        errors = pd.DataFrame(errors)
        csv_buffer = StringIO()
        errors.to_csv(csv_buffer)
        content = csv_buffer.getvalue()
        s3_cl.put_object(Bucket=bucket_name, Body=content, Key=qaqcdir+"errors_{}_{}.csv".format(network, end_api)) # Make sure error files save to correct directory


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run function
if __name__ == "__main__":
    network = "NDBC"
    rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
    print(cleandir, qaqcdir) # testing
    whole_station_qaqc(network, cleandir, qaqcdir)

# To do:
# flag as attribute?  only files that pass get saved?
# check the h5netcdf vs. netcdf4 engine
# delete testing notes
