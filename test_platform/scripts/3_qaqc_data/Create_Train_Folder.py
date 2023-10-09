"""
This script performs qa/qc protocols for cleaned station data for ingestion into the Historical Observations Platform, and is
independent of network. 
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
import datetime
import pandas as pd
import xarray as xr
import boto3
import s3fs
from io import BytesIO, StringIO
import argparse
import time

## Import qaqc stage calc functions
try:
    from calc_qaqc import *
except:
    print("Error importing calc_qaqc.py")

if __name__ == "__main__":

    # Create parser
    parser = argparse.ArgumentParser(
        prog="Create_Train_Folder",
        description="""This script creates Train_Files folder and downloads a pre-defined list of files to it""",
    )
    
    # Define arguments for the program
    parser.add_argument('-n', '--nFiles', default=46, help="Number of files to download (18 networks, 3 of each --if there are 3 available--)", type=int)
    parser.add_argument('-v', '--verbose', default=True, help="Verbose", type=bool)
    
    # Parse arguments
    args = parser.parse_args()
    nFiles = args.nFiles
    verbose = args.verbose

    # Read file list
    if verbose:
        print("Reading local file list")
    files = pd.read_csv("local_training_station.csv")

    # Create directory for local train files
    if verbose:
        print("Creating Folder")
    train_dir = "./Train_Files"
    if not os.path.exists(train_dir):
        os.makedirs(train_dir)

    # Setting up AWS
    if verbose:
        print("AWS setup")
        
    ## Set AWS credentials
    s3 = boto3.resource("s3")
    s3_cl = boto3.client('s3') # for lower-level processes
    bucket = s3.Bucket(bucket_name)
    
    ## Set relative paths to other folders and objects in repository.
    bucket_name = "wecc-historical-wx"
    
    if verbose:
        print("Downloading files")
    for station,network in zip(files['era.id'], files['network']):
        rawdir, cleandir, qaqcdir, mergedir = get_file_paths(network)
        
        file_name = cleandir+station+".nc"
        local_file_name = '{}/{}.nc'.format(train_dir, station)

        t0 = time.time()
        bucket.download_file(file_name, local_file_name)
        print("{} downloaded. Time ellapsed: {:.2f}".format(local_file_name, time.time()-t0))
