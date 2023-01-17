"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
from datetime import datetime
import numpy as np
import pandas as pd
import xarray as xr
import boto3
import s3fs


## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


#----------------------------------------------------------------------
## QA/QC Helper functions
# Given a network name, return all relevant AWS filepaths for other functions
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir

#----------------------------------------------------------------------
## Part 1 functions (whole station/network)
## missing spatial coords (lat-lon)
def qaqc_missing_latlon(df):
    """
    Checks if latitude and longitude is missing for a station.
    If missing, station is flagged to not proceed through QA/QC.
    Args:
        df: dataframe of cleaned station data to check
    Returns:

    """
    # first check if latitude is present in the dataframe
    if 'lat' in df.columns:
        if df['lat'].isnull().values.any() == True:
            print('missing latitude, skipping')
            # errors['File'].append(file)
            # errors['Time'].append(end_api)
            # errors['Error'].append('Missing latitude, skipping qa/qc.')
            # continue
    else:
        print('df has no latitude coordinate')

    # check if longitude is present in the dataframe
    if 'lon' in df.columns:
        if df['lon'].isnull().values.any() == True:
            print('missing longitude, skipping')
            # errors['File'].append(file)
            # errors['Time'].append(end_api)
            # errors['Error'].append('Missing longitude, skipping qa/qc.')
            # continue
    else:
        print('df has no longitude coordinate')

    return df








## in bounds of WECC
def qaqc_within_wecc():
    """
    Checks if station is within terrestrial & marine WECC boundaries.
    If outside of boundaries, station is flagged to not proceed through QA/QC.
    """

## missing/out of range elevation
# will involve DEM
def qaqc_elev_check(df):
    """
    Checks if elevation is outside of range of reasonable values for WECC region.
    If elevation is NA/missing, fill in elevation from DEM.
    """
    # Step 1: identify range of reasonable values
    # Step 2: DEM
    # Step 3: missing elevation check
    if df['elevation'].isnull().values.any() == True:
        print('missing elevation, will be processed via DEM')
    else:
        print('all valid elevation values')

#----------------------------------------------------------------------
## Part 2 functions
## Time conversions
## Need function to calculate sub-hourly to hourly


#----------------------------------------------------------------------
# To do
# establish false positive rate
# unit tests on each of these functions(?)
