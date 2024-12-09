"""
This is a script where Stage 4: merge data related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import numpy as np
import pandas as pd
import requests
import urllib
import datetime
import math
import shapely
import xarray as xr
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import scipy.stats as stats

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## QA/QC helper functions
#-----------------------------------------------------------------------------
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir

#-----------------------------------------------------------------------------
def precip_hourly_standardization(df):
    
    """
    Ensures that precipitation accumulation amounts are consistent with reporting time frame.
    Only needs to be applied when 2 or more precipitation duration specific variables are present (pr_5min, pr_1h, pr_24h)
    For example: pr_5min should not be larger than pr_1h

    Rules:
    ------
        1) pr_5min < pr_1h < pr_24h
        2) pr_localmid should never exceed pr_24h

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if success:
            df [pd.DataFrame]: QAQC dataframe with additional column "pr_1hr"
        if failure:
            None
    """

    print("Running: precip_hourly_standardization", flush=True)

    # identify which precipitation vars are reported by a station
    vars_to_remove = ['qc', 'duration', 'method', 'depth']
    all_pr_vars = [var for var in df.columns if 'pr' in var] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if not any(True for item in vars_to_remove if item in var)] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    
    # TODO: add to log file
    try:        
        # if station does not report any precipitation values, or only one, bypass
        if len(pr_vars) == 0 # or len(pr_vars) == 1:
            print('Station does not report precipitation variables - bypassing hourly aggregation', flush=True)
            return df
        if 'pr_1h' in pr_vars:
            print('Station already reports hourly precipitation - bypassing hourly aggregation', flush=True)
            return df
        if 'pr_24h' in pr_vars:
            print('Station reports daily precipitation - bypassing hourly aggregation', flush=True)
            return df
        if 'pr_5min' in pr_vars and 'pr_24h' not in pr_vars and 'pr_1h' not in pr_vars:
            
        return df

    except Exception as e:
        print("precip_hourly_standardization failed with Exception: {0}".format(e), flush=true)
        return None