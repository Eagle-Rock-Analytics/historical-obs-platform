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
def hourly_standardization(df):
    """
    
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
            df [pd.DataFrame]: QAQC dataframe with all columns resampled to one hour
        if failure:
            None
    """

    print("Running: hourly_standardization", flush=True)

    ##### define the variables for each sub-dataframe #####
    constant_vars = ["time","station", "lat", "lon", "elevation",
                    "anemometer_height_m","thermometer_height_m",
                    'sfcWind_method',
                    'pr_duration',
                    "hour","day","month","year","date"]

    qaqc_vars = ["time","tas_qc", "tas_eraqc",
                "pr_5min_eraqc","pr_1h_eraqc","pr_5min_qc",'pr_eraqc','pr_depth_qc', 
                "ps_qc",'ps_altimeter_qc','ps_eraqc','ps_altimeter_eraqc', 
                'psl_qc','psl_eraqc',
                'tdps_qc','tdps_eraqc',
                'sfcWind_qc','sfcWind_dir_qc','sfcWind_eraqc','sfcWind_dir_eraqc'
                "elevation_eraqc", "qaqc_process"]

    #precipitatino and solar radiation
    sum_vars = ["time","tas",
                "pr","pr_localmid","pr_24h","pr_5min","pr_1h",
                "rsds"]

    #temperature, dewpoint, wind speed, wind direction, relative humidity, air pressure
    instant_vars = ["time",
                "tdps","tdps_derived",
                "ps","psl","ps_altimeter",
                "hurs",
                "sfcwind","sfcwind_dir",
                "total"]
    
    ##### subset the dataframe
    print("generating subsets", flush=True)

    constant_df = df[[col for col in constant_vars if col in df.columns]]

    qaqc_df = df[[col for col in qaqc_vars if col in df.columns]]
    qaqc_vars_subset = qaqc_df.columns.tolist()
    qaqc_vars_subset.remove("time")
    qaqc_df[qaqc_vars_subset] = qaqc_df[qaqc_vars_subset].astype(str)

    sum_df = df[[col for col in sum_vars if col in df.columns]]

    instant_df = df[[col for col in instant_vars if col in df.columns]]
    
    ##### 
    print("checking if dataset contains sub-hourly precipitation data", flush=True)
    # TODO: add to log file
    try:        
        # if station does not report any precipitation values, or only one, bypass
        if len(df.columns) == 0:
            print('Empty dataset - bypassing hourly aggregation', flush=True)
            return df
        else: 
            print('performing hourly aggregation', flush=True)
            constant_result =  constant_df.resample('1h',on='time').first()
            qaqc_result = qaqc_df.resample('1h',on='time').apply(lambda x: ','.join(x))
            instant_result =  instant_df.resample('1h',on='time').first()
            sum_result =  sum_df.resample('1h',on='time').sum().rename(columns={'pr_5min': 'pr_1hr_test'})

            print('producing final result', flush=True)
            result_list = [sum_result,instant_result,constant_result,qaqc_result]
            result = reduce(lambda  left,right: pd.merge(left,right,on=['time'],
                                            how='outer'), result_list)
            return result
        
    
    except Exception as e:
        print("precip_hourly_standardization failed with Exception: {0}".format(e), flush=True)
        return None