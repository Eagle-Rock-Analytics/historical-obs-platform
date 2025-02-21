"""
This is a script where Stage 4: merge data related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
from functools import reduce
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

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"


## QA/QC helper functions
# -----------------------------------------------------------------------------
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir


## Standardization hlper functions
# -----------------------------------------------------------------------------
def custom_sum(df):
    return df.apply(lambda x: np.nan if x.isna().all() else x.sum())


# -----------------------------------------------------------------------------
def merge_hourly_standardization(df, verbose=verbose):
    """

    Resamples meteorological variables to hourly timestep according to standard conventions.

    Rules
    ------


    Parameters
    ------
    df: pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    verbose: boolean
        input for printf() to print to log file - set in script initialization

    Returns
    -------
    if success:
        df [pd.DataFrame]
            QAQC dataframe with all columns resampled to one hour (column name retained)
    if failure:
        None

    Notes
    -----
    1. Top of the hour: take the first value in each hour. Standard convention for temperature, dewpoint, wind speed, direction, relative humidity, air pressure
    2. Summation across hour: sum observations within each hour. Standard convention for precipitation and solar radiation
    3. Constant across the hour: take the first value in each hour. This applies to variables, like station name and location, that do not change within each hour
    """

    # Variables that remain constant within each hour
    constant_vars = [
        "time",
        "station",
        "lat",
        "lon",
        "elevation",
        "anemometer_height_m",
        "thermometer_height_m",
        "sfcWind_method",
        "pr_duration",
        "hour",
        "day",
        "month",
        "year",
        "date",
    ]

    # Aggregation across hour variables, standard meteorological convention: precipitation and solar radiation
    sum_vars = ["time", "pr", "pr_localmid", "pr_24h", "pr_5min", "pr_1h", "rsds"]

    # Top of the hour variables, standard meteorological convention: temperature, dewpoint temperature, pressure, humidity, winds
    instant_vars = [
        "time",
        "tas",
        "tdps",
        "tdps_derived",
        "ps",
        "psl",
        "ps_altimeter",
        "ps_derived",
        "hurs",
        "sfcwind",
        "sfcwind_dir",
    ]

    # QAQC flags, which remain constants within each hour
    qaqc_vars = [
        "tas_qc",
        "tas_eraqc",
        "pr_5min_eraqc",
        "pr_1h_eraqc",
        "pr_5min_qc",
        "pr_eraqc",
        "pr_depth_qc",
        "ps_qc",
        "ps_altimeter_qc",
        "ps_eraqc",
        "ps_altimeter_eraqc",
        "psl_qc",
        "psl_eraqc",
        "tdps_qc",
        "tdps_eraqc",
        "sfcWind_qc",
        "sfcWind_dir_qc",
        "sfcWind_eraqc",
        "sfcWind_dir_eraqc",
        "elevation_eraqc",
        "qaqc_process",
    ]

    # All variables, necessary for producing columns with hourly counts for each variable
    all_vars = constant_vars + sum_vars + instant_vars + qaqc_vars

    ##### Subset the dataframe according to rules
    constant_df = df[[col for col in constant_vars if col in df.columns]]

    qaqc_df = df[[col for col in qaqc_vars if col in df.columns if col != "time"]]
    qaqc_df = qaqc_df.astype(str)
    qaqc_df.insert(0, "time", df["time"])

    sum_df = df[[col for col in sum_vars if col in df.columns]]

    instant_df = df[[col for col in instant_vars if col in df.columns]]

    #####

    try:
        # if station does not report any variable, bypass
        if len(df.columns) == 0:
            printf(
                "Empty dataset - bypassing hourly aggregation",
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            return df
        else:
            # Performing hourly aggregation
            constant_result = constant_df.resample("1h", on="time").first()
            instant_result = instant_df.resample("1h", on="time").first()
            sum_result = sum_df.resample("1h", on="time").apply(
                lambda x: np.nan if x.isna().all() else x.sum(skipna=True)
            )
            qaqc_result = qaqc_df.resample("1h", on="time").apply(
                lambda x: ",".join(x.unique())
            )  # adding unique flags

            # Aggregating and outputting reduced dataframe
            result_list = [sum_result, instant_result, constant_result, qaqc_result]
            result = reduce(
                lambda left, right: pd.merge(left, right, on=["time"], how="outer"),
                result_list,
            )
            return result

    except Exception as e:
        # printf(
        #     "hourly_standardization failed with Exception: {0}".format(e),
        #     verbose=verbose,
        #     log_file=log_file,
        #     flush=True,
        # )
        # conver to logger version
        return None
