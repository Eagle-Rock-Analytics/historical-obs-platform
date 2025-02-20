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
import scipy.stats as stats

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## QA/QC helper functions
# -----------------------------------------------------------------------------
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir


# -----------------------------------------------------------------------------
def open_log_file_merge(file):
    global log_file
    log_file = file


# -----------------------------------------------------------------------------
# Log print auxiliary functions
def printf(*args, verbose=True, log_file=None, **kwargs):
    import datetime

    tLog = lambda: datetime.datetime.utcnow().strftime("%m-%d-%Y %H:%M:%S") + " : \t"
    args = [str(a) for a in args]

    if verbose:
        if log_file is not None:
            print(" ".join([tLog(), *args]), **kwargs) or print(
                " ".join([tLog(), *args]), file=log_file, **kwargs
            )
        else:
            print(" ".join([tLog(), *args]), **kwargs)
    else:
        if log_file is not None:
            print(" ".join([tLog(), *args]), file=log_file, **kwargs)
        else:
            pass


# -----------------------------------------------------------------------------
def custom_sum(df):
    return df.apply(lambda x: np.nan if x.isna().all() else x.sum())


# -----------------------------------------------------------------------------
def hourly_standardization(df, verbose=verbose):
    """

    Resamples meteorological variables to hourly timestep according to standard conventions.

    Rules
    ------
        1.) top of the hour: take the first value in each hour
            - standard convention for temperature, dewpoint, wind speed, direction, relative humidity, air pressure
        2.) summation across hour: sum observations within each hour
            - standard convention for precipitation and solar radiation
        3.) constant across the hour: take the first value in each hour
            - this applies to variables, like station name and location, that do not change within each hour

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
    """

    printf(
        "Running: hourly_standardization",
        verbose=verbose,
        log_file=log_file,
        flush=True,
    )

    ##### define the variables for each sub-dataframe #####

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
        "sfcWind_dir_eraqc" "elevation_eraqc",
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
            result_list = []
            # Performing hourly aggregation
            if len(constant_df.columns) == 0:
                pass
            else:
                constant_result = constant_df.resample("1h", on="time").first()
                result_list.append(constant_result)

            if len(instant_df.columns) == 0:
                pass
            else:
                instant_result = instant_df.resample("1h", on="time").first()
                result_list.append(instant_result)

            if len(sum_df.columns) == 0:
                pass
            else:
                sum_result = sum_df.resample("1h", on="time").apply(
                    lambda x: np.nan if x.isna().all() else x.sum(skipna=True)
                )
                result_list.append(sum_result)

            if len(qaqc_df.columns) == 0:
                pass
            else:
                qaqc_result = qaqc_df.resample("1h", on="time").apply(
                    lambda x: ",".join(x.unique())
                )  # adding unique flags
                result_list.append(qaqc_result)

            # Aggregating and outputting reduced dataframe

            #result_list = [sum_result, instant_result, constant_result, qaqc_result]

            result = reduce(
                lambda left, right: pd.merge(left, right, on=["time"], how="outer"),
                result_list,
            )
            return result

    except Exception as e:
        printf(
            "hourly_standardization failed with Exception: {0}".format(e),
            verbose=verbose,
            log_file=log_file,
            flush=True,
        )
        return None
