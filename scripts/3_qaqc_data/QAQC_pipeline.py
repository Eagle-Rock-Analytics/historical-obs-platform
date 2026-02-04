"""
QAQC_pipeline.py

This script performs qa/qc protocols for cleaned station data for ingestion into the Historical Observations Platform, and is
independent of network.
Approach:
(1) Remove duplicate stations
(2) Handle variables that report at different intervals and/or change frequency over time (convert to hourly?)
(3) QA/QC testing, including consistency checks, gaps, checks against climatological distributions, and cross variable checks.

Intended Use
-------------
This script runs the full QA/QC pipeline for a single weather station dataset. It processes cleaned output and applies a
    sequence of QA/QC tests to generate a QA/QC'd station file.

Inputs
------
Cleaned data for an individual network

Outputs
-------
QA/QC-processed data for an individual network, priority variables, all times. Organized by station as .zarr.
"""

import os
import sys
import datetime
import pandas as pd
import numpy as np
import xarray as xr
import boto3
import s3fs
from io import StringIO
import time
import tempfile
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from paths import BUCKET_NAME, RAW_WX, CLEAN_WX, QAQC_WX, MERGE_WX, STATIONS_CSV_PATH

try:
    from log_config import setup_logger
    from qaqc_plot import *
    from qaqc_utils import *
    from qaqc_wholestation import *
    from qaqc_logic_checks import *
    from qaqc_buoy_check import *
    from qaqc_frequent import *
    from qaqc_unusual_gaps import *
    from qaqc_unusual_large_jumps import *
    from qaqc_climatological_outlier import *
    from qaqc_unusual_streaks import *
    from qaqc_deaccumulate import *
except Exception as e:
    print(f"Error importing qaqc script: {e}")


s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Make local directories to save files temporarily and save timing
dirs = ["./qaqc_logs/"]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)


def setup_error_handling() -> tuple[dict[str, list], str, str]:
    """
    Sets-up error handling.

    Params
    ------
    None

    Returns
    -------
    errors : dict[str, list]
        dictionary of file, timing, and error message
    end_api : datetime
        time at beginnging of data download
    tiemstamp: datetime
        time at runtime
    """
    errors = {"File": [], "Time": [], "Error": []}  # Set up error handling
    end_api = datetime.datetime.now().strftime(
        "%Y%m%d%H%M"
    )  # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")
    return errors, end_api, timestamp


def print_qaqc_failed(
    errors: dict[str, list],
    station: str | None = None,
    end_api: str | None = None,
    message: str | None = None,
    test: str | None = None,
) -> dict[str, list]:
    """
    QAQC failure messaging

    Parameters
    ----------
    errors : dict
        dictionary of file, timing, and error message
    station : str, optional
        station name
    end_api : datetime, optional
        time at beginning of data download
    message : str, optional
        error message
    test : str, optional
        QAQC test name to include in error message

    Returns
    -------
    errors : dict[str, list]
        updated error dictionary with failure entry appended
    """
    logger.info(f"{station} {message}")
    errors["File"].append(station)
    errors["Time"].append(end_api)
    errors["Error"].append(f"Failure on {test}")
    return errors


def file_on_s3(df: pd.DataFrame, zarr: bool) -> pd.Series:
    """
    Check if network files are in s3 bucket

    Parameters
    ----------
    df : pd.DataFrame
        Table with information about each network and station
    zarr : bool
        Search the folder for zarr stores (zarr=True) or netcdfs (zarr=False)?

    Returns
    -------
    substring_in_filepath : bool
        True/False: Is the file in the s3 bucket?
    """

    files = []  # Get files
    for item in s3.Bucket(BUCKET_NAME).objects.filter(Prefix=df["cleandir"].iloc[0]):
        file = str(item.key)
        files += [file]

    # Depending on file type, modify the filepaths differently
    # The goal is to check if the era-id is contained in each substring
    if zarr == False:  # Get netcdf files
        file_st = [f.split(".nc")[0].split("/")[-1] for f in files if f.endswith(".nc")]
    elif zarr == True:  # Get zarrs
        # We just want to get the top directory for each station, i.e. "VALLEYWATER_6001.zarr/"
        # Since each station has a bunch of individual zarr stores, the split() function returns many copies of the same string
        # We just want one filename per station
        file_st = [
            file.split(".zarr/")[0].split("/")[-1] for file in files if ".zarr/" in file
        ]
        file_st = [x for i, x in enumerate(file_st) if x not in file_st[:i]]

    substring_in_filepath = df["era-id"].isin(file_st)
    return substring_in_filepath


def read_network_files(network: str, zarr: bool) -> pd.DataFrame:
    """
    Read files for a network from AWS

    Parameters
    ----------
    network : str
        Name of network
    zarr : bool
        Search the folder for zarr stores (zarr=True) or netcdfs (zarr=False)?

    Returns
    -------
    final_df : pd.DataFrame
        dataframe of cleaned data
    """

    # Read csv from s3
    full_df = pd.read_csv(STATIONS_CSV_PATH).loc[:, ["era-id", "network"]]

    # Add path info as new columns
    full_df["rawdir"] = full_df["network"].apply(lambda row: f"{RAW_WX}/{row}/")
    full_df["cleandir"] = full_df["network"].apply(lambda row: f"{CLEAN_WX}/{row}/")
    full_df["qaqcdir"] = full_df["network"].apply(lambda row: f"{QAQC_WX}/{row}/")
    full_df["mergedir"] = full_df["network"].apply(lambda row: f"{MERGE_WX}/{row}/")

    # If its a zarr store, use the zarr file extension (".zarr")
    if zarr == True:
        full_df["key"] = full_df.apply(
            lambda row: row["cleandir"] + row["era-id"] + ".zarr", axis=1
        )

    # If its a netcdf, use the netcdf file extension (".nc")
    elif zarr == False:
        full_df["key"] = full_df.apply(
            lambda row: row["cleandir"] + row["era-id"] + ".nc", axis=1
        )
    full_df["exist"] = np.zeros(len(full_df)).astype("bool")

    # To use the full dataset for specific sample stations
    df = full_df.copy()[full_df["network"] == network]

    # subset for specific network
    for n in df["network"].unique():
        ind = df["network"] == n
        df.loc[ind, "exist"] = file_on_s3(df[ind], zarr=zarr)
    df = df[df["exist"]]

    df["file_size"] = df["key"].apply(
        lambda row: s3_cl.head_object(Bucket=BUCKET_NAME, Key=row)["ContentLength"]
    )
    df = df.sort_values(by=["file_size", "network", "era-id"]).drop(columns="exist")

    # Evenly distribute df by size to help with memory errors
    # Number of groups is the total number of stations divided by the node size
    num_groups = len(df) // (72 * 3)
    total_size = df["file_size"].sum()
    target_size = total_size / num_groups

    groups = []
    current_group = []
    current_group_size = 0

    # Sort DataFrame by size to improve grouping efficiency
    df_sorted = df.sort_values(by="file_size", ascending=False)

    for index, row in df_sorted.iterrows():
        if current_group_size + row["file_size"] > target_size and current_group:
            groups.append(pd.DataFrame(current_group))
            current_group = []
            current_group_size = 0

        current_group.append(row)
        current_group_size += row["file_size"]

    if current_group:
        groups.append(pd.DataFrame(current_group))

    # Create a new DataFrame to hold the groups
    final_df = (
        pd.concat([df.assign(Group=i) for i, df in enumerate(groups)])
        .reset_index(drop=True)
        .drop(columns="Group")
    )

    return final_df


def process_output_ds(
    df: pd.DataFrame,
    attrs: dict[str, str],
    var_attrs: dict[str, str],
    network: str,
    timestamp: datetime.datetime,
    station: str,
    qaqcdir: str,
    zarr: bool,
):
    """
    Processes the final dataset for export to AWS.

    Parameters
    ----------
    df : pd.DataFrame
        data that has completed QAQC
    attrs : list of str
        attributes to be reset on xr.Dataset
    var_attrs: list of str
        variable attributes to be reset on xr.Dataset
    network: str
        network name
    timestamp: datetime.datetime
        time at runtime
    station: str
        station name
    qaqcdir : str
        path to QAQC AWS directory
    zarr : bool
        if True, output is a .zarr. if False, output is a .nc

    Returns
    -------
    None

    To Do
    -----
    1. Add metadata for _eraqc variables.
    """

    # Convert back to dataset
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        ds = df.to_xarray()

    # Inherit variable attributes
    ds = ds.assign_attrs(attrs)
    for var, value in var_attrs.items():
        ds[var] = ds[var].assign_attrs(value)

    # Update ancillary variables

    for eraqc_var in list(ds.data_vars.keys()):
        if "_eraqc" in eraqc_var:
            var = eraqc_var.split("_eraqc")[0]
            # sfcWind_eraqc is added to dataset by (`qaqc_sensor_height_w`) even if sfcWind is not.
            # We need to account this to avoid errors in ds[var] for sfcWind
            if var in list(ds.data_vars.keys()):
                # Only if var was originally present in dataset
                if "ancillary_variables" in list(ds[var].attrs.keys()):
                    ds[var].attrs["ancillary_variables"] = (
                        ds[var].attrs["ancillary_variables"] + f", {eraqc_var}"
                    )
                else:
                    ds[var].attrs["ancillary_variables"] = f"{eraqc_var}"
    # Overwrite file title
    ds = ds.assign_attrs(title=network + " quality controlled")

    # Append qaqc to the file history and comments (https://docs.unidata.ucar.edu/netcdf-c/current/attribute_conventions.html)
    ds.attrs["history"] = (
        ds.attrs["history"] + f" \nQAQC_pipeline.py script run on {timestamp} UTC"
    )
    ds.attrs["comment"] = (
        ds.attrs["comment"]
        + " \nAn intermediate data product: subject to cleaning but may not be subject to full QA/QC processing."
    )

    # Write station file
    if zarr == False:
        filename = station + ".nc"  # Make file name
    elif zarr == True:
        filename = station + ".zarr"
    filepath = qaqcdir + filename  # Writes file path

    # Push file to AWS with correct file name
    t0 = time.time()
    logger.info(
        f"Saving/pushing {filename} with dims {ds.dims} to {BUCKET_NAME}/{qaqcdir}"
    )
    if zarr == False:  # Upload as netcdf
        s3.Bucket(BUCKET_NAME).upload_file(filename, filepath)
    elif zarr == True:
        filepath_s3 = f"s3://{BUCKET_NAME}/{qaqcdir}{filename}"

        # Delete existing zarr store before writing
        # This avoids inconsistencies if there are differences between old vs. new
        # mode="w" only overwrites arrays with matching names, leaving old arrays behind
        fs = s3fs.S3FileSystem()
        if fs.exists(filepath_s3):
            logger.info(
                f"Found existing zarr at path {filepath_s3}. Deleting existing file and writing new zarr to the same path."
            )
            fs.rm(filepath_s3, recursive=True)

        # Write data to zarr
        ds.to_zarr(
            filepath_s3,
            consolidated=True,  # https://docs.xarray.dev/en/stable/internals/zarr-encoding-spec.html
            mode="w",  # Write & overwrite if file with same name exists already
        )

        ttime = time.time() - t0
        logger.info(f"Done saving/pushing file to AWS. Ellapsed time: {ttime:.2f} s.")
    ds.close()
    del ds

    return None


def qaqc_ds_to_df(
    ds: xr.Dataset,
) -> tuple[pd.DataFrame, pd.Index, dict[str, str], dict[str, str], list[str]]:
    """Converts xarray ds for a station to pandas df in the format needed for the pipeline

    Parameters
    ----------
    ds : xr.Dataset
        input data from the clean step

    Returns
    -------
    df : pd.DataFrame
        converted xr.Dataset into dataframe
    df_multi_idx : pd.Index
        multi-index of station and time
    attrs : dict[str, str]
        attributes from xr.Dataset
    var_attrs : dict[str, str]
        variable attributes from xr.Dataset
    era_qc_vars : list of str
        QAQC variables
    """

    # Add qc_flag variable for all variables, including elevation;
    # defaulting to nan for fill value that will be replaced with qc flag
    for key, val in ds.variables.items():
        if val.dtype == object:
            if key == "station":
                if str in [type(v) for v in ds[key].values]:
                    ds[key] = ds[key].astype(str)
            else:
                if str in [type(v) for v in ds.isel(station=0)[key].values]:
                    ds[key] = ds[key].astype(str)

    # lat, lon have different qc check
    exclude_qaqc = [
        "time",
        "station",
        "lat",
        "lon",
        "qaqc_process",
        "sfcWind_method",
        "pr_duration",
        "pr_depth",
        "PREC_flag",
        "rsds_duration",
        "rsds_flag",
        "hurs_temp",
        "hurs_temp_flag",
        "hurs_duration",
        "hurs_flag",
        "anemometer_height_m",
        "thermometer_height_m",
    ]

    raw_qc_vars = []  # qc_variable for each data variable, will vary station to station
    era_qc_vars = []  # our ERA qc variable

    for var in ds.data_vars:
        if "q_code" in var:
            # raw qc variable, need to keep for comparison, then drop
            raw_qc_vars.append(var)
            # raw qc variables, need to keep for comparison, then drop
        if "_qc" in var:
            raw_qc_vars.append(var)

    logger.info(f"Existing observation and QC variables: {list(ds.keys())}")

    # only in-fill nans for valid variables
    for var in ds.data_vars:
        if var not in exclude_qaqc and var not in raw_qc_vars and "_eraqc" not in var:
            qc_var = var + "_eraqc"  # variable/column label

            # if qaqc var does not exist, adds new variable in shape of original variable with designated nan fill value
            if qc_var not in era_qc_vars:
                ds = ds.assign({qc_var: xr.ones_like(ds[var]) * np.nan})
                era_qc_vars.append(qc_var)
                logger.info(f"nans created for {qc_var}")
                ds = ds.assign({qc_var: xr.ones_like(ds[var]) * np.nan})

    n_qc = len(era_qc_vars)  # determine length of eraqc variables per station
    logger.info(f"Created {n_qc} era_qc variables: {era_qc_vars}")

    # Save attributes to inheret them to the QAQC'ed file
    attrs = ds.attrs
    var_attrs = {var: ds[var].attrs for var in list(ds.data_vars.keys())}

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        df = ds.to_dataframe()

    # instrumentation heights
    if "anemometer_height_m" not in df.columns:
        try:
            df["anemometer_height_m"] = (
                np.ones(ds["time"].shape) * ds.anemometer_height_m
            )
        except:
            logger.info("Filling anemometer_height_m with NaN.")
            df["anemometer_height_m"] = np.ones(len(df)) * np.nan

    if "thermometer_height_m" not in df.columns:
        try:
            df["thermometer_height_m"] = (
                np.ones(ds["time"].shape) * ds.thermometer_height_m
            )
        except:
            logger.info("Filling thermometer_height_m with NaN.")
            df["thermometer_height_m"] = np.ones(len(df)) * np.nan

    # De-duplicate time axis
    df = df[~df.index.duplicated()].sort_index()

    # Save station/time multiindex
    df_multi_idx = df.index
    station = df.index.get_level_values(0)
    df["station"] = station

    # Station pd.Series to str
    station = station.unique().values[0]

    # Convert time/station index to columns and reset index
    df = df.droplevel(0).reset_index()

    # Add time variables needed by multiple functions
    df["hour"] = pd.to_datetime(df["time"]).dt.hour
    df["day"] = pd.to_datetime(df["time"]).dt.day
    df["month"] = pd.to_datetime(df["time"]).dt.month
    df["year"] = pd.to_datetime(df["time"]).dt.year
    df["date"] = pd.to_datetime(df["time"]).dt.date

    return df, df_multi_idx, attrs, var_attrs, era_qc_vars


def run_qaqc_pipeline(
    ds: xr.Dataset,
    network: str,
    errors: dict,
    station: str,
    end_api: datetime.datetime,
    rad_scheme: str,
) -> tuple[pd.DataFrame, list[str], list[str], list[str]]:
    """
    Runs all QAQC functions.

    Parameters
    ----------
    ds : xr.Dataset
        input data
    network : str
        network name
    errors : dict
        errors dictionary
    station : str
        station name
    end_api : datetime
        time at beginnging of data download
    rad_scheme : str
        radiation handling scheme for qaqc_frequent

    Returns
    -------
    stn_to_qaqc : pd.DataFrame
        dataframe of QAQC data
    attrs : list of str
        attributes from original xr.Dataset
    var_attrs : list of str
        variable attributes from original xr.Dataset
    era_qc_vars : list of str
        QAQC variables

    Notes
    -----
    1. Order of operations
    - Data converted to pd.DataFrame for processing
    - Part 1a: Whole station checks - if failure, entire station does not proceed through QA/QC
    - Part 1b: Whole station checks - if failure, entire station does proceed through QA/QC
    - Part 2: Logic checks
    - Part 3: Distribution & time series checks
    """

    # First identify if the station has valid data variables
    ds = qaqc_eligible_vars(ds)
    if ds is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="does not report any data variables. Station failed QAQC.",
            test="qaqc_eligible_vars",
        )
        return [None] * 4  # whole station failure, skip to next station

    else:
        ds = ds

    # Convert from xarray ds to pandas df in the format needed for qaqc pipeline
    df, df_multidx, attrs, var_attrs, era_qc_vars = qaqc_ds_to_df(ds)

    # Close ds file, netCDF,HDF5 unclosed files can sometimes cause issues during the mpi4py cleanup phase.
    ds.close()
    del ds

    # =========================================================
    ## START QA/QC ASSESSMENT
    # ---------------------------------------------------------
    ## Part 1a: Whole station checks - if failure, entire station does not proceed through QA/QC

    t0 = time.time()
    logger.info("QA/QC whole station tests")
    # ---------------------------------------------------------
    ## Missing values -- does not proceed through qaqc if failure
    stn_to_qaqc = df.copy()  # Need to define before qaqc_pipeline, in case
    new_df = qaqc_missing_vals(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="has an unchecked missing value or does not report any observation variables. Station failed QAQC.",
            test="qaqc_missing_vals",
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_missing_vals")

    # ---------------------------------------------------------
    ## Lat-lon -- does not proceed through qaqc if failure
    new_df = qaqc_missing_latlon(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="missing lat-lon. Station failed QAQC.",
            test="qaqc_missing_latlon",
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_missing_latlon")

    # ---------------------------------------------------------
    ## Within WECC -- does not proceed through qaqc if failure
    new_df = qaqc_within_wecc(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="lat-lon is out of range for WECC. Station failed QAQC.",
            test="qaqc_within_wecc",
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_within_wecc")

    # ---------------------------------------------------------
    ## Elevation -- if DEM in-filling fails, does not proceed through qaqc
    new_df = qaqc_elev_infill(stn_to_qaqc)  # nan infilling must be before range check
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="DEM in-filling failed. Station failed QAQC.",
            test="DEM in-filling, may not mean station does not pass qa/qc -- check",
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_elev_infill")

    # ---------------------------------------------------------
    ## Elevation -- range within WECC
    new_df = qaqc_elev_range(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="elevation out of range for WECC. Station failed QAQC.",
            test="qaqc_elev_range",
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_elev_range")

    # ---------------------------------------------------------
    ## Elevation -- range consistency check
    new_df = qaqc_elev_internal_range_consistency(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="internal elevation range is inconsistent. Station failed QAQC.",
            test="qaqc_elev_internal_range_consistency",
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_elev_internal_range_consistency")

    # =========================================================
    ## Part 1b: Whole station checks - if failure, entire station does proceed through QA/QC
    # ---------------------------------------------------------
    ## Pressure units fix (temporary)
    new_df = qaqc_pressure_units_fix(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with world record check.",
            test="qaqc_pressure_units_fix",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_pressure_units_fix",
        )

    # ---------------------------------------------------------
    ## Precipitation de-accumulation
    new_df = qaqc_deaccumulate_precip(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with precip deaccumulation",
            test="qaqc_deaccumulate_precip",
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_deaccumulate_precip")

    ttime = time.time() - t0
    logger.info(f"Done whole station tests, Ellapsed time: {ttime:.2f} s.\n")

    # ---------------------------------------------------------
    ## World record checks: air temperature, dewpoint, wind, pressure
    new_df = qaqc_world_record(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with world record check",
            test="qaqc_world_record",
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_world_record")

    ttime = time.time() - t0
    logger.info("Done whole station tests, Ellapsed time: {ttime:.2f} s.\n")

    # =========================================================
    ## Part 2: Variable logic checks
    t0 = time.time()
    logger.info("QA/QC logic checks")
    # ---------------------------------------------------------
    ## dew point temp cannot exceed air temperature
    new_df = qaqc_crossvar_logic_tdps_to_tas_supersat(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with temperature cross-variable logic check",
            test="qaqc_crossvar_logic_tdps_to_tas_supersat",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_crossvar_logic_tdps_to_tas_supersat",
        )

    # ---------------------------------------------------------
    ## dew point temp cannot exceed air temperature (wet bulb drying)
    new_df = qaqc_crossvar_logic_tdps_to_tas_wetbulb(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with temperature cross-variable logic check",
            test="qaqc_crossvar_logic_tdps_to_tas_wetbulb",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_crossvar_logic_tdps_to_tas_wetbulb",
        )

    # ---------------------------------------------------------
    ## precipitation is not negative
    new_df = qaqc_precip_logic_nonegvals(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with negative precipitation values",
            test="qaqc_precip_logic_nonegvals",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_precip_logic_nonegvals",
        )

    # ---------------------------------------------------------
    ## precipitation duration logic
    new_df = qaqc_precip_logic_accum_amounts(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with precip duration logic check",
            test="qaqc_precip_logic_accum_amounts",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_precip_logic_accum_amounts",
        )
    # ---------------------------------------------------------
    ## wind direction should be 0 when wind speed is also 0
    new_df = qaqc_crossvar_logic_calm_wind_dir(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with wind cross-variable logic check",
            test="qaqc_crossvar_logic_calm_wind_dir",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_crossvar_logic_calm_wind_dir",
        )

    ttime = time.time() - t0
    logger.info("Done logic checks, Ellapsed time: {ttime:.2f} s.\n")

    # =========================================================
    ## Part 3: Distribution and timeseries checks - order matters!
    # ---------------------------------------------------------
    ## Buoys with known issues with specific qaqc flags
    ## NDBC and MARITIME only
    if network == "MARITIME" or network == "NDBC":
        t0 = time.time()
        logger.info("QA/QC bouy check")

        new_df = spurious_buoy_check(stn_to_qaqc, era_qc_vars)
        if new_df is None:
            errors = print_qaqc_failed(
                errors,
                station,
                end_api,
                message="Flagging problematic buoy issue",
                test="spurious_buoy_check",
            )
        else:
            stn_to_qaqc = new_df
            logger.info(
                "pass spurious_buoy_check",
            )

    ttime = time.time() - t0
    logger.info(f"Done QA/QC bouy check, Ellapsed time: {ttime:.2f} s.\n")

    # ---------------------------------------------------------
    # frequent values
    t0 = time.time()
    logger.info("QA/QC frequent values")

    new_df = qaqc_frequent_vals(stn_to_qaqc, rad_scheme=rad_scheme)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with frequent values function",
            test="qaqc_frequent_vals",
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_frequent_vals")

    ttime = time.time() - t0
    logger.info(f"Done QA/QC frequent values, Ellapsed time: {ttime:.2f} s.\n")

    # ---------------------------------------------------------
    # distribution / unusual gaps
    t0 = time.time()
    logger.info("QA/QC unusual gaps")

    new_df = qaqc_unusual_gaps(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with unusual gap distribution function",
            test="qaqc_unusual_gaps",
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_unusual_gaps")

    ttime = time.time() - t0
    logger.info(f"Done QA/QC unusual gaps, Ellapsed time: {ttime:.2f} s.\n")

    # ---------------------------------------------------------
    # climatological outliers
    t0 = time.time()
    logger.info("QA/QC climatological outliers")

    new_df = qaqc_climatological_outlier(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with climatological outlier check",
            test="qaqc_climatological_outlier",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_climatological_outlier",
        )

    ttime = time.time() - t0
    logger.info(f"Done QA/QC climatological outliers, Ellapsed time: {ttime:.2f} s.\n")

    # ---------------------------------------------------------
    # unusual streaks (repeated values)
    t0 = time.time()
    logger.info("QA/QC unsual repeated streaks")

    new_df = qaqc_unusual_repeated_streaks(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with unusual streaks (repeated values) check",
            test="qaqc_unusual_repeated_streaks",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_unusual_repeated_streaks",
        )

    ttime = time.time() - t0
    logger.info(f"Done QA/QC unsual repeated streaks, Ellapsed time: {ttime:.2f} s.\n")

    # ---------------------------------------------------------
    # unusual large jumps (spikes)
    t0 = time.time()
    logger.info("QA/QC unsual large jumps")

    new_df = qaqc_unusual_large_jumps(stn_to_qaqc)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with unusual large jumps (spike check) check",
            test="qaqc_unusual_large_jumps",
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_unusual_large_jumps",
        )

    ttime = time.time() - t0
    logger.info(f"Done QA/QC unsual large jumps, Ellapsed time: {ttime:.2f} s.\n")

    ## END QA/QC ASSESSMENT
    # =========================================================
    # Re-index to original time/station values
    # Calculate flag coverage per variable
    logger.info(
        "Summary of QA/QC flags set per variable",
    )
    flag_summary(stn_to_qaqc)

    stn_to_qaqc = stn_to_qaqc.set_index(df_multidx).drop(
        columns=["time", "hour", "day", "month", "year", "date", "station"]
    )

    # Sort by time and remove any overlapping timesteps
    # Check back to see if this can or needs to be removed
    stn_to_qaqc = stn_to_qaqc[~stn_to_qaqc.index.duplicated()].sort_index()

    return stn_to_qaqc, attrs, var_attrs, era_qc_vars


def run_qaqc_one_station(
    station: str, verbose: bool = False, rad_scheme: str = "remove_zeros"
):
    """
    Runs the full QA/QC pipeline on a single weather station dataset.

    This function handles file path setup, reads cleaned station data
    from AWS S3 (NetCDF or Zarr), runs a series of QA/QC tests, writes
    the processed output, and logs errors and status messages. It is
    intended to be used as part of the Historical Observations Platform.

    Parameters
    ----------
    station : str
        Unique identifier for the weather station (e.g., "LOXWFO_CBGC1").
    verbose : bool, optional
        If True, enables verbose logging. Default is False.
    rad_scheme : str, optional
        Strategy for handling solar radiation data. Default is "remove_zeros".

    Returns
    -------
    None
    """

    ## ======== BASIC SETUP ========

    # Read in csv file containing information about each station
    stations_df = pd.read_csv(STATIONS_CSV_PATH)
    station_row = stations_df[stations_df["era-id"] == station]

    # Check that the input station exists in the station list :)
    if len(station_row) == 0:
        print(f"No file found in records for station {station}.")
        return None

    # Get the network and directories for network data in AWS
    network = station_row["network"].item()
    raw_data_dir, cleaned_data_dir, qaqc_dir, merge_dir = get_file_paths(network)

    # Set zarr argument
    zarrified_networks = ["VALLEYWATER", "CW3E"]  # Networks with zarrified data
    if network in zarrified_networks:
        zarr = True
    else:
        zarr = False

    # Set up error handling
    t0 = time.time()
    errors, end_api, t0_str = setup_error_handling()

    # Initialize the logger with the specific log file name
    tsplit = t0_str.split(", ")[0]
    log_fname = f"qaqc_logs/qaqc_{station}.{tsplit}.log"
    logger = setup_logger(log_file=log_fname, verbose=verbose)
    logger.info(
        f"Starting QAQC for station: {station}\n",
    )

    ## ======== READ FILE FROM S3 ========

    # Create an s3fs filesystem object
    fs = s3fs.S3FileSystem()

    # Define path to file using correct file extension
    aws_url_no_extension = f"s3://{BUCKET_NAME}/{cleaned_data_dir}{station}"
    aws_url = aws_url_no_extension + ".zarr" if zarr else aws_url_no_extension + ".nc"

    # Open the file
    logger.info(f"Reading file from AWS S3...")
    try:
        if zarr:  # Open zarr
            ds = xr.open_zarr(aws_url)
        else:  # Open netcdf
            with fs.open(aws_url) as fileObj:
                # Now we use the open file handle with xarray, without closing it prematurely
                ds = xr.open_dataset(fileObj).load()
    except:
        print(
            f"{station} did not pass QA/QC because the file could not be opened and/or found in AWS. File path: {aws_url}"
        )
        exit()  # End script here

    logger.info(
        f"Done reading. Ellapsed time: {time.time() - t0} s.\n",
    )

    ## ======== MANUAL PREPROCESSING =========

    # Drop time duplicates, if they were not properly dropped in the cleaning process
    if ("station" in list(ds.dims.keys())) and ("time" in list(ds.dims.keys())):
        ds = ds.drop_duplicates(dim="time")
    # There are stations without time/station dimensions
    elif ("time" in list(ds.data_vars.keys())) and (
        "station" in list(ds.data_vars.keys())
    ):
        tt = ds["time"]
        ss = [ds["station"].values[0]]
        ds = (
            ds.drop(["time", "station"])
            .rename_dims(index="time")
            .expand_dims({"station": ss})
            .rename(index="time")
            .assign_coords(time=tt.values)
        )
    else:
        ds = ds.drop_duplicates(dim="time")

    ## ======== RUN FULL QAQC PIPELINE =========
    logger.info(f"Running QA/QC on: {station}\n")
    try:
        # Run main QAQC functions
        test = "run_qaqc_pipeline"  # Used for making error message
        df, attrs, var_attrs, era_qc_vars = run_qaqc_pipeline(
            ds,
            network,
            errors,
            station,
            end_api,
            rad_scheme,
        )
        if df is None:
            # No data is returned by qaqc_pipeline
            # Error handling should have happened within run_qaqc_pipeline
            # Thus, just skip right to the finally section
            return

        # Save file
        # Attributes are assigned to dataset
        # File is uploaded as a zarr/nc locally or to AWS
        test = "process_output_ds"
        process_output_ds(
            df,
            attrs,
            var_attrs,
            network,
            t0_str,
            station,
            qaqc_dir,
            zarr=True,  # Default to always write to zarr
        )

    except Exception as e:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message=f"{test} failed with error: {e}",
            test=test,
        )

    ## ======== FINISH =========
    finally:

        # Print elapsed time
        ttime = time.time() - t0
        logger.info(f"Script complete. Elapsed time: {ttime:.2f} s.\n")

        # Convert errors to DataFrame and
        errors_df = pd.DataFrame(errors)
        errors_s3_filepath = (
            f"s3://{BUCKET_NAME}/{qaqc_dir}qaqc_errs/errors_{station}_{end_api}.csv"
        )
        errors_df.to_csv(errors_s3_filepath)

        # Print error file location
        logger.info(f"errors saved to {errors_s3_filepath}\n")

        # Save log file to s3 bucket
        logfile_s3_filepath = f"s3://{BUCKET_NAME}/{qaqc_dir}{log_fname}"
        logger.info(
            f"Saving log file to {logfile_s3_filepath}\n",
        )
        s3.Bucket(BUCKET_NAME).upload_file(log_fname, f"{qaqc_dir}{log_fname}")

        # Close logging handlers manually
        for handler in logger.handlers:
            handler.close()
            logger.removeHandler(handler)

    return None
