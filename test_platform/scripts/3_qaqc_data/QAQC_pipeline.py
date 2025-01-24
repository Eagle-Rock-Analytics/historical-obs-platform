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
from io import StringIO
import time
import tempfile
from mpi4py import MPI
import logging
from simplempi import simpleMPI

# from simplempi.parfor import parfor, pprint

# Import all qaqc script functions
try:
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
except Exception as e:
    print("Error importing qaqc script: {}".format(e))

# Set up directory to save files temporarily and save timing, if it doesn't already exist.
dirs = ["./temp/", "./timing/", "./local_qaqced_files/", "./qaqc_logs/"]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

from log_config import setup_logger

os.environ["HDF5_USE_FILE_LOCKING"] = "TRUE"
# ----------------------------------------------------------------------------
## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"

# ============================================================================
# Define global functions and variables


# ----------------------------------------------------------------------------
def setup_error_handling():
    """ """
    errors = {"File": [], "Time": [], "Error": []}  # Set up error handling
    end_api = datetime.datetime.now().strftime(
        "%Y%m%d%H%M"
    )  # Set end time to be current time at beginning of download: for error handling csv
    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")
    return errors, end_api, timestamp


# ----------------------------------------------------------------------------
def print_qaqc_failed(
    errors, station=None, end_api=None, message=None, test=None, verbose=False
):
    """ """
    logger.info(
        "{0} {1}, skipping station".format(station, message),
    )
    errors["File"].append(station)
    errors["Time"].append(end_api)
    errors["Error"].append("Failure on {}".format(test))
    return errors


# ----------------------------------------------------------------------------
## Check if network files are in s3 bucket
def file_on_s3(df, zarr):
    """Check if network files are in s3 bucket

    Parameters
    ----------
    df: pd.DataFrame
        Table with information about each network and station
    zarr: boolean
        Search the folder for zarr stores (zarr=True) or netcdfs (zarr=False)?

    Returns
    -------
    substring_in_filepath: boolean
        True/False: Is the file in the s3 bucket?

    """

    files = []  # Get files
    for item in s3.Bucket(bucket_name).objects.filter(Prefix=df["cleandir"].iloc[0]):
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


# ----------------------------------------------------------------------------
## Read network nc files
def read_network_files(network, zarr):
    """Read files for a network from AWS

    Parameters
    ----------
    network: str
        Name of network
    zarr: boolean
        Search the folder for zarr stores (zarr=True) or netcdfs (zarr=False)?
    """
    # Read csv from local drive
    # THIS FILE IS NOT CURRENT
    # csv_filepath_local = "temp_clean_all_station_list.csv"
    # full_df = pd.read_csv(csv_filepath_local).loc[:,['era-id','network']]

    # Read csv from s3
    csv_filepath_s3 = (
        "s3://wecc-historical-wx/2_clean_wx/temp_clean_all_station_list.csv"
    )
    full_df = pd.read_csv(csv_filepath_s3).loc[:, ["era-id", "network"]]

    # Add path info as new columns
    full_df["rawdir"] = full_df["network"].apply(lambda row: "1_raw_wx/{}/".format(row))
    full_df["cleandir"] = full_df["network"].apply(
        lambda row: "2_clean_wx/{}/".format(row)
    )
    full_df["qaqcdir"] = full_df["network"].apply(
        lambda row: "3_qaqc_wx/{}/".format(row)
    )
    full_df["mergedir"] = full_df["network"].apply(
        lambda row: "4_merge_wx/{}/".format(row)
    )

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

    # Setting up the QAQC training station list to match temp_clean_all_station_list
    if network == "TRAINING":
        logger.info("Using training station list!")
        df = pd.read_csv("qaqc_training_station_list.csv")
        df["rawdir"] = df["network"].apply(lambda row: "1_raw_wx/{}/".format(row))
        df["cleandir"] = df["network"].apply(lambda row: "2_clean_wx/{}/".format(row))
        df["qaqcdir"] = df["network"].apply(lambda row: "3_qaqc_wx/{}/".format(row))
        df["mergedir"] = df["network"].apply(lambda row: "4_merge_wx/{}/".format(row))
        df["key"] = df.apply(
            lambda row: row["cleandir"] + row["era-id"] + ".nc", axis=1
        )
        df["exist"] = np.zeros(len(df)).astype("bool")

    # If it's a network (not training) run, keep it fast by only checking that network files on s3
    else:
        df = full_df.copy()[
            full_df["network"] == network
        ]  # To use the full dataset for specific sample stations

    # subset for specific network
    for n in df["network"].unique():
        ind = df["network"] == n
        df.loc[ind, "exist"] = file_on_s3(df[ind], zarr=zarr)
    df = df[df["exist"]]

    # If it's a network (not training) run, return df as is
    if network != "TRAINING":
        return df
    else:
        df["file_size"] = df["key"].apply(
            lambda row: s3_cl.head_object(Bucket=bucket_name, Key=row)["ContentLength"]
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


# ----------------------------------------------------------------------------
## Assign ds attributes and save
def process_output_ds(
    df,
    attrs,
    var_attrs,
    network,
    timestamp,
    station,
    qaqcdir,
    errors,
    end_api,
    zarr,
    verbose=False,
    local=False,
):
    """
    DOCUMENTATION NEEDED
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
            if var in list(
                ds.data_vars.keys()
            ):  # Only if var was originally present in dataset
                if "ancillary_variables" in list(ds[var].attrs.keys()):
                    ds[var].attrs["ancillary_variables"] = ds[var].attrs[
                        "ancillary_variables"
                    ] + ", {}".format(eraqc_var)
                else:
                    ds[var].attrs["ancillary_variables"] = "{}".format(eraqc_var)
    # Overwrite file title
    ds = ds.assign_attrs(title=network + " quality controlled")

    # Append qaqc to the file history and comments (https://docs.unidata.ucar.edu/netcdf-c/current/attribute_conventions.html)
    ds.attrs["history"] = ds.attrs[
        "history"
    ] + " \nALLNETWORKS_qaqc.py script run on {} UTC".format(timestamp)
    ds.attrs["comment"] = ds.attrs[
        "comment"
    ] + " \nAn intermediate data product: subject to cleaning but may not be subject to full QA/QC processing.".format(
        timestamp
    )
    # # Flag meaninng attribute
    # ds = ds.assign_attrs(flags_meaning = flags_attrs)
    # --------------------------------------------------------
    # TO DO:
    # Add metadata to `_eraqc` variables

    # --------------------------------------------------------

    # Write station file to netcdf format
    try:
        if zarr == False:
            filename = station + ".nc"  # Make file name
        elif zarr == True:
            filename = station + ".zarr"
        filepath = qaqcdir + filename  # Writes file path

        tmpFile = tempfile.NamedTemporaryFile(
            dir="./temp/", prefix="_" + station, suffix=".nc", delete=False
        )

        # Push file to AWS with correct file name
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            ds.to_netcdf(tmpFile.name)  # Save station file.
        logger.info(
            "Saving/pushing {0} with dims {1} to {2}".format(
                filename, ds.dims, bucket_name + "/" + qaqcdir
            ),
        )
        if zarr == False:  # Upload as netcdf
            s3.Bucket(bucket_name).upload_file(tmpFile.name, filepath)
        elif zarr == True:
            filepath_s3 = "s3://{0}/{1}{2}".format(bucket_name, qaqcdir, filename)
            ds.to_zarr(
                filepath_s3,
                consolidated=True,  # https://docs.xarray.dev/en/stable/internals/zarr-encoding-spec.html
                mode="w",  # Write & overwrite if file with same name exists already
            )
        logger.info(
            "Done saving/pushing file to AWS. Ellapsed time: {:.2f} s.".format(
                time.time() - t0
            ),
        )
        ds.close()
        del ds

        if local:
            t0 = time.time()
            logger.info(
                "Saving local file temp/{}.nc".format(station),
            )
            # Write locally
            os.system("mv {} local_qaqced_files/{}.nc".format(tmpFile.name, station))
            logger.info(
                "Done saving local file. Ellapsed time: {:.2f} s.".format(
                    time.time() - t0
                ),
            )
        else:
            os.system("rm {}".format(tmpFile.name))

    except Exception as e:
        logger.info(
            "netCDF writing failed for {} with Error: {}".format(filename, e),
        )
        errors = print_qaqc_failed(
            errors,
            filename,
            end_api,
            message="Error saving ds as .nc file to AWS bucket: {}".format(e),
            test="process_output_ds",
            verbose=verbose,
        )
        ds.close()
        del ds

        return


# --------------------------------------------------------------------------------
## xarray ds for a station to pandas df in the format needed for the pipeline
def qaqc_ds_to_df(ds, verbose=False):
    ## Add qc_flag variable for all variables, including elevation;
    ## defaulting to nan for fill value that will be replaced with qc flag

    for key, val in ds.variables.items():
        if val.dtype == object:
            if key == "station":
                if str in [type(v) for v in ds[key].values]:
                    ds[key] = ds[key].astype(str)
            else:
                if str in [type(v) for v in ds.isel(station=0)[key].values]:
                    ds[key] = ds[key].astype(str)

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
        "anemometer_height_m",
        "thermometer_height_m",
    ]  # lat, lon have different qc check

    raw_qc_vars = []  # qc_variable for each data variable, will vary station to station
    era_qc_vars = []  # our ERA qc variable
    old_era_qc_vars = []  # our ERA qc variable

    for var in ds.data_vars:
        if "q_code" in var:
            raw_qc_vars.append(
                var
            )  # raw qc variable, need to keep for comparison, then drop
        if "_qc" in var:
            raw_qc_vars.append(
                var
            )  # raw qc variables, need to keep for comparison, then drop
        if "_eraqc" in var:
            era_qc_vars.append(
                var
            )  # raw qc variables, need to keep for comparison, then drop
            old_era_qc_vars.append(var)

    n_qc = len(era_qc_vars)

    for var in ds.data_vars:
        if var not in exclude_qaqc and var not in raw_qc_vars and "_eraqc" not in var:
            qc_var = var + "_eraqc"  # variable/column label

            # if qaqc var does not exist, adds new variable in shape of original variable with designated nan fill value
            if qc_var not in era_qc_vars:
                ds = ds.assign({qc_var: xr.ones_like(ds[var]) * np.nan})
                era_qc_vars.append(qc_var)
                logger.info(
                    "nans created for {}".format(qc_var),
                )
                ds = ds.assign({qc_var: xr.ones_like(ds[var]) * np.nan})
                era_qc_vars.append(qc_var)

    logger.info(
        "{} created era_qc variables".format(len(era_qc_vars) - len(old_era_qc_vars))
    )
    logger.info(old_era_qc_vars)
    logger.info(era_qc_vars)
    if len(era_qc_vars) != n_qc:
        logger.info("{}".format(np.setdiff1d(old_era_qc_vars, era_qc_vars)))
    exit
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
        finally:
            pass
    if "thermometer_height_m" not in df.columns:
        try:
            df["thermometer_height_m"] = (
                np.ones(ds["time"].shape) * ds.thermometer_height_m
            )
        except:
            logger.info("Filling thermometer_height_m with NaN.")
            df["thermometer_height_m"] = np.ones(len(df)) * np.nan
        finally:
            pass

    # De-duplicate time axis
    df = df[~df.index.duplicated()].sort_index()

    # Save station/time multiindex
    MultiIndex = df.index
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

    return df, MultiIndex, attrs, var_attrs, era_qc_vars


# ----------------------------------------------------------------------------
## Run full QA/QC pipeline
def run_qaqc_pipeline(
    ds,
    network,
    file_name,
    errors,
    station,
    end_api,
    rad_scheme,
    verbose=False,
    local=False,
    log_file=log_file,
):
    """ """
    # Convert from xarray ds to pandas df in the format needed for qaqc pipeline
    df, MultiIndex, attrs, var_attrs, era_qc_vars = qaqc_ds_to_df(ds, verbose=verbose)

    # Close ds file, netCDF,HDF5 unclosed files can sometimes cause issues during the mpi4py cleanup phase.
    ds.close()
    del ds

    ##########################################################
    ## QAQC Functions
    # Order of operations
    # Part 1a: Whole station checks - if failure, entire station does not proceed through QA/QC
    # Part 1b: Whole station checks - if failure, entire station does proceed through QA/QC
    # Part 2: Logic checks
    # Part 3: Distribution & time series checks

    # =========================================================
    ## Part 1a: Whole station checks - if failure, entire station does not proceed through QA/QC

    t0 = time.time()
    logger.info("QA/QC whole station tests")
    # ---------------------------------------------------------
    ## Missing values -- does not proceed through qaqc if failure
    stn_to_qaqc = df.copy()  # Need to define before qaqc_pipeline, in case
    new_df = qaqc_missing_vals(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="has an unchecked missing value",
            test="qaqc_missing_vals",
            verbose=verbose,
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_missing_vals")

    # ### DEBUG - REMOVE
    # return stn_to_qaqc

    # ---------------------------------------------------------
    ## Lat-lon -- does not proceed through qaqc if failure
    new_df = qaqc_missing_latlon(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="missing lat-lon",
            test="qaqc_missing_latlon",
            verbose=verbose,
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_missing_latlon")

    # ---------------------------------------------------------
    ## Within WECC -- does not proceed through qaqc if failure
    new_df = qaqc_within_wecc(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="lat-lon is out of range for WECC",
            test="qaqc_within_wecc",
            verbose=verbose,
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_within_wecc")

    # ---------------------------------------------------------
    ## Elevation -- if DEM in-filling fails, does not proceed through qaqc
    new_df = qaqc_elev_infill(
        stn_to_qaqc, verbose=verbose
    )  # nan infilling must be before range check
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="DEM in-filling failed",
            test="DEM in-filling, may not mean station does not pass qa/qc -- check",
            verbose=verbose,
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_elev_infill")

    # ---------------------------------------------------------
    ## Elevation -- range within WECC
    new_df = qaqc_elev_range(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="elevation out of range for WECC",
            test="qaqc_elev_range",
            verbose=verbose,
        )
        return [None] * 4  # whole station failure, skip to next station
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_elev_range")

    # =========================================================
    ## Part 1b: Whole station checks - if failure, entire station does proceed through QA/QC
    # ---------------------------------------------------------
    ## Pressure units fix (temporary)
    new_df = qaqc_pressure_units_fix(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with world record check",
            test="qaqc_pressure_units_fix",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_pressure_units_fix",
        )

    # ---------------------------------------------------------
    ## World record checks: air temperature, dewpoint, wind, pressure
    new_df = qaqc_world_record(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with world record check",
            test="qaqc_world_record",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_world_record")

    logger.info(
        "Done whole station tests, Ellapsed time: {:.2f} s.\n".format(time.time() - t0),
    )
    # =========================================================
    ## Part 2: Variable logic checks

    t0 = time.time()
    logger.info("QA/QC logic checks")
    # ---------------------------------------------------------
    ## dew point temp cannot exceed air temperature
    new_df = qaqc_crossvar_logic_tdps_to_tas_supersat(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with temperature cross-variable logic check",
            test="qaqc_crossvar_logic_tdps_to_tas_supersat",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_crossvar_logic_tdps_to_tas_supersat",
        )

    # ---------------------------------------------------------
    ## dew point temp cannot exceed air temperature (wet bulb drying)
    new_df = qaqc_crossvar_logic_tdps_to_tas_wetbulb(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with temperature cross-variable logic check",
            test="qaqc_crossvar_logic_tdps_to_tas_wetbulb",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_crossvar_logic_tdps_to_tas_wetbulb",
        )

    # ---------------------------------------------------------
    ## precipitation is not negative
    new_df = qaqc_precip_logic_nonegvals(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with negative precipitation values",
            test="qaqc_precip_logic_nonegvals",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_precip_logic_nonegvals",
        )

    # ---------------------------------------------------------
    ## precipitation duration logic
    new_df = qaqc_precip_logic_accum_amounts(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with precip duration logic check",
            test="qaqc_precip_logic_accum_amounts",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_precip_logic_accum_amounts",
        )

    # ---------------------------------------------------------
    ## wind direction should be 0 when wind speed is also 0
    new_df = qaqc_crossvar_logic_calm_wind_dir(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with wind cross-variable logic check",
            test="qaqc_crossvar_logic_calm_wind_dir",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_crossvar_logic_calm_wind_dir",
        )

    logger.info(
        "Done logic checks, Ellapsed time: {:.2f} s.\n".format(time.time() - t0),
    )
    # =========================================================
    ## Part 3: Distribution and timeseries checks - order matters!
    # buoy check
    # frequent values check
    # distributional check (unusual gaps)
    # climatological outliers check
    # unusual streaks check
    # unusual large jumps check (spike)
    ####
    # ---------------------------------------------------------
    ## Buoys with known issues with specific qaqc flags
    ## NDBC and MARITIME only
    if network == "MARITIME" or network == "NDBC":
        t0 = time.time()
        logger.info("QA/QC bouy check")

        new_df = spurious_buoy_check(stn_to_qaqc, era_qc_vars, verbose=verbose)
        if new_df is None:
            errors = print_qaqc_failed(
                errors,
                station,
                end_api,
                message="Flagging problematic buoy issue",
                test="spurious_buoy_check",
                verbose=verbose,
            )
        else:
            stn_to_qaqc = new_df
            logger.info(
                "pass spurious_buoy_check",
            )

        logger.info(
            "Done QA/QC bouy check, Ellapsed time: {:.2f} s.\n".format(
                time.time() - t0
            ),
        )
    # ---------------------------------------------------------
    # frequent values
    t0 = time.time()
    logger.info("QA/QC frequent values")

    new_df = qaqc_frequent_vals(stn_to_qaqc, rad_scheme=rad_scheme, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with frequent values function",
            test="qaqc_frequent_vals",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_frequent_vals")

    logger.info(
        "Done QA/QC frequent values, Ellapsed time: {:.2f} s.\n".format(
            time.time() - t0
        ),
    )
    # ---------------------------------------------------------
    # distribution / unusual gaps
    t0 = time.time()
    logger.info("QA/QC unusual gaps")

    new_df = qaqc_unusual_gaps(stn_to_qaqc, verbose=verbose, local=local)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with unusual gap distribution function",
            test="qaqc_unusual_gaps",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info("pass qaqc_unusual_gaps")

    logger.info(
        "Done QA/QC unusual gaps, Ellapsed time: {:.2f} s.\n".format(time.time() - t0),
    )
    # ---------------------------------------------------------
    # climatological outliers
    t0 = time.time()
    logger.info("QA/QC climatological outliers")

    new_df = qaqc_climatological_outlier(stn_to_qaqc, verbose=verbose)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with climatological outlier check",
            test="qaqc_climatological_outlier",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_climatological_outlier",
        )

    logger.info(
        "Done QA/QC climatological outliers, Ellapsed time: {:.2f} s.\n".format(
            time.time() - t0
        ),
    )
    # ---------------------------------------------------------
    # unusual streaks (repeated values)
    t0 = time.time()
    logger.info("QA/QC unsual repeated streaks")

    new_df = qaqc_unusual_repeated_streaks(stn_to_qaqc, verbose=verbose, local=local)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with unusual streaks (repeated values) check",
            test="qaqc_unusual_repeated_streaks",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_unusual_repeated_streaks",
        )

    logger.info(
        "Done QA/QC unsual repeated streaks, Ellapsed time: {:.2f} s.\n".format(
            time.time() - t0
        ),
    )
    # ---------------------------------------------------------
    # unusual large jumps (spikes)
    t0 = time.time()
    logger.info("QA/QC unsual large jumps")

    new_df = qaqc_unusual_large_jumps(stn_to_qaqc, verbose=verbose, local=local)
    if new_df is None:
        errors = print_qaqc_failed(
            errors,
            station,
            end_api,
            message="Flagging problem with unusual large jumps (spike check) check",
            test="qaqc_unusual_large_jumps",
            verbose=verbose,
        )
    else:
        stn_to_qaqc = new_df
        logger.info(
            "pass qaqc_unusual_large_jumps",
        )

    logger.info(
        "Done QA/QC unsual large jumps, Ellapsed time: {:.2f} s.\n".format(
            time.time() - t0
        ),
    )

    ## END QA/QC ASSESSMENT
    # =========================================================
    # Re-index to original time/station values

    # Calculate flag coverage per variable
    logger.info(
        "Summary of QA/QC flags set per variable",
    )
    flag_summary(stn_to_qaqc, verbose=verbose, local=local)

    stn_to_qaqc = stn_to_qaqc.set_index(MultiIndex).drop(
        columns=["time", "hour", "day", "month", "year", "date", "station"]
    )

    # Sort by time and remove any overlapping timesteps
    # TODO: Is this necessary? Probably done in the cleaning step
    # Check back to see if this can or needs to be removed
    stn_to_qaqc = stn_to_qaqc[~stn_to_qaqc.index.duplicated()].sort_index()

    return stn_to_qaqc, attrs, var_attrs, era_qc_vars


# ==============================================================================
## Function: Conducts whole station qa/qc checks (lat-lon, within WECC, elevation)
def whole_station_qaqc(
    network,
    cleandir,
    qaqcdir,
    rad_scheme,
    zarr,
    verbose=False,
    local=False,
    sample=None,
):
    """
    -----------------------------------
    for station in stations: # full run
    -----------------------------------
    """
    smpi = simpleMPI()

    specific_sample = None
    # ------------------------------------------
    # How to run on a specific station
    # Uncomment "specific_sample" and input desired station id as a list of strings
    # Example: ["ASOSAWOS_74948400395"]
    # specific_sample = ["ASOSAWOS_74948400395", "ASOSAWOS_74509023244", "ASOSAWOS_72494523293"]
    # ------------------------------------------

    # Read in network files
    # This needs to be done only in rank 0, otherwise it gets run by every thread and will overwrite results
    if smpi.rank == 0:
        try:
            files_df = read_network_files(network, zarr)
        except Exception as e:
            errors = print_qaqc_failed(
                errors,
                station="Whole network",
                end_api=end_api,
                message="Error in whole network:",
                test=e,
            )
        # import pdb; pdb.set_trace()
        # # DEBUG -- REMOVE
        # files_df_new = files_df.copy()

        # When "sample" argument is passed to ALLNETWORKS, implements a smaller subset to test
        # Subsetting for a specific set of stations in a single network
        if specific_sample:
            # logger.info(f"Running on specific stations: {specific_sample}")
            stations_sample = specific_sample

        # "all" for no restrictions on sample size
        elif sample == "all":
            stations_sample = list(files_df["era-id"].values)

        # DOCUMENTATION NEEDED
        elif all(char.isnumeric() for char in sample):
            nSample = int(sample)
            files_df = files_df.sample(nSample)
            stations_sample = list(files_df["era-id"])

        # DOCUMENTATION NEEDED
        else:
            files_df = files_df[files_df["era-id"] == sample]
            if len(files_df) == 0:
                logger.info(
                    "Sample station '{}' not in network/stations_df. Please double-check names".format(
                        sample
                    ),
                )
                exit()
                stations_sample = list(files_df["era-id"])
            stations_sample = [sample]

        # print(stations_sample)
        # print(files_df)
        # raise
        # ### DEBUG - REMOVE
        # files_df = files_df_new[files_df_new["era-id"] == "RAWS_AATC1"]
        # # print(files_df)
        # stations_sample = ["RAWS_AATC1"]

        print(
            "Running a sample of {} files on {} network \n Stations: {}".format(
                len(stations_sample), network, stations_sample
            ),
            flush=True,
        )
    else:
        stations_sample = None
        files_df = None

    # DOCUMENTATION NEEDED
    files_df = smpi.comm.bcast(files_df, root=0)
    stations_sample = smpi.comm.bcast(stations_sample, root=0)
    stations_sample_scatter = smpi.scatterList(stations_sample)

    # Loop over stations

    for station in stations_sample_scatter:
        try:
            # ----------------------------------------------------------------------------
            # Set up error handling
            errors, end_api, timestamp = setup_error_handling()

            # ----------------------------------------------------------------------------
            # Set log file
            # DOCUMENTATION NEEDED
            ts = datetime.datetime.utcnow().strftime("%m-%d-%Y")
            # for station in stations_sample:
            #     log_fname = "qaqc_logs/qaqc_{}.{}.log".format(station, ts)
            #     # Initialize the logger with the specific log file name
            #     logger = setup_logger(log_file=log_fname, verbose=verbose)

            log_fname = "qaqc_logs/qaqc_{}.{}.log".format(station, ts)
            # Initialize the logger with the specific log file name
            logger = setup_logger(log_file=log_fname, verbose=verbose)

            # ----------------------------------------------------------------------------
            file_name = files_df.loc[files_df["era-id"] == station, "key"].values[0]
            qaqcdir = files_df.loc[files_df["era-id"] == station, "qaqcdir"].values[0]
            network_ds = files_df.loc[files_df["era-id"] == station, "network"].values[
                0
            ]

            ###################################################################################################
            ## The file_df dataframe must have already checked if file exist in clean directory
            # if file_name not in files: # dont run qa/qc on a station that isn't cleaned
            #     logger.info("{} was not cleaned - skipping qa/qc".format(station))
            #     message = "No cleaned data for this station, does not proceed to qa/qc: see cleaned station list for reason"
            #     errors = print_qaqc_failed(errors, station="Whole network", end_api=end_api,
            #                                message=message, test="whole_station_qaqc")
            #     continue
            # else:
            ## The file_df dataframe must have already checked if file exist in clean directory
            ###################################################################################################
            T0 = time.time()
            # logger.info('Running QA/QC on: {}\n'.format(station)) # testing

            # =====================================================================================
            # Testing speed-up re-order in case file is locally found
            # TODO: DELETE LOCAL READING FOR FINAL VERSION
            fs = s3fs.S3FileSystem()
            aws_url = "s3://wecc-historical-wx/" + file_name
            t0 = time.time()
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=RuntimeWarning)
                    ds = xr.open_dataset("Train_Files/{}.nc".format(station)).load()
            except:
                if zarr == False:  # Read netcdf file
                    with fs.open(aws_url) as fileObj:
                        try:
                            logger.info("Reading {}".format(aws_url))
                            with warnings.catch_warnings():
                                warnings.filterwarnings(
                                    "ignore", category=RuntimeWarning
                                )
                                ds = xr.open_dataset(fileObj).load()
                        except Exception as e:
                            logger.info(
                                "{} did not pass QA/QC because the file could not be opened and/or found in AWS - station not saved.".format(
                                    station
                                )
                            )
                elif zarr == True:  # Or, read zarr
                    try:
                        ds = xr.open_zarr(aws_url)
                    except Exception as e:
                        logger.info(
                            "{} did not pass QA/QC because the file could not be opened and/or found in AWS - station not saved.".format(
                                station
                            ),
                        )
            # Testing speed-up re-order in case file is locally found
            # =====================================================================================

            try:
                # TODO:
                # Same issue than in the pipeline:
                # Probably not needed to drop time duplicates here, if they were properly
                # dropped in the cleaning process?
                # Drop time duplicates

                # There are stations without time/station dimensions
                if ("station" in list(ds.dims.keys())) and (
                    "time" in list(ds.dims.keys())
                ):
                    ds = ds.drop_duplicates(dim="time")
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

                logger.info(
                    f"Done reading. Ellapsed time: {time.time() - t0} s.\n",
                )

                # CHECK THE ENGINE HERE:
                # setting to default which operates on best with dependencies, previously 'h5netcdf'

                # Run full QA/QC pipeline
                logger.info(
                    "Running QA/QC on: {}\n".format(station),
                )
                df, attrs, var_attrs, era_qc_vars = run_qaqc_pipeline(
                    ds,
                    network_ds,
                    file_name,
                    errors,
                    station,
                    end_api,
                    rad_scheme,
                    verbose=verbose,
                    local=local,
                )

                ## Assign ds attributes and save .nc file
                if df is not None:
                    t0 = time.time()

                    process_output_ds(
                        df,
                        attrs,
                        var_attrs,
                        network_ds,
                        timestamp,
                        station,
                        qaqcdir,
                        errors,
                        end_api,
                        zarr,
                        verbose=verbose,
                        local=local,
                    )

            except Exception as e:
                logger.info(
                    "run_qaqc_pipeline failed with error: {}".format(e),
                )
                errors = print_qaqc_failed(
                    errors,
                    station,
                    end_api,
                    message="run_qaqc_pipeline failed with error: {}".format(e),
                    test="run_qaqc_pipeline",
                    verbose=verbose,
                )

        except Exception as e:
            logger.info(
                "QAQC failed\n\n{}\n{}\n\n".format(station, e),
            )
        # Write errors to csv
        finally:
            # pass
            errors = pd.DataFrame(errors)
            csv_buffer = StringIO()
            errors.to_csv(csv_buffer)
            content = csv_buffer.getvalue()

            # Done with station qaqc
            logger.info(
                "Done full QAQC for {}. Ellapsed time: {:.2f} s.\n".format(
                    station, time.time() - T0
                ),
            )

            # Make sure error files save to correct directory
            s3_cl.put_object(
                Bucket=bucket_name,
                Body=content,
                Key=qaqcdir + "qaqc_errs/errors_{}_{}.csv".format(station, end_api),
            )
            # Print error file location
            logger.info(
                "errors_{0}_{1}.csv saved to {2}qaqc_errs/\n".format(
                    network_ds, end_api, bucket_name + "/" + qaqcdir
                ),
            )

            # Save log file to s3 bucket
            logger.info(
                "Saving log file to s3://{0}/{1}{2}\n".format(
                    bucket_name, qaqcdir, log_fname
                ),
            )
            s3.Bucket(bucket_name).upload_file(log_fname, f"{qaqcdir}{log_fname}")

            # Close logging handlers manually
            for handler in logger.handlers:
                handler.close()
                logger.removeHandler(handler)

    return None
