"""
This script performs the final merge protocols for cleaned and quality controlled station data for ingestion into the Historical Observations Platform, 
and is independent of network. 
Approach:
(1) Derive any missing variables
(2) Standardize sub-hourly observations to hourly
(3) Homogenize ASOSAWOS stations where there are historical jumps
(4) Remove duplicate stations
(5) Re-orders variables into final preferred order
(6) Drops raw _qc variables (DECISION TO MAKE) OR PROVIDE CODE TO FILTER
(7) Exports final station file as a .nc file (or .zarr)

Inputs: QA/QC-processed data for an individual network
Outputs: Merged data for an individual network, priority variables, all times. Organized by station as .nc file.
"""

# Import libraries
import os
import datetime
import xarray as xr
import boto3
import s3fs
from io import BytesIO, StringIO
import time
import tempfile
from merge_log_config import logger

# Import all merge script functions
try:
    from merge_utils import merge_hourly_standardization
except Exception as e:
    logger.debug("Error importing merge script: ".format(e))

# Set up directory to save files temporarily and save timing, if it doesn't already exist.
# TODO: Decide if we also merge all log files into a single one too
dirs = ["./temp/", "./local_merged_files/", "./merge_logs/"]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)

# Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

# Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"

# ----------------------------------------------------------------------------
# Global functions and variables


# POTENTIALLY MOVE INTO UTILS
def setup_error_handling():
    """DOCUMETNATION NEEDED"""

    errors = {"File": [], "Time": [], "Error": []}  # Set up error handling
    end_api = datetime.datetime.now().strftime(
        "%Y%m%d%H%M"
    )  # Set end time to be current time at beginning of download: for error handling csv

    timestamp = datetime.datetime.utcnow().strftime("%m-%d-%Y, %H:%M:%S")

    return errors, end_api, timestamp


# ----------------------------------------------------------------------------
def print_merge_failed(
    errors, station=None, end_api=None, message=None, test=None, verbose=False
):
    """DOCUMETNATION NEEDED"""
    logger.info("{0} {1}, skipping station".format(station, message))
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
    for item in s3.Bucket(bucket_name).objects.filter(Prefix=df["qaqcdir"].iloc[0]):
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
## COULD BE PUT IN UTILS
def merge_ds_to_df(ds, verbose=verbose):
    """Converts xarray ds for a station to pandas df in the format needed for processing.

    Parameters
    ----------
    ds: xr.Dataset
        Data object with information about each network and station
    verbose: boolean
        Flag as to whether to print runtime statements to terminal. Default is False. Set in ALLNETWORKS_merge.py run.

    Returns
    -------
    df: pd.DataFrame
        Table object with information about each network and station
    MultiIndex: pd.DataFrame (I think)
        Original multi-index of station and time, to be used on conversion back to ds
    attrs:
        Save ds attributes to inherent to the final merged file
    var_attrs:
        Save variable attributes to inherent to the final merged file
    """

    # Save attributes to inherent them to the final merged file
    attrs = ds.attrs
    var_attrs = {var: ds[var].attrs for var in list(ds.data_vars.keys())}

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        df = ds.to_dataframe()

    # Save instrumentation heights
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

    return df, MultiIndex, attrs, var_attrs


# ----------------------------------------------------------------------------
# Run full merge pipeline
def run_merge_pipeline(
    ds,
    network,
    file_name,
    errors,
    station,
    end_api,
    verbose=verbose,
    log_file=None,
):
    """Runs all final merge standardization functions, and exports final station file.

    Parameters
    ----------
    ds: xr.Dataset
        Data object with information about each network and station
    network: str
        Network identifier
    file_name: str
        Station filename
    errors: str
        Path to the errors file location -- CHECK
    station: str
        Staiton identifier
    end_api: str
        Script end time, for error handling csv

    Returns
    -------
    None
        This function does not return a value
    """

    # Convert to working dataframe
    df, MultiIndex, attrs, var_attrs = merge_ds_to_df(ds, verbose=verbose)

    # Close ds file, netCDF, HDF5 unclosed files can cause issues during mpi4py run
    ds.close()
    del ds

    # =========================================================
    # Set up timing and logging for script runtime
    t0 = time.time()
    logger.info("Beginning final merge step...")

    # Set up final dataframe, in case
    stn_to_merge = df.copy()

    ##########################################################
    ## Merge Functions: Order of operations
    # Part 1: Derive any missing variables
    # Part 2: Standardize sub-hourly observations to hourly
    # Part 3: Homogenize ASOSAWOS, VALLEYWATER, NDBC stations where there are historical jumps
    # Part 4: Remove duplicate stations
    # Part 5: Re-orders variables into final preferred order
    # Part 6: Drops raw _qc variables (DECISION TO MAKE) OR PROVIDE CODE TO FILTER
    # Part 7: Exports final station file as a .nc file (or .zarr)

    # =========================================================
    # Part 1: Derive any missing variables
    # TODO: Do this only for variables which the station has no sensor for (do not mix observed & calculated values)
    # Will require meteorological formulae -- some in calc_clean.py, some in climakitae?
    # dew point temperature
    # relative humidity
    # Not started

    # ----------------------------------------------------------
    # Part 2: Standardize sub-hourly observations to hourly
    new_df = merge_hourly_standardization(df, verbose=verbose)
    if new_df is None:
        errors = print_merge_failed(
            errors,
            station,
            end_api,
            message="hourly standardization failed",
            test="merge_hourly_standardization",
        )
        return [None]
    else:
        stn_to_merge = new_df
        # Update attributes
        
        logger.info("pass merge_hourly_standardization")

    # ----------------------------------------------------------
    # Part 3: Homogenize ASOSAWOS stations where there are historical jumps
    # In progress -- need to read in csv file of suspect stations
    new_df = merge_concat_jump_stns(df, verbose=verbose)
    if new_df is None:
        errors = print_merge_failed(
            errors,
            station,
            end_api,
            message="station concatenation failed",
            test="merge_concat_jump_stns",
        )
        return [None]
    else:
        stn_to_merge = new_df
        logger.info("pass merge_concat_jump_stns")

    # ----------------------------------------------------------
    # Part 4: Remove duplicate stations
    # TODO:
    # name string matching
    # stations within a certain distance (vs. lat-lon matching)
    # Not started

    # ----------------------------------------------------------
    # Part 5: Re-orders variables into final preferred order
    # TODO:
    # Not started

    # ----------------------------------------------------------
    # Part 6: Drops raw _qc variables (DECISION TO MAKE) OR PROVIDE CODE TO FILTER
    # TODO: Decision needs to be made as to whether we keep raw qc variables and/or eraqc variables in final data product
    # Not started

    # ----------------------------------------------------------
    # Part 7: Exports final station file as a .zarr file (or .nc)
    # AE preference would be zarrs
    # Not started
    # Assign ds attributes and save .zarr
    # process output ds
    # ensure that each variable is the right datatype!!
    # Close and save log file
    # Write errors to csv
    # Make sure error files save to correct directory

    # for testing!
    return None
