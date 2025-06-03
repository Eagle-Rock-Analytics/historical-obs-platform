"""
MERGE_pipeline.py

This script runs the full merge pipeline for a single weather station dataset.
It processes cleaned QA/QC output and applies a sequence of standardization,
homogenization, and export operations to generate a finalized station file.

Input:
    - Cleaned station Zarr file from S3 (from 3_qaqc_wx/ directory)
    - Station metadata CSV from S3

Output:
    - Final merged and standardized Zarr file written to S3 (4_merge_wx/ directory)
    - Logfile documenting processing steps (also uploaded to S3)

Example:
    Run from command line or as part of a larger batch process for all stations.

"""

from datetime import datetime, timedelta, timezone
import time
import inspect
from typing import Dict

import pandas as pd
import xarray as xr
import logging

from merge_log_config import setup_logger, upload_log_to_s3
from merge_hourly_standardization import merge_hourly_standardization
from merge_derive_missing import merge_derive_missing_vars
from merge_clean_vars import merge_reorder_vars, merge_drop_vars
from merge_eraqc_counts import (
    eraqc_counts_original_timestep,
    eraqc_counts_hourly_timestep,
)


def read_station_metadata(s3_path: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Read station metadata from a CSV file located on S3.

    Parameters
    ----------
    s3_path : str
        Full S3 URI to the CSV metadata file.
    logger : logging.Logger
        Logger instance for logging errors.

    Returns
    -------
    pd.DataFrame
        Station metadata as a DataFrame.

    Raises
    ------
    Exception
        If the file cannot be read.
    """
    try:
        return pd.read_csv(s3_path)
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Could not read file: {s3_path}"
        )
        raise e


def validate_station(
    station: str, stations_df: pd.DataFrame, logger: logging.Logger
) -> str:
    """
    Validate that the station exists in the metadata and return its network.

    Parameters
    ----------
    station : str
        Station identifier to check.
    stations_df : pd.DataFrame
        DataFrame containing station metadata.
    logger : logging.Logger
        Logger for reporting errors.

    Returns
    -------
    str
        Network name associated with the station.

    Raises
    ------
    ValueError
        If the station is not found.
    """
    if station not in stations_df["era-id"].values:
        logger.error(f"Station {station} not found in station list.")
        raise ValueError(f"Station {station} not found in station list.")
    return stations_df.loc[stations_df["era-id"] == station, "network"].item()


def read_zarr_dataset(
    bucket: str, qaqc_dir: str, network: str, station: str, logger: logging.Logger
) -> xr.Dataset:
    """
    Load a Zarr dataset from S3.

    Parameters
    ----------
    bucket : str
        Name of the S3 bucket.
    qaqc_dir : str
        Directory within the bucket containing QAQC datasets.
    network : str
        Network name for the station.
    station : str
        Station identifier.
    logger : logging.Logger
        Logger for error reporting.

    Returns
    -------
    xr.Dataset
        Loaded xarray dataset.

    Raises
    ------
    Exception
        If the dataset cannot be opened.
    """
    s3_uri: str = f"s3://{bucket}/{qaqc_dir}/{network}/{station}.zarr/"
    try:
        return xr.open_zarr(s3_uri)
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Could not open Zarr dataset at {s3_uri}"
        )
        raise e


def get_var_attrs(
    ds: xr.Dataset, network: str, logger: logging.Logger
) -> Dict[str, dict]:
    """
    Extracts variable attributes from all data variables in an xarray Dataset.

    For the "ASOSAWOS" network, if the "pr" variable is present, its "units"
    attribute is explicitly set to "mm".

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset containing climate variables.
    network : str
        Network identifier (e.g., "ASOSAWOS").
    logger : logging.Logger
        Logger instance used for logging.

    Returns
    -------
    var_attrs : dict[str, dict]
        A dictionary mapping variable names to their attribute dictionaries.
    """

    try:
        var_attrs = {var: ds[var].attrs.copy() for var in ds.data_vars}

        if network == "ASOSAWOS" and "pr" in ds.data_vars:
            var_attrs["pr"]["units"] = "mm"
            logger.info(
                f"{inspect.currentframe().f_code.co_name}: Set 'pr' units to 'mm' for ASOSAWOS network to resolve error in units attribute."
            )

        return var_attrs

    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to retrieve variable attributes from xr.Dataset."
        )
        raise e


def convert_xr_to_df(
    ds: xr.Dataset, logger: logging.Logger
) -> tuple[pd.DataFrame, pd.MultiIndex]:
    """
    Convert an xarray Dataset into a flat pandas DataFrame, while preserving the original MultiIndex.

    Parameters
    ----------
    ds : xr.Dataset
        Input xarray Dataset.
    logger : logging.Logger
        Logger for error reporting.

    Returns
    -------
    pd.DataFrame
        The flattened DataFrame with index reset.

    Raises
    ------
    Exception
        If the conversion from xarray to pandas DataFrame fails.
    """
    try:
        df = ds.to_dataframe()
        df.reset_index(inplace=True)  # Flatten to remove MultiIndex
        return df
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to convert xarray Dataset to DataFrame."
        )
        raise e


def convert_df_to_xr(
    df: pd.DataFrame, ds_attrs: dict, var_attrs: dict, logger: logging.Logger
) -> xr.Dataset:
    """
    Converts a DataFrame to an xarray Dataset and assigns global and variable attributes.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame to convert.
    ds_attrs : dict
        Global attributes to assign to the resulting Dataset.
    var_attrs : dict
        Dictionary mapping variable names to their attributes.
    logger : logging.Logger
        Logger for error and info reporting.

    Returns
    -------
    xr.Dataset
        Dataset with assigned global and variable-level attributes.

    Raises
    ------
    Exception
        If setting the index, conversion, or attribute assignment fails.
    """
    try:
        df.to_csv("test.csv")
        df.set_index(["station", "time"], inplace=True)
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to set DataFrame MultiIndex with station and time. This is required for conversion from pd.DataFrame --> xr.Dataset object with correct cooridnates."
        )
        raise e

    try:
        ds = df.to_xarray()
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to convert pd.DataFrame to xr.Dataset."
        )
        raise e

    try:
        # Add history into Dataset attributes
        timestamp = datetime.now(timezone.utc).strftime("%m-%d-%Y, %H:%M:%S")
        ds_attrs["history"] += f"\nMERGE_pipeline run on {timestamp} UTC"
        ds = ds.assign_attrs(ds_attrs)

        # Assign attributes for each variable
        for var, attrs in var_attrs.items():
            ds[var] = ds[var].assign_attrs(attrs)
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to assign attributes to final xr.Dataset object."
        )
        raise e

    return ds


def write_zarr_to_s3(
    ds: xr.Dataset,
    bucket_name: str,
    merge_dir: str,
    network: str,
    station: str,
    logger: logging.Logger,
) -> None:
    """
    Writes the xarray Dataset to a Zarr file and uploads to S3.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to write.
    bucket_name : str
        S3 bucket name.
    merge_dir : str
        S3 directory prefix.
    network : str
        Network identifier.
    station : str
        Station identifier.
    logger : logging.Logger
        Logger for status and errors.

    Raises
    ------
    Exception
        If writing to Zarr fails.
    """
    zarr_s3_path = f"s3://{bucket_name}/{merge_dir}/{network}/{station}.zarr"
    try:
        ds.to_zarr(
            zarr_s3_path,
            consolidated=True,
            mode="w",
        )
        logger.info(f"Uploaded Zarr to: {zarr_s3_path}")
    except Exception as e:
        logger.error(
            f"{inspect.currentframe().f_code.co_name}: Failed to write dataset to Zarr format on S3"
        )
        raise e


def run_merge_one_station(
    station: str,
    verbose: bool = False,
) -> None:
    """
    Main entry point for running the merge pipeline for a single station.

    Parameters
    ----------
    station : str
        Unique identifier for the weather station (e.g., "LOXWFO_CBGC1").
    verbose : bool, optional
        If True, enables verbose logging. Default is False.

    Returns
    -------
    None

    """

    bucket_name = "wecc-historical-wx"
    stations_csv_path = f"s3://{bucket_name}/2_clean_wx/temp_clean_all_station_list.csv"
    qaqc_dir = "3_qaqc_wx"
    merge_dir = "4_merge_wx"

    # Log start time
    start_time = time.time()

    ## ======== SETUP ========

    # Set up logger
    logger, log_filepath = setup_logger(station, verbose=verbose)

    # Load station metadata
    stations_df = read_station_metadata(stations_csv_path, logger)

    # Validate station and get network name
    network_name = validate_station(station, stations_df, logger)

    try:

        ## ======== READ IN AND REFORMAT DATA ========

        # Load Zarr dataset from S3
        ds = read_zarr_dataset(bucket_name, qaqc_dir, network_name, station, logger)

        # Get variable attributes from dataset
        var_attrs = get_var_attrs(ds, network_name, logger)

        # Convert dataset to DataFrame
        df = convert_xr_to_df(ds, logger)

        # ======== MERGE FUNCTIONS ========

        # Part 1: Construct and export table of raw QAQC counts per variable
        # For success report
        eraqc_counts_original_timestep(df, network_name, station, logger)

        # Part 2: Derive any missing variables
        df, var_attrs = merge_derive_missing_vars(df, var_attrs, logger)

        # Part 3: Standardize sub-hourly observations to hourly
        df, var_attrs = merge_hourly_standardization(df, var_attrs, logger)

        # Part 3b: Construct and export table of raw QAQC counts per variable post-hourly standardization
        eraqc_counts_hourly_timestep(df, network_name, station, logger)

        # Part 4: Drops raw _qc variables (DECISION TO MAKE) or provide code to filter
        df, var_attrs = merge_drop_vars(df, var_attrs)

        # Part 5: Re-orders variables into final preferred order
        df = merge_reorder_vars(df)

        # ======== CLEANUP & UPLOAD DATA TO S3 ========

        # Convert the cleaned DataFrame to an xarray.Dataset and assign global + variable-level metadata
        ds_merged = convert_df_to_xr(df, ds.attrs, var_attrs, logger)

        # Write the xarray Dataset as a Zarr file to the specified S3 path
        write_zarr_to_s3(
            ds_merged, bucket_name, merge_dir, network_name, station, logger
        )

        # Done! Print elapsed time
        logger.info(f"Finished processing station: {station}")

    except Exception as e:
        logger.info(f"Error traceback: {type(e).__name__}: {e}")
        logger.info(f"Terminating merge script.")

    finally:  # Even in case of failure, upload logfile to s3

        # Compute elapsed time
        elapsed_time_seconds = time.time() - start_time
        formatted_elapsed = str(timedelta(seconds=round(elapsed_time_seconds)))
        logger.info(f"Elapsed time: {formatted_elapsed}")

        # Save log file to s3 bucket
        upload_log_to_s3(bucket_name, merge_dir, network_name, log_filepath, logger)

        print("Script complete.")
        print(f"Elapsed time: {formatted_elapsed}")
