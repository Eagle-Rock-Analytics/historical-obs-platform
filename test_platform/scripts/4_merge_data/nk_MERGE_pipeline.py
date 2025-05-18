"""MERGE_pipeline.py

Run merge functions for a single station.
"""

from datetime import datetime, timedelta, timezone
import time

import pandas as pd
import xarray as xr
import logging
import boto3

from merge_log_config import init_logger
from merge_hourly_standardization import merge_hourly_standardization


def setup_logger(
    station: str, logs_dir: str = "qaqc_logs", verbose: bool = True
) -> tuple[logging.Logger, str]:
    """
    Initialize a timestamped logger for the current station.

    Parameters
    ----------
    station : str
        Station identifier for log file naming.
    logs_dir : str
        Directory to save logs. Defaults to "qaqc_logs".
    verbose : bool
        If True, also logs to console.

    Returns
    -------
    logger : logging.Logger
        Configured logger instance.
    log_filepath : str
        Full path to the created log file.
    """
    timestamp = datetime.now().strftime("%Y%m%d%H%M")
    log_filename = f"merge_{station}.{timestamp}.log"
    log_filepath = f"{logs_dir}/{log_filename}"

    logger = init_logger(log_filename, logs_dir=logs_dir, verbose=verbose)
    logger.info(f"Starting merge script for station: {station}")

    return logger, log_filepath


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
        logger.error(f"Could not read file: {s3_path}")
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
        logger.error(f"Could not open Zarr dataset at {s3_uri}")
        raise e


def convert_to_dataframe(
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
        logger.error("Failed to convert xarray Dataset to DataFrame.")
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
            "Failed to set DataFrame MultiIndex with station and time. This is required for conversion from pd.DataFrame --> xr.Dataset object with correct cooridnates."
        )
        raise e

    try:
        ds = df.to_xarray()
    except Exception as e:
        logger.error("Failed to convert pd.DataFrame to xr.Dataset.")
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
        logger.error("Failed to assign attributes to final xr.Dataset object.")
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
        logger.error("Failed to write dataset to Zarr format on S3")
        raise e


def upload_log_to_s3(
    bucket_name: str,
    merge_dir: str,
    network: str,
    logfile_fpath: str,
    logger: logging.Logger,
) -> None:
    """
    Upload a local log file to an S3 bucket and log the operation.

    Parameters
    ----------
    bucket_name : str
        Name of the S3 bucket.
    merge_dir : str
        Prefix (directory path) inside the bucket to store the log file.
    network: str
        Name of network. Used to construct upload path for logfile
    logfile_fpath : str
        Path to local logfile
    logger : logging.Logger
        Logger instance for logging actions and errors.

    Raises
    ------
    Exception
        If upload to S3 fails, the original exception is raised.
    """

    # Ensure all log data is flushed to disk before uploading
    close_logger(logger)

    # Construct full S3 URI for logging
    s3_key = f"{merge_dir}/{network}/{logfile_fpath}"
    logfile_s3_uri = f"s3://{bucket_name}/{s3_key}"

    # Upload to s3
    s3 = boto3.resource("s3")
    try:
        s3.Bucket(bucket_name).upload_file(logfile_fpath, s3_key)
    except Exception as e:
        logger.error(f"Failed to upload log file to S3: {e}")
        raise

    logger.info(f"Saved log file to {logfile_s3_uri}")


def close_logger(logger: logging.Logger) -> None:
    """
    Closes all handlers associated with the given logger.

    Parameters
    ----------
    logger : logging.Logger
        The logger instance to close.
    """
    for handler in logger.handlers[:]:  # Copy the list to avoid modification issues
        handler.close()
        logger.removeHandler(handler)


def main() -> None:
    """
    Main entry point for running the merge pipeline for a single station.
    """
    STATION = "ASOSAWOS_69007093217"
    VERBOSE = True

    bucket_name = "wecc-historical-wx"
    stations_csv_path = f"s3://{bucket_name}/2_clean_wx/temp_clean_all_station_list.csv"
    qaqc_dir = "3_qaqc_wx"
    merge_dir = "4_merge_wx"

    # Log start time
    start_time = time.time()

    try:
        # Set up logger
        logger, log_filepath = setup_logger("ASOSAWOS_690070932172", verbose=VERBOSE)

        # Load station metadata
        stations_df: pd.DataFrame = read_station_metadata(stations_csv_path, logger)

        # Validate station and get network name
        network: str = validate_station(STATION, stations_df, logger)

        # Load Zarr dataset from S3
        ds = read_zarr_dataset(bucket_name, qaqc_dir, network, STATION, logger)
        var_attrs = {
            var: ds[var].attrs for var in list(ds.data_vars.keys())
        }  # Attributes from each variable

        # Convert dataset to DataFrame
        df = convert_to_dataframe(ds, logger)

        # Perform hourly standardization
        # df, var_attrs = merge_hourly_standardization(df, var_attrs, logger)

        # Convert the cleaned DataFrame to an xarray.Dataset and assign global + variable-level metadata
        ds_merged = convert_df_to_xr(df, ds.attrs, var_attrs, logger)

        # Write the xarray Dataset as a Zarr file to the specified S3 path
        write_zarr_to_s3(ds_merged, bucket_name, merge_dir, network, STATION, logger)

        # Done! Print elapsed time
        logger.info(f"Finished processing station: {STATION}")

    except Exception as e:
        logger.error(f"Terminating merge script.")

    finally:  # Even in case of failure, upload logfile to s3

        # Compute elapsed time
        elapsed_time_seconds = time.time() - start_time
        formatted_elapsed = str(timedelta(seconds=round(elapsed_time_seconds)))
        logger.info(f"Elapsed time: {formatted_elapsed}")

        # Save log file to s3 bucket
        upload_log_to_s3(bucket_name, merge_dir, network, log_filepath, logger)

        print("Script complete.")
        print(f"Elapsed time: {formatted_elapsed}")


if __name__ == "__main__":
    main()
