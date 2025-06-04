"""
merge_eraqc_counts.py

Generate and exports CSVs with counts of unique QAQC flag values per variable, in native
and hourly timesteps. These counts are used to produce QAQC flag statistics for the QAQC success report.

Functions
---------
- eraqc_counts_original_timestep: QAQC flag counts in the original timestep
- eraqc_counts_hourly_timestep: QAQC flag counts in the hourly (post standardization) timestep

Intended Use
------------
Import into merge workflows to generate information for the QAQC success report.
"""

import pandas as pd
import logging
import inspect


# -----------------------------------------------------------------------------
def eraqc_counts_original_timestep(
    df: pd.DataFrame, network: str, station: str, logger: logging.Logger
) -> None:
    """
    Generates a dataframe of raw qaqc flag value counts for every variable,
    in their native timestep, before hourly standardization.
    Exports the dataframe as a csv to AWS.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    network: str
        network name
    station: str
        station name
    logger : logging.Logger
        Logger instance for recording messages during processing.

    Returns
    -------
    None
    """
    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    try:
        # identify _eraqc variables
        eraqc_vars = [var for var in df.columns if "_eraqc" in var]

        # filter df for only qaqc columns
        # also replace Nan values with 'no_flag' for two reasons:
        #   1. to enable us to count total observations for the success report
        #   2. to clarify what the Nan value indicates
        df = df[eraqc_vars].fillna("no_flag")

        # generate df of counts of each unique flag for each variable
        # fill all Nan values with 0, since Nan = no observations counted
        flag_counts = df.apply(pd.Series.value_counts).fillna(0)

        # rename columns
        flag_counts.columns = flag_counts.columns.str.replace("_eraqc", "", regex=True)

        # rename index (i.e. eraqc values) and then reset index
        flag_counts = flag_counts.rename_axis("eraqc_flag_values")

        # set all counts to integers, for readability
        flag_counts = flag_counts.astype(int)

        # send file to AWS
        csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/{network}/eraqc_counts/{station}_flag_counts_native_timestep.csv"
        flag_counts.to_csv(csv_s3_filepath, index=True)

        # Update logger
        logger.info(f"Uploaded file to: {csv_s3_filepath}")
        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e


# -----------------------------------------------------------------------------
def eraqc_counts_hourly_timestep(
    df: pd.DataFrame, network: str, station: str, logger: logging.Logger
) -> None:
    """
    Generates a dataframe of raw qaqc flag value counts for every variable, for the hourly
    timestep, after hourly standardization. Includes the total observation count.
    Exports the dataframe as a CSV to AWS.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    network: str
        network name
    station: str
        station name
    logger : logging.Logger
        Logger instance for recording messages during processing.

    Returns
    -------
    None
    """
    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    try:
        # filter out rows that were infilled during hourly standardization
        df = df[df["standardized_infill"] == "n"]

        # identify _eraqc variables
        eraqc_vars = [var for var in df.columns if "_eraqc" in var]

        # filter df for only qaqc columns
        # also replace Nan values with 'no_flag' for two reasons:
        #   1. to enable us to count total observations for the success report
        #   2. to clarify what the Nan value indicates
        df_qaqc = df[eraqc_vars]

        # generate df of counts of each unique flag for each variable
        # fill all Nan values with 0, since Nan = no observations counted
        flag_counts = df_qaqc.apply(
            lambda x: x.str.split(",", expand=True).stack().value_counts()
        ).fillna(0)

        # rename columns
        flag_counts.columns = flag_counts.columns.str.replace("_eraqc", "", regex=True)

        # rename index (i.e. eraqc values) and then reset index
        flag_counts = flag_counts.rename_axis("eraqc_flag_values")

        # replace 'nan' (a string) with 'no_flag', for clarity
        flag_counts = flag_counts.rename(index={"nan": "no_flag"})

        # add row with total observation count
        total_obs_count = len(df)
        flag_counts.loc["total_obs_count"] = [total_obs_count] * flag_counts.shape[1]

        # set all counts to integers, for readability
        flag_counts = flag_counts.astype(int)

        # send file to AWS
        csv_s3_filepath = f"s3://wecc-historical-wx/4_merge_wx/{network}/eraqc_counts/{station}_flag_counts_hourly_standardized.csv"
        flag_counts.to_csv(csv_s3_filepath, index=True)

        # Update logger
        logger.info(f"Uploaded file to: {csv_s3_filepath}")
        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e
