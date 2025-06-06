"""merge_hourly_standardization.py

Resamples data to an hourly timestep according to standard conventions.

"""

from functools import reduce
import numpy as np
import pandas as pd
import logging
import inspect


# -----------------------------------------------------------------------------
def qaqc_flag_fcn(flags: str) -> str:
    """
    Used for resampling QAQC flag columns. Ensures that the final standardized dataframe
    does not contain any empty strings by returning 'nan' when given an empty input (i.e. in time gaps).

    Parameters
    -----------
    flags : array_like
        sub-hourly timestep data

    Returns
    -------
    str : final flag value

    """
    if len(flags) == 0:
        return "nan"
    else:
        return ",".join(flags.unique())


# -----------------------------------------------------------------------------
def _modify_infill(df: pd.DataFrame, constant_vars: list) -> pd.DataFrame:
    """
    This function does two things:
    1. Flags rows ('y' under 'standardized_infill') that were infilled by resampling in the hourly standardization
        process, where there were time gaps in the input dataframe. These infilled rows will NOT count towards
        the total observations count when calculating flag rates for the success report
    2. Infills constant variables (ie those in "constant_vars") observations that were left empty because
        they were in a time gap. They are infilled with the first non-nan value of each column, and set to
        np.nan if there are no non-nan values.

    Parameters
    -----------
    df : pd.Dataframe
        hourly standardized dataframe
    constant_vars: list
        variables that are constant throughout time

    Returns
    -------
    df : pd.Dataframe
        dataframe with updates added to rows infilled by hourly standardization

    """
    # Mask for rows where station is None (or np.nan)
    mask = df["station"].isnull()

    # Initialize dict to hold first non-NaN values
    first_valids = {}

    # Populate first_valids only for existing columns
    for col in constant_vars:
        if col not in df.columns or col == "time":
            # skip if constant var not in df cols or if var == "time"
            continue

        first_valids[col] = (
            df[col].dropna().iloc[0] if df[col].notna().any() else np.nan
        )

    # Update values in masked rows for existing columns
    for col, val in first_valids.items():
        df.loc[mask, col] = val

    # Add or update 'standardized_infill' column
    # 'y' indicates that the data was infilled during resampling (ie was not present in the input dataframe)
    # 'n' indicates that the data was not infilled (ie present in the input dataframe)
    df["standardized_infill"] = np.where(mask, "y", "n")

    return df


# -----------------------------------------------------------------------------
def merge_hourly_standardization(
    df: pd.DataFrame, var_attrs: dict, logger: logging.Logger
) -> tuple[pd.DataFrame, dict]:
    """Resamples meteorological variables to hourly timestep according to standard conventions.

    Parameters
    -----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    var_attrs: library
        attributes for sub-hourly variables
    logger : logging.Logger
        Logger instance for recording messages during processing.

    Returns
    -------
    df : pd.DataFrame | None
        returns a dataframe with all columns resampled to one hour (column name retained)
    var_attrs : dict | None
        returns variable attributes dictionary updated to note that sub-hourly variables are now hourly

    Notes
    -----
    Regarding "accum_pr":
    A properly standardized "accum_pr" would keep the last value within each hour. This is not how the variable is standardized
    in this function which, rather, takes the FIRST value within each hour. This choice is inconsequantial, though, since "accum_pr" 
    is dropped in part 4 of the merge pipeline.

    Rules:
    1. Top of the hour: take the first value in each hour. Standard convention for temperature, dewpoint, wind speed, direction, relative humidity, air pressure.
    2. Summation across the hour: sum observations within each hour. Standard convention for precipitation and solar radiation.
    3. Constant across the hour: take the first value in each hour. This applied to variables that do not change.
    """

    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    # Variables that remain constant within each hour
    constant_vars = [
        "time",
        "station",
        "lat",
        "lon",
        "elevation",
        "anemometer_height_m",
        "thermometer_height_m",
    ]

    # Aggregation across hour variables, standard meteorological convention: precipitation and solar radiation
    sum_vars = [
        "time",
        "accum_pr",
        "pr",
        "pr_localmid",
        "pr_24h",
        "pr_1h",
        "pr_15min",
        "pr_5min",
        "rsds",
    ]

    # Top of the hour variables, standard meteorological convention: temperature, dewpoint temperature, pressure, humidity, winds
    instant_vars = [
        "hurs_derived",
        "time",
        "tas",
        "tas_derived",
        "tdps",
        "tdps_derived",
        "ps",
        "psl",
        "ps_altimeter",
        "ps_derived",
        "hurs",
        "sfcWind",
        "sfcWind_dir",
        "accum_pr",  # NOT handled properly here, intentionally - see note in docstring above
    ]

    # QAQC variable suffixes
    # the QAQC variables contain these words (ex: ps_eraqc)
    qaqc_var_suffixes = [
        "qc",
        "eraqc",
        "duration",
        "method",
        "flag",
        "depth",
        "process",
    ]

    # QAQC variables, which we will concatenate within each hour
    qaqc_vars = [
        var for var in df.columns if any(item in var for item in qaqc_var_suffixes)
    ]

    try:
        # Subset the dataframe according to rules
        constant_df = df[[col for col in constant_vars if col in df.columns]]

        qaqc_df = df[[col for col in qaqc_vars if col in df.columns if col != "time"]]
        qaqc_df = qaqc_df.astype(str)
        qaqc_df.insert(0, "time", df["time"])

        sum_df = df[[col for col in sum_vars if col in df.columns]]

        instant_df = df[[col for col in instant_vars if col in df.columns]]

        # Performing hourly aggregation, only if subset contains more than one (ie more than the 'time' time) column
        # This is to account for input dataframes that do not contain ALL subsets of variables defined above - just a subset of them.
        result_list = []
        if len(constant_df.columns) > 1:
            constant_result = constant_df.resample("1h", on="time").first()
            result_list.append(constant_result)

        if len(instant_df.columns) > 1:
            instant_result = instant_df.resample("1h", on="time").first()
            result_list.append(instant_result)

        if len(sum_df.columns) > 1:
            sum_result = sum_df.resample("1h", on="time").apply(
                lambda x: np.nan if x.isna().all() else x.sum(skipna=True)
            )
            result_list.append(sum_result)

        if len(qaqc_df.columns) > 1:
            qaqc_result = qaqc_df.resample("1h", on="time").apply(
                lambda x: qaqc_flag_fcn(x)
            )  # concatenating unique flags
            result_list.append(qaqc_result)

        # Aggregate and output reduced dataframe - this merges all dataframes defined
        # This function sets "time" to the index; reset index to return to original index
        result = reduce(
            lambda left, right: pd.merge(left, right, on=["time"], how="outer"),
            result_list,
        )
        result.reset_index(inplace=True)  # Convert time index --> column

        # Infill constant values and flag rows added through resampling
        result = _modify_infill(result, constant_vars)

        # Update attributes for sub-hourly variables
        sub_hourly_vars = [i for i in df.columns if "min" in i and "qc" not in i]
        for var in sub_hourly_vars:
            var_attrs[var]["standardization"] = (
                "{} has been standardized to an hourly timestep, but will retain its original name".format(
                    var
                )
            )
        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")
        return result, var_attrs

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e
