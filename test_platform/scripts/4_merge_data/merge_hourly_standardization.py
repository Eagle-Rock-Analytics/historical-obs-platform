"""
This is a script where meteorological variables are resampled to an hourly timestep according to standard conventions.
"""

## Import Libraries
from functools import reduce
import boto3
import numpy as np
import pandas as pd

# New logger function
from merge_log_config import logger


# -----------------------------------------------------------------------------
def merge_hourly_standardization(
    df: pd.DataFrame, var_attrs: dict
) -> tuple[pd.DataFrame, dict]:
    """Resamples meteorological variables to hourly timestep according to standard conventions.

    Parameters
    -----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    var_attrs: library
        attributes for sub-hourly variables

    Returns
    -------
    df : pd.DataFrame | None
        returns a dataframe with all columns resampled to one hour (column name retained)
    var_attrs : dict | None
        returns variable attributes dictionary updated to note that sub-hourly variables are now hourly

    Notes
    -----
    Rules:
    1. Top of the hour: take the first value in each hour. Standard convention for temperature, dewpoint, wind speed, direction, relative humidity, air pressure.
    2. Summation across the hour: sum observations within each hour. Standard convention for precipitation and solar radiation.
    3. Constant across the hour: take the first value in each hour. This applied to variables that do not change.
    """

    printf(
        "Running: hourly_standardization",
        verbose=verbose,
        log_file=log_file,
        flush=True,
    )

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
    vars_to_remove = ["qc", "eraqc", "duration", "method", "flag", "depth", "process"]
    qaqc_vars = [
        var
        for var in df.columns
        if not any(True for item in vars_to_remove if item in var)
    ]

    # All variables, necessary for producing columns with hourly counts for each variable
    # all_vars = constant_vars + sum_vars + instant_vars + qaqc_vars

    # Subset the dataframe according to rules
    constant_df = df[[col for col in constant_vars if col in df.columns]]

    qaqc_df = df[[col for col in qaqc_vars if col in df.columns if col != "time"]]
    qaqc_df = qaqc_df.astype(str)
    qaqc_df.insert(0, "time", df["time"])

    sum_df = df[[col for col in sum_vars if col in df.columns]]

    instant_df = df[[col for col in instant_vars if col in df.columns]]

    try:
        # If station does not report any variable, bypass
        if len(df.columns) == 0:
            printf(
                "Empty dataset - bypassing hourly aggregation",
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            return df, var_attrs
        else:
            result_list = []

            # Performing hourly aggregation, only if subset contains more than one (ie 'time') column
            # This is to account for input dataframes that do not contain all subsets of variables defined above.
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
                    lambda x: ",".join(x.unique())
                )  # adding unique flags
                result_list.append(qaqc_result)

            # Aggregating and outputting reduced dataframe
            result = reduce(
                lambda left, right: pd.merge(left, right, on=["time"], how="outer"),
                result_list,
            )

            # Update attributes for sub-hourly variables
            sub_hourly_vars = [i for i in df.columns if "min" in i and "qc" not in i]
            for var in sub_hourly_vars:
                var_attrs[var]["standardization"] = (
                    "{} has been standardized to an hourly timestep, but will retain its original name".format(
                        var
                    )
                )

            return result, var_attrs

    except Exception as e:
        printf(
            "hourly_standardization failed with Exception: {0}".format(e),
            verbose=verbose,
            log_file=log_file,
            flush=True,
        )
        # convert to logger version
        return None
