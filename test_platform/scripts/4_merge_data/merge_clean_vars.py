"""
merge_clean_vars.py

Drop unwanted variables from station-level climate data.

This script performs a key preprocessing step in the merge pipeline:
1. Filtering out unnecessary or intermediate variables (e.g., QAQC flags, method indicators),
   while preserving essential columns for downstream analysis.

Functions
---------
- filter_columns(df, logger): Filters a DataFrame's columns based on inclusion and exclusion rules.

"""

import pandas as pd
import inspect
import logging


def filter_columns(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Filters and DataFrame columns based on variable inclusion, exclusion,
    and QA/QC suffix conventions.

    Rules
    -----
    1. Keep only columns that contain any substring in `desired_variables`.
    2. Exclude columns that contain any substring in `delete_variables`.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to process.
    logger : logging.Logger
        Logger for tracking function execution.

    Returns
    -------
    pd.DataFrame
        A DataFrame with filtered and reordered columns.
    """

    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    try:
        # Define variable inclusion and exclusion criteria
        desired_variables = [
            "ps",
            "tas",
            "tdps",
            "pr",
            "hurs",
            "rsds",
            "sfcWind",
            "elevation",
        ]

        delete_variables = [
            "qaqc_process",
            "pr_depth",
            "PREC_flag",
            "rsds_flag",
            "_qc",
            "duration",
            "method",
            "accum_pr",
        ]

        # Always-include variables, regardless of filtering rules
        keep_always = ["station", "time", "lat", "lon"]

        # Filter columns:
        # - Include if the column name contains any substring in desired_variables
        # - Exclude if it contains any substring in delete_variables
        # - Always include if it's in keep_always
        filtered_cols = []
        for col in df.columns:
            is_desired = any(substr in col for substr in desired_variables)
            is_deleted = any(substr in col for substr in delete_variables)
            if (is_desired and not is_deleted) or col in keep_always:
                filtered_cols.append(col)

        # Subset the DataFrame to filtered columns
        df = df[filtered_cols]

        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

        return df

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e
