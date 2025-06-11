"""
merge_clean_vars.py

Drop unwanted variables and reorder columns.

This script combines two merge pipeline steps:
1. Dropping unnecessary or intermediate variables such as QAQC flags and method indicators.
2. Reordering the remaining columns so that core variables and their QAQC counterparts are grouped consistently.

Functions
---------
- merge_drop_vars: Removes temporary and QAQC-related variables from the dataset.
- merge_reorder_vars: Reorders the columns in a standardized format for downstream use.

Intended Use
------------
Import into merge workflows to clean and organize station-level climate data.

"""

import pandas as pd
import inspect
import logging


def filter_and_reorder_columns(
    df: pd.DataFrame, logger: logging.Logger
) -> pd.DataFrame:
    """
    Filters and reorders DataFrame columns based on variable inclusion, exclusion,
    and QA/QC suffix conventions.

    Rules
    -----
    1. Keep only columns that contain any substring in `desired_variables`.
    2. Exclude columns that contain any substring in `delete_variables`.
    3. Reorder so that variables ending in '_eraqc' appear last.

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
        # Keep only columns that match any desired_variables substring,
        # and do NOT match any delete_variables substring
        desired_variables = [
            "ps",
            "tas",
            "tdps",
            "pr",
            "hurs",
            "rsds",
            "sfcWind",
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
        df = df[
            [
                col
                for col in df.columns
                if any(substr in col for substr in desired_variables)
                and not any(substr in col for substr in delete_variables)
            ]
        ]

        # Separate columns that end with '_eraqc' and those that don't
        normal_cols = [col for col in df.columns if not col.endswith("_eraqc")]
        eraqc_cols = [col for col in df.columns if col.endswith("_eraqc")]

        # Reorder the DataFrame with _eraqc columns at the end
        df = df[normal_cols + eraqc_cols]

        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

        return df

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e
