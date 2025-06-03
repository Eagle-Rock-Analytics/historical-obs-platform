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


def merge_reorder_vars(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """
    Reorders input dataframe columns

    Rules
    ------
        1.) Non-qaqc variables that start with the strings in "desired_order" come first,
            followed by their associated qaqc variables, followed by all remaining variables

    Parameters
    ------
    df: pd.DataFrame

    Returns
    -------
    df : pd.DataFrame | None
        Returns a dataframe with reordered columns

    Notes
    -------

    """
    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    try: 
        desired_order = [
            "ps",
            "tas",
            "tdps",
            "pr",
            "hurs",
            "rsds",
            "sfcWind",
            "pvp",
            "svp",
        ]

        # Select variables with names that start with those in "desired_order"
        new_order = [
            i for keyword in desired_order for i in df.columns if i.startswith(keyword)
        ]

        # Now split them into qaqc and non-qaqc variables
        qaqc_vars = [i for i in new_order if "qc" in i]
        nonqaqc_vars = [i for i in new_order if i not in qaqc_vars]

        # Now store all remaining columns
        rest_of_vars = [i for i in list(df.columns) if i not in new_order]

        # Generate the complete list of variables, in the correct order
        final_order = nonqaqc_vars + qaqc_vars + rest_of_vars

        # Use that list to reorder the columns in "df"
        df = df[final_order]

        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

        return df
    
    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e


def merge_drop_vars(df: pd.DataFrame, var_attrs: dict, logger: logging.Logger) -> tuple[pd.DataFrame, dict]:
    """
    Keep “_eraqc” vars and drop the following variables
        - qaqc_process
        - pr_duration
        - pr_depth
        - PREC_flag
        - rsds_duration
        - rsds_flag
        - q_code
        - any "_qc" or "method "variable

    Parameters
    ------
    df: pd.DataFrame
        station data
    var_attrs: dict
        variable attributes

    Returns
    -------
    if success:
        df: pd.DataFrame
        var_attrs: dict

    if failure:
        None
    """

    logger.info(f"{inspect.currentframe().f_code.co_name}: Starting...")

    try: 
        drop_vars_keywords = [
            "qaqc_process",
            "pr_depth",
            "PREC_flag",
            "rsds_flag",
            "_qc",
            "duration",
            "method",
        ]

        # Select variables that contain the keywords defined above
        drop_vars = [
            i for keyword in drop_vars_keywords for i in df.columns if keyword in i
        ]

        # Remove those variables
        df = df.drop(columns=drop_vars)

        # Remove the attributes of the dropped variables
        for key in drop_vars:
            if key in var_attrs:
                del var_attrs[key]

        logger.info(f"{inspect.currentframe().f_code.co_name}: Completed successfully")

        return df, var_attrs

    except Exception as e:
        logger.error(f"{inspect.currentframe().f_code.co_name}: Failed")
        raise e
