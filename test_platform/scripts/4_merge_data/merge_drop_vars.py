"""
This is a script where unwanted variables are dropped, for step 6 of the merge pipeline
"""

## Import Libraries
import pandas as pd

# New logger function
from merge_log_config import logger


# ----------------------------------------------------------------------
def delete_vars(df: pd.DataFrame, var_attrs: dict) -> tuple[pd.DataFrame, dict]:
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
    drop_vars_keywords = [
        "qaqc_process",
        "pr_duration",
        "pr_depth",
        "PREC_flag",
        "rsds_duration",
        "rsds_flag",
        "_qc",
        "method",
    ]

    # Select variables that contain the keywords defined above
    drop_vars = [
        i for keyword in drop_vars_keywords for i in df.columns if keyword in i
    ]

    # Remove those variables
    df = df.drop(columns=drop_vars)

    # Remove the attributes of the dropped variables
    for key in drop_vars_keywords:
        if key in var_attrs:
            del var_attrs[key]

    return df, var_attrs
