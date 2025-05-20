"""
This is a script where columns are reordered into a standard order.
"""

## Import Libraries
import pandas as pd

# New logger function
from merge_log_config import logger


# -----------------------------------------------------------------------------
def merge_reorder_vars(df: pd.DataFrame) -> pd.DataFrame:
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
    ##### Reorder variables
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

    return df
