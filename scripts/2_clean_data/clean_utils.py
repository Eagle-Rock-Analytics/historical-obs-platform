"""
clean_utils.py

Functions
---------
- var_to_unique_list: Given a variable column in an xarray object, this function returns a string with all of the unique variables
    present in that column.
- get_file_paths: Given a network name, return all relevant AWS filepaths for other functions.

Intended Use
------------
Support utility functions for cleaning processes, as a part of the clean pipeline.
"""

import pandas as pd
import numpy as np
import xarray as xr


def var_to_unique_list(ds: xr.Dataset, column: str) -> str:
    """
    Given a variable column in an xarray object, this function returns a string with all of the unique variables
    present in that column. This is used to generate lists of qaqc flag values from existing data in the cleaning stage.

    Parameters
    ----------
    ds : xr.Dataset
        the xarray Dataset containing the variable/column to process
    column : str
        column is the name of the column

    Returns
    -------
    str : flagvals
        string that can be provided as the flag_values attribute for a QA/QC flag.
    """
    flagvals = ds[column].values.tolist()[0]
    flagvals = [x for x in flagvals if pd.isnull(x) == False]  # Remove nas
    flagvals = list(np.unique(flagvals))  # Get unique values
    flagvals = " ".join(flagvals)
    return flagvals


#
def get_file_paths(network: str) -> tuple[str, str, str]:
    """
    Given a network name, return all relevant AWS filepaths for other functions.

    Parameters
    ----------
    network : str
        name of the weather station network

    Returns
    -------
    tuple[str, str, str]
        tuple list of rawdir, cleandir, qaqcdir
    """
    rawdir = f"1_raw_wx/{network}/"
    cleandir = f"2_clean_wx/{network}/"
    qaqcdir = f"3_qaqc_wx/{network}/"
    return rawdir, cleandir, qaqcdir
