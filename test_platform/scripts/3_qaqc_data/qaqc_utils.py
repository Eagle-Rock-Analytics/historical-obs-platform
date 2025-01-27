"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import numpy as np
import time
import pandas as pd
import requests
import urllib
import datetime
import math
import shapely
import xarray as xr
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import scipy.stats as stats
import sys

# # New logger function
# from log_config import logger

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"


## QA/QC helper functions
# -----------------------------------------------------------------------------
def get_bin_size_by_var(var):
    """Get bin size for a given variable

    Parameters
    ----------
    var: str
        String variable name

    Returns
    -------
    float

    """
    # bin sizes: using 1 degC for tas/tdps, and 1 hPa for ps vars
    # all of our pressure vars are in Pa, convert to 100 Pa bin size
    bin_size_by_var = {
        "default": 1,  # Bin size for all other variables
        "ps": 100,  # Pa
        "ps_altimeter": 100,  # Pa
        "psl": 100,  # Pa
        "ps_derived": 100,  # Pa
        "pr_5min": 0.1,  # mm
        "pr_15min": 0.1,  # mm
        "pr_1h": 0.1,  # mm
        "pr_24h": 0.1,  # mm
        "pr_localmid": 0.1, # mm
        "rsds": 50,  # W/m2
    }
    return bin_size_by_var[var]

# -----------------------------------------------------------------------------
def create_bins_frequent(df, var, bin_size=None):
    """Create bins from data covering entire data range
    Used in frequent value check and qaqc plot of frequent values

    Parameters
    ----------
    df: pd.DataFrame
        Table of the data
    var: str
        String name of the variable
    bin_size: float, optional
        Size of the bins
        See function get_bin_size_by_var for default value

    Returns
    -------
    bins: np.array

    """

    # Get bin size per variable
    # Default bin size is defined in get_bin_size_by_var function
    if bin_size is None:
        bin_size = get_bin_size_by_var("default")
    else:
        bin_size = get_bin_size_by_var(var)

    # Get data
    data = df[var]

    # Compute bins
    b_min = np.floor(
        np.nanmin(data)
    )  # Get the minimum of the data; get closest integer to min (floor)
    b_max = np.ceil(
        np.nanmax(data)
    )  # Get the maximum of the data; get closest integer to max (ceil)
    bins = np.arange(
        b_min, b_max + bin_size, bin_size
    )  # Arange the bins; largest bin should be maximum + bin size.

    return bins


# -----------------------------------------------------------------------------
def create_bins(data, bin_size=0.25):
    """Create bins from data covering entire data range"""

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data)) + bin_size
    bins = np.arange(b_min, b_max, bin_size)

    return bins

# -----------------------------------------------------------------------------
def progressbar(it, prefix="", size=60, out=sys.stdout):
    """
    Print a progress bar to console

    Parameters
    ----------
    it: int
        iternation of list
    size: int, optional
        size (length) of progress bar

    Returns
    -------
    progress bar printed to console

    Example
    -------
    >>> for i in progressbar(10): # Progress bar of length 10 is printed after each iteration i
    >>> # Loop does something

    References
    ----------
    https://stackoverflow.com/questions/3160699/python-progress-bar

    """
    count = len(it)
    start = time.time()  # time estimate start

    def show(j):
        x = int(size * j / count)
        # time estimate calculation and string
        remaining = ((time.time() - start) / j) * (count - j)
        mins, sec = divmod(remaining, 60)  # limited to minutes
        time_str = f"{int(mins):02}:{sec:03.1f}"
        print(
            f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} Est wait {time_str}",
            end="\r",
            file=out,
            flush=True,
        )

    show(0.1)  # avoid div/0
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    print("\n", flush=True, file=out)


# -----------------------------------------------------------------------------
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir


# -----------------------------------------------------------------------------
def get_filenames_in_s3_folder(bucket, folder):
    """Get a list of files in s3 bucket.
    Make sure you follow the naming rules exactly for the two function arguments.
    See example in the function docstrings for more details.

    Parameters
    ---------
    bucket: str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc...
    folder: str
        Folder within the bucket that you want the filenames from
        MAKE SURE folder doesn't have a trailing "/"
        i.e. it should be "[folder]", not "[folder]/"

    Returns
    -------
    files_in_s3: list of str
        List of filenames in the bucket

    Example
    -------
    You want to get all the filenames in a s3 bucket with the following path:
    s3 URI: "s3://wecc-historical-wx/1_raw_wx/VALLEYWATER/"
    >>> get_filenames_in_s3_folder(
    >>>    bucket = "wecc-historical-wx",
    >>>    folder = "1_raw_wx/VALLEYWATER"
    >>> )
    ['ValleyWater_6001_1900-01-01_2024-11-11.csv','ValleyWater_6004_1900-01-01_2024-11-11.csv']

    References
    ----------
    https://stackoverflow.com/questions/59225939/get-only-file-names-from-s3-bucket-folder

    """

    s3 = boto3.resource("s3")
    s3_bucket = s3.Bucket(bucket)

    # Get all the filenames
    # Just get relative path (f.key.split(folder + "/")[1])
    files_in_s3 = [
        f.key.split(folder + "/")[1]
        for f in s3_bucket.objects.filter(Prefix=folder).all()
    ]

    # Delete empty filenames
    # I think the "empty" filename/s is just the bucket path, which isn't a file but is read as an object by the objects.filter function
    files_in_s3 = [f for f in files_in_s3 if f != ""]

    return files_in_s3


# -----------------------------------------------------------------------------
def get_wecc_poly(terrpath, marpath):
    """
    Identifies a bbox of WECC area to filter stations against
    Input vars: shapefiles for maritime and terrestrial WECC boundaries
    Returns: spatial objects for each shapefile, and bounding box for their union.
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)  ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds  ## Combine polygons and get bounding box of union.
    return t, m, bbox


# # -----------------------------------------------------------------------------
# # Log print auxiliary functions
# def printf(*args, verbose=True, log_file=None, **kwargs):
#     import datetime

#     tLog = lambda: datetime.datetime.utcnow().strftime("%m-%d-%Y %H:%M:%S") + " : \t"
#     args = [str(a) for a in args]

#     if verbose:
#         if log_file is not None:
#             print(" ".join([tLog(), *args]), **kwargs) or print(
#                 " ".join([tLog(), *args]), file=log_file, **kwargs
#             )
#         else:
#             print(" ".join([tLog(), *args]), **kwargs)
#     else:
#         if log_file is not None:
#             print(" ".join([tLog(), *args]), file=log_file, **kwargs)
#         else:
#             pass


# -----------------------------------------------------------------------------
def create_bins(data, bin_size=0.25):
    """Create bins from data covering entire data range"""

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data)) + bin_size
    bins = np.arange(b_min, b_max, bin_size)

    return bins


# -----------------------------------------------------------------------------
def pdf_bounds(df, mu, sigma, bins):
    """Calculate pdf distribution, return pdf and threshold bounds"""
    y = stats.norm.pdf(bins, mu, sigma)

    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1).squeeze()
    if len(pdf_bounds) == 0:
        logger.info(
            "PDF distribution warning: there is a bad value present causing issues with pdf y=0.1 determination",
        )
        print(
            "PDF distribution warning: there is a bad value present causing issues with pdf y=0.1 determination",
            flush=True,
        )
        return (
            y,
            int(0),
            int(len(y) - 1),
        )  # returning furthest edge cases, return to in V2
        # return (y, bnds[0], bnds[-1]) # returning furthest edge cases, return to in V2

    else:
        # find first index
        # left_bnd = round(bins[pdf_bounds[0] - 1])
        # right_bnd = round(bins[pdf_bounds[-1] + 1])
        # rounds +1 and -1 is giving out of bounds error, using ceil, floor, and clip instead
        try:
            bnds = np.clip(
                [np.floor(pdf_bounds[0]), np.ceil(pdf_bounds[-1])], 0, len(bins) - 1
            ).astype("int")
            bnds = bins[bnds]
            return (y, bnds[0] - 1, bnds[-1] + 1)
        except:
            left_bnd = round(bins[pdf_bounds[0] - 1])
            right_bnd = round(bins[pdf_bounds[-1] + 1])
            return (y, left_bnd, right_bnd)


# -----------------------------------------------------------------------------
def qaqc_dist_whole_stn_bypass_check(df, vars_to_check, min_num_months=5):
    """
    Checks the number of valid observation months in order to proceed through monthly distribution checks.
    Identifies whether a station record has too few months and produces a fail pass flag.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        vars_to_check [list]: list of variables to run whole station bypass check on
        min_num_months [int]: minimum number of months required to pass check, default is 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        pass_flag [str]: pass flag indicating whether whole station pass/fail number of minimum months

    Flag meaning:
    -------------
        19,Yellow flag,Warning: Whole station has too few monthly observations to proceed through any monthly distribution check. Observations were not assessed for quality
    """
    # This piece will return a dictionary with the var name as key, and values are pd.Series with the
    # month and the number of years of data
    global stn_length

    stn_length = map(
        qaqc_var_length_bypass_check, [df] * len(vars_to_check), vars_to_check
    )
    stn_length = {k: v for k, v in zip(vars_to_check, stn_length)}

    nYears = np.array([v.max() for k, v in stn_length.items()])

    for var in vars_to_check:
        if stn_length[var].max() < min_num_months:
            df.loc[:, [var + "_eraqc"]] = 19  # see era_qaqc_flag_meanings.csv

    return df, stn_length


# -----------------------------------------------------------------------------
def qaqc_dist_var_bypass_check(df, var, min_num_months=5):
    """
    Checks the number of valid observation months in order to proceed through monthly distribution checks.
    Identifies whether a station record has too few months and produces a fail pass flag.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        var [string]: variable to run bypass check on
        min_num_months [int]: minimum number of months required to pass check, default is 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        pass_flag [str]: pass flag indicating whether whole station pass/fail number of minimum months

    Flag meaning:
    -------------
        20,Yellow flag,Warning: Variable has too few monthly non-nan observations to proceed through any monthly distribution check. Observations were not assessed for quality
    """

    df = df.copy()

    # if all values are null for that month across years
    if df[var].isnull().all() == True:
        df[var + "_eraqc"] = 20  # see era_qaqc_flag_meanings.csv

    # if more than min_num_months have invalid obs
    if (
        df.groupby(by=["year", "month"])[var]
        .aggregate("median")
        .transform(np.isnan)
        .sum()
        > min_num_months
    ):
        df[var + "_eraqc"] = 20  # see era_qaqc_flag_meanings.csv

    return df


# -----------------------------------------------------------------------------
def qaqc_var_length_bypass_check(df, var):
    return (
        df.loc[:, [var, "month", "year"]]
        .groupby(by=["month"])["year"]
        .unique()
        .apply(len)
    )


# -----------------------------------------------------------------------------
# Red vs. Yellow flagging
def grab_valid_obs(df, var, var2=None, kind="keep"):
    """
    Observations that have been flagged by QA/QC test should not proceed through any
    other QA/QC test.

    Exception is yellow flag (flags: 19, 20), where the entire obs record is too short
    to evaluate for quality. Yellow flag is based on the qaqc_dist_var_bypass_check and
    qaqc_dist_whole_stn_bypass_check checks (in qaqc_unusual_gaps.py)

    kind [str]: options are "keep" to include 19 and 20 flags (yellow setting)
                       and "drop" to exclude 19 and 20 flags (red setting)
    """

    # grab obs with no flags
    df_noflags = df.loc[df[var + "_eraqc"].isnull() == True]

    # retains yellow flagged obs for QA/QC checks where distribution not assessed
    if kind == "keep":
        df_yellowflag = df.loc[
            (df[var + "_eraqc"] == 19) | (df[var + "_eraqc"] == 20)
        ]  # grab obs with yellow flags
        df_merge = pd.concat([df_noflags, df_yellowflag])  # merge together
        df_valid = df_merge.sort_values(by="time")  # sort by time

    # yellow flag set to red flag for distibution checks (gaps + climatological outliers)
    # these obs do not pass through specific QA/QC functions
    elif kind == "drop":  #
        df_valid = df_noflags

    # only applies to some logic checks
    if var2 != None:
        df_valid = df.loc[
            (df[var + "_eraqc"].isnull() == True)
            & (df[var2 + "_eraqc"].isnull() == True)
        ]

    return df_valid
