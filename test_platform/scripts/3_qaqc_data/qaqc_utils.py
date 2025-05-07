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
import scipy.stats as stats
import sys

## New logger function
from log_config import logger

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

## Following functions
## 1. QA/QC statistical helper functions (e.g., histogram bins)
## 2. QA/QC runtime functions (e.g., progress bar)


## QA/QC statistical helper functions
# -----------------------------------------------------------------------------
# fns to return histogram bins
def get_bin_size_by_var(var: str) -> float:
    """Get bin size for a given variable

    Parameters
    ----------
    var : str
        String variable name

    Returns
    -------
    float
        bin size per variable
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
        "pr_localmid": 0.1,  # mm
        "rsds": 50,  # W/m2
    }
    return bin_size_by_var[var]


# -----------------------------------------------------------------------------
def create_bins_frequent(
    df: pd.DataFrame, var: str, bin_size: float | None = None
) -> np.array:
    """Create bins from data covering entire data range
    Used in frequent value check and qaqc plot of frequent values

    Parameters
    ----------
    df : pd.DataFrame
        Table of the data
    var : str
        String name of the variable
    bin_size : float, optional
        Size of the bins
        See function get_bin_size_by_var for default value

    Returns
    -------
    bins : np.array
        binsizes per variable, frequent qaqc test
    """

    # Get bin size per variable
    # Default bin size is defined in get_bin_size_by_var function
    try:
        bin_size = get_bin_size_by_var(var)
    except:
        bin_size = get_bin_size_by_var("default")

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
def create_bins(data: pd.DataFrame, bin_size: float = 0.25) -> list:
    """Create bins from data covering entire data range.

    Parameters
    ----------
    data : pd.DataFrame
        QAQC data to check
    bin_size : float, optional
        size of bin width

    Returns
    -------
    bins : list of floats
        histogram bins

    Notes
    -----
    To be replaced by above bin functionality to reduce redundancy.
    """

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data)) + bin_size
    bins = np.arange(b_min, b_max, bin_size)

    return bins


# -----------------------------------------------------------------------------
def pdf_bounds(
    df: pd.DataFrame, mu: float, sigma: float, bins: list
) -> tuple[np.array, float, float]:
    """Calculate pdf distribution, return pdf and threshold bounds.

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    mu : float
        mean
    sigma : float
        standard deviation
    bins : list of floats
        histogram bins

    Returns
    -------
    y : list
        pdf distribution values
    left_bnd : float
        leftmost bin
    right_bnd : float
        rightmost bin
    """
    y = stats.norm.pdf(bins, mu, sigma)

    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1).squeeze()
    if len(pdf_bounds) == 0:
        logger.info(
            "PDF distribution warning: there is a bad value present causing issues with pdf y=0.1 determination, returning furthest bins",
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
def qaqc_dist_whole_stn_bypass_check(
    df: pd.DataFrame, vars_to_check: list, min_num_months: int = 5
) -> tuple[pd.DataFrame, int]:
    """
    Checks the number of valid observation months in order to proceed through monthly distribution checks.
    Identifies whether a station record has too few months and produces a fail pass flag.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    vars_to_check : list of str
        list of variables to run whole station bypass check on
    min_num_months : int, optional
        minimum number of months required to pass check, default is 5

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with flagged values (see below for flag meaning)
    stn_length : int
        length of station obs record to determine if bypass is required

    Notes
    -----
    Flag meaning : 19,Yellow flag,Warning: Whole station has too few monthly observations to proceed through any monthly distribution check. Observations were not assessed for quality
    """
    # Copy df to avoid pandas warning
    new_df = df.copy()

    # This piece will return a dictionary with the var name as key, and values are pd.Series with the
    # month and the number of years of data
    global stn_length

    stn_length = map(
        qaqc_var_length_bypass_check, [new_df] * len(vars_to_check), vars_to_check
    )
    stn_length = {k: v for k, v in zip(vars_to_check, stn_length)}

    nYears = np.array([v.max() for k, v in stn_length.items()])

    vars_to_flag = []
    for var in vars_to_check:
        if stn_length[var].max() < min_num_months:
            vars_to_flag.append(var)
    if len(vars_to_flag) > 0:
        # Avoid pandas warning
        # Future-proof way:
        # The safest way to write this is using direct column access by name:
        # This avoids ambiguity and is the most stable across Pandas versions.
        for col in vars_to_flag:
            new_df[col] = 19  # see era_qaqc_flag_meanings.csv

    return new_df, stn_length


# -----------------------------------------------------------------------------
def qaqc_dist_var_bypass_check(
    df: pd.DataFrame, var: str, min_num_months: int = 5
) -> pd.DataFrame:
    """
    Checks the number of valid observation months in order to proceed through monthly distribution checks.
    Identifies whether a station record has too few months and produces a fail pass flag.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    var : str
        variable to run bypass check on
    min_num_months : int, optional
        minimum number of months required to pass check, default is 5

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with flagged values (see below for flag meaning)

    Notes
    -----
    Flag meaning : 20,Yellow flag,Warning: Variable has too few monthly non-nan observations to proceed through any monthly distribution check. Observations were not assessed for quality
    """
    # Copy df to avoid pandas warning
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
def qaqc_var_length_bypass_check(df: pd.DataFrame, var: str) -> pd.DataFrame:
    """Bypass function based on variable obs record length

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    var : str
        variable name

    Returns
    -------
    df : pd.DataFrame
    """
    return (
        df.loc[:, [var, "month", "year"]]
        .groupby(by=["month"])["year"]
        .unique()
        .apply(len)
    )


# -----------------------------------------------------------------------------
# Red vs. Yellow flagging
def grab_valid_obs(
    df: pd.DataFrame, var: str, var2: str | None = None, kind: str = "keep"
) -> pd.DataFrame:
    """Observations that have been flagged by QA/QC test should not proceed through any
    other QA/QC test.

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    var : str
        variable name
    var2 : str
        variable 2 name
    kind : str, optional
        option to keep yellow-flagged data, or drop red-flagged data

    Returns
    -------
    df_valid : pd.DataFrame
        QAQC dataframe with only valid observations returned

    Notes
    -----
    1. Exception is yellow flag (flags: 19, 20), where the entire obs record is too short
    to evaluate for quality. Yellow flag is based on the qaqc_dist_var_bypass_check and
    qaqc_dist_whole_stn_bypass_check checks (in qaqc_unusual_gaps.py)
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


## QA/QC other helper functions
# -----------------------------------------------------------------------------
def progressbar(it: int, prefix: str = "", size: int = 60, out: str = sys.stdout):
    """Print a progress bar to console

    Parameters
    ----------
    it : int
        iternation of list
    prefix : str
        idk
    size : int, optional
        size (length) of progress bar

    Returns
    -------
    None

    Example
    -------
    >>> for i in progressbar(10): # Progress bar of length 10 is printed after each iteration i
    >>> # Loop does something

    References
    ----------
    [1] https://stackoverflow.com/questions/3160699/python-progress-bar
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
def get_file_paths(network: str) -> tuple[str, str, str, str]:
    """Returns AWS filepaths for historical data platform bucket.

    Parameters
    ----------
    network : str
        name of network to find paths for

    Returns
    -------
    rawdir : str
        Path to the raw data bucket
    cleandir : str
        Path to the cleaned data bucket
    qaqcdir : str
        Path to the QAQC'd data bucket
    mergedir : str
        Path to the final merged data bucket
    """
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir


# -----------------------------------------------------------------------------
def get_filenames_in_s3_folder(bucket: str, folder: str) -> list:
    """Get a list of files in s3 bucket.
    Make sure you follow the naming rules exactly for the two function arguments.
    See example in the function docstrings for more details.

    Parameters
    ---------
    bucket : str
        Simply, the name of the bucket, with no slashes, prefixes, suffixes, etc...
    folder : str
        Folder within the bucket that you want the filenames from
        MAKE SURE folder doesn't have a trailing "/"
        i.e. it should be "[folder]", not "[folder]/"

    Returns
    -------
    files_in_s3 : list of str
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
    [1] https://stackoverflow.com/questions/59225939/get-only-file-names-from-s3-bucket-folder
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
def get_wecc_poly(
    terrpath: str, marpath: str
) -> tuple[gp.polygon, gp.polygon, gp.polygon]:
    """Identifies a bbox of WECC area to filter stations against.

    Parameters
    ----------
    terrpath : str
        Path to the WECC terrestrial boundary shapefile
    marpath : str
        Path to the WECC maritime boundary shapefile

    Returns
    -------
    t : gp.polygon
        spatial object for terrestrial boundary
    m : gp.polygon
        spatial object for maritime boudary
    bbox : gp.polygon
        bounding box for union of t and m
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)  ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds  ## Combine polygons and get bounding box of union.
    return t, m, bbox
