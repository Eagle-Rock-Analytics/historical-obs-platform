"""
This is a script where Stage 3: QA/QC function(s) on unusual gaps / gaps within the monthly distribution with data observations are flagged. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import numpy as np
import pandas as pd
import scipy.stats as stats

# New logger function
from log_config import logger

## Import plotting functions
try:
    from qaqc_utils import *
except Exception as e:
    logger.debug("Error importing qaqc_utils: {}".format(e))

try:
    # from qaqc_plot import standardized_median_bounds, dist_gap_part1_plot, dist_gap_part2_plot
    from qaqc_plot import *
except Exception as e:
    logger.debug("Error importing qaqc_plot: {}".format(e))


# -----------------------------------------------------------------------------
#
def is_precip_accumulated(pr):
    """
    Determines whether a precipitation time series is accumulated by analyzing its autocorrelation.

    Parameters
    ----------
    pr : pandas.Series
        A time series of precipitation values, which may contain missing (NaN) values.

    Returns
    -------
    bool
        `True` if the mean autocorrelation of the filtered series is greater than 0.9, indicating
        that the precipitation data is likely accumulated.
        `False` otherwise.

    Notes
    -----
    - The function filters out non-positive and missing values from `pr` before computing autocorrelation.
    - The Pearson autocorrelation is computed using `pandas.Series.autocorr()`, which measures the
      correlation of the series with a lag of 1.
    - If the mean autocorrelation exceeds 0.9, the function assumes that the data is accumulated precipitation.
    - Missing values (`NaN`) in the autocorrelation calculation are handled using `np.nanmean()`
      to avoid bias in the decision.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> pr = pd.Series([0, 0, 5, 10, 15, 20, np.nan, 25, 30, 35])
    >>> is_precip_accumulated(pr)
    True
    """
    test_pr = pr[(pr > 0) & (~pr.isnull())]
    autocorr = test_pr.autocorr()
    if np.nanmean(autocorr) > 0.9:
        return True
    else:
        return False


# -----------------------------------------------------------------------------
#
def flag_ringing(series, window=3, threshold=None):
    """
    Flags values exhibiting ringing (back-and-forth oscillations).

    Parameters
    ----------
    series : pandas.Series
        The input time series.
    window : int, optional
        The window size for detecting frequent oscillations (default is 3).
    threshold : float, optional
        A custom threshold for detecting large fluctuations (default is automatic detection).

    Returns
    -------
    pandas.Series
        A boolean series where `True` indicates a flagged (bad) ringing value.
    """
    diff_series = series.diff()  # Compute first difference
    sign_changes = np.sign(diff_series).diff().fillna(0).abs()  # Detect sign changes

    # Rolling sum of sign changes (to detect frequent oscillations in a window)
    ringing_flags = sign_changes.rolling(window=window, center=True).sum() > (
        window - 1
    )

    # Optional threshold for large alternating fluctuations
    if threshold is None:
        threshold = (
            diff_series.abs().median() * 2
        )  # Use 2x median difference as default threshold

    large_fluctuations = diff_series.abs() > threshold
    ringing_flags = (
        ringing_flags & large_fluctuations
    )  # Flag only if oscillations are significant

    return ringing_flags


# -----------------------------------------------------------------------------
#
def de_accumulate(original_series, reset_threshold=None, window=3, threshold=None):
    """
    Compute incremental values from an accumulated time series while handling resets and filtering artifacts.

    This function converts an accumulated time series (e.g., precipitation, energy usage)
    into incremental values by computing first-order differences. It detects and removes 
    oscillatory ringing values and handles cases where the accumulation resets to zero or a 
    lower value based on a specified threshold.

    Parameters
    ----------
    original_series : pandas.Series
        The accumulated time series. The index should represent time or sequence.
    reset_threshold : float, optional
        Threshold to detect a reset event. If `None`, a reset is assumed when the 
        accumulated value drops to zero. If provided, resets are detected when 
        the difference is smaller than `-reset_threshold`.
    window : int, optional
        Window size for detecting ringing behavior (default is 3).
    threshold : float, optional
        Threshold for identifying ringing values. If `None`, it is determined automatically.

    Returns
    -------
    deaccumulated_series : pandas.Series
        The incremental time series with resets handled and ringing values removed.
    final_flags : pandas.Series
        A boolean series indicating which values were flagged as ringing or reset events.

    Notes
    -----
    - The function first computes the first-order difference of the input series.
    - It applies `flag_ringing()` to detect and remove oscillatory artifacts.
    - Negative differences (except resets) are set to zero, assuming accumulation should not decrease naturally.
    - If `reset_threshold` is provided, it is used to detect unnatural drops; otherwise, any drop to zero is treated as a reset.
    - The cleaned differences are returned with a flagging series to indicate removed values.

    Examples
    --------
    >>> import pandas as pd
    >>> series = pd.Series([0, 5, 10, 20, 0, 3, 7, 15, 25, 0, 4, 9],
    ...                   index=pd.date_range("2024-01-01", periods=12, freq="D"))
    >>> deaccumulated, flags = de_accumulate(series, reset_threshold=15)
    """
    
    series = original_series.copy()
    diff_series = series.copy().diff()

    # print(len(series))
    series = series.loc[original_series.dropna().index]
    diff_series = diff_series.loc[original_series.dropna().index]

    # First check (to see if there is any ringing/up-down measurements, only take the highest
    flags = flag_ringing(diff_series, window=window, threshold=threshold)
    diff_series2 = diff_series.copy().loc[flags[~flags].index].dropna()

    # Only take difference values higher than zero
    flags2 = diff_series2 < 0

    # Merge together >0 flags and original flags
    flags = pd.concat([flags, flags2], axis=1)

    # Rename columns (optional)
    flags.columns = ["flags1", "flags2"]

    # Apply OR operation (handling NaNs as False)
    flags = flags.fillna(False)["flags1"] | flags.fillna(False)["flags2"]

    # Fix the flags where original series is zero
    flags[series == 0] = False

    # print(len(flags))
    # De-accumulate clean series (without ringing)
    # clean_series = original_series.copy().loc[flags[~flags]]
    clean_series = original_series.copy().loc[flags.dropna().index]
    clean_diff_series = clean_series.diff()

    if reset_threshold is not None:
        resets = clean_diff_series < -reset_threshold  # Detect large drops
    else:
        resets = (clean_series.diff() < 0) & (
            clean_series.shift(-1) == 0
        )  # Drop to zero

    clean_diff_series[resets] = np.nan  # Mark reset points
    clean_diff_series.fillna(
        0, inplace=True
    )  # Replace NaNs with zero (or use interpolation if needed)

    # Re-flag negative differences
    flags2 = clean_diff_series < 0

    # Merge together >0 flags and original flags
    flags = pd.concat([flags, flags2], axis=1)

    # Rename columns (optional)
    flags.columns = ["flags1", "flags2"]

    # Apply OR operation (handling NaNs as False)
    flags = flags.fillna(False)["flags1"] | flags.fillna(False)["flags2"]

    # Return the complete difference series
    diff_series = original_series.copy().diff() * np.nan
    diff_series.loc[flags.dropna().index] = clean_diff_series

    final_flags = (diff_series * 0).astype(bool)
    final_flags.loc[flags.dropna().index] = flags
    return diff_series, final_flags


# -----------------------------------------------------------------------------
#
def qaqc_deaccumulate_precip(
    df, var="pr", reset_threshold=50, threshold=10, window=3, plot=True, local=False
):
    """
    Performs quality control (QAQC) and de-accumulation on a precipitation time series.

    This function checks if a given precipitation variable is accumulated and, if so,
    de-accumulates it while handling data quality issues such as ringing. It flags
    affected values accordingly.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing the precipitation data.
    var : str, optional
        The column name of the precipitation variable to be checked and de-accumulated (default is 'pr').
    reset_threshold : float, optional
        The threshold to detect accumulation resets (default is 50). If precipitation drops
        by more than this value, it is considered a reset.
    threshold : float, optional
        The threshold for detecting ringing values in the time series (default is 10).
    window : int, optional
        The rolling window size for detecting ringing values (default is 3).

    Returns
    -------
    pandas.DataFrame or None
        - If precipitation is accumulated, a DataFrame with de-accumulated precipitation and
          additional QAQC flags is returned.
        - If precipitation is not accumulated, the original DataFrame is returned unchanged.
        - If an exception occurs, `None` is returned.

    Notes
    -----
    - The function first determines whether `var` (e.g., 'pr') represents accumulated precipitation
      using `is_precip_accumulated()`.
    - If accumulation is detected:
        - It applies `de_accumulate()` to compute incremental precipitation values.
        - Flags ringing values in `var + "_eraqc"` with flag 34.
        - Saves the original accumulated precipitation in `accum_var` (e.g., `accum_pr`).
        - Flags `accum_var + "_eraqc"` with flag 35 to indicate de-accumulation.
    - If `var` is **not** accumulated, the function **bypasses processing** and returns `df` unchanged.
    - If an error occurs, it logs the exception and returns `None`.
    Flag meaning : 34, qaqc_deaccumulate_precip, Value flagged has a period of oscillating values (probably sensor malfunction, ringing-like data) in the accumulated precipitation (that would result on incorrect/false de-accumulated values that look correct): cannot determine which value (higher or lower) is the real one (it has been also converted to nan in the deaccumulated series)
    Flag meaning : 35, qaqc_deaccumulate_precip, Original precipitation data, deaccumulation process has been applied

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'pr': [0, 5, 10, 20, 0, 3, 7, 15, 25, 0, 4, 9]
    ... })
    >>> df = qaqc_deaccumulate_precip(df)
    >>> df[['pr', 'accum_pr', 'accum_pr_eraqc']]

    """
    vars_for_deacummulation = ["pr"]

    vars_to_check = [var for var in df.columns if var in vars_for_deacummulation]

    df = df.copy()

    try:
        # if True:
        logger.info(
            "Running {} on {}".format("qaqc_deaccumulate_precip", vars_to_check),
        )

        # First, determine if the series is accumulated or instantaneous
        if is_precip_accumulated(df[var]) and len(vars_to_check) > 0:

            # Calculate de-accumulated precip and flag "ringing" values
            diff_series, flags = de_accumulate(
                df[var], reset_threshold=50, window=3, threshold=10
            )

            df.loc[:, var + "_eraqc"] = np.nan
            df.loc[flags, var + "_eraqc"] = 34  # see era_qaqc_flag_meanings.csv

            # Save original accumulated precip into new variable, and de-accumulated into original pr
            # and save de-accumulated precip into original pr
            tmp_var, tmp_index = df[var].values, df[var].index

            # I wonder if this is neccessary, in my opinion it is
            # We should flag those oscillating/ringing values in the de-accumulated
            diff_series.loc[flags] = np.nan

            # Re-assign de-accumulated to pr
            df[var] = diff_series.values
            df["accum_" + var] = tmp_var

            # Flag the new accum_pr to acknowledge that precip was de-accumulated
            df.loc[:, "accum_" + var + "_eraqc"] = 35  # see era_qaqc_flag_meanings.csv
            logger.info(
                "{} on {} done".format("qaqc_deaccumulate_precip", vars_to_check),
            )

            # --------------------------------------------------------
            if plot:
                precip_deaccumulation_plot(
                    df.loc[:, ["pr", "accum_pr", "time", "station"]], flags, local=local
                )
                logger.info("plot produced for precip de-accumulation"),
            return df

        else:  # If it's not accumulated, bypass and return original df
            logger.info("qaqc_deaccumulate_precip bypassed: Precip is not accumulated")

            return df

    except Exception as e:
        # else:
        logger.info(
            "qaqc_deaccumulate_precip failed with Exception: {}".format(e),
        )
        return None


# -----------------------------------------------------------------------------
