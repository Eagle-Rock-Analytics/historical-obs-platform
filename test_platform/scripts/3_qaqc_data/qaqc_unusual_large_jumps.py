"""
qaqc_unusual_large_jumps.py

This is a script where Stage 3: QA/QC function(s) on unusual large jumps / spikes within the monthly distribution with data observations are flagged.
For use within the PIR-19-006 Historical Obsevations Platform.

Functions
---------
- qaqc_unusual_large_jumps: Test for unusual large jumps or spikes, given the statistics of the series. Analysis for each individual month in
    time series to account for seasonal cycles in different regions.
- potential_spike_check: Checks for neccessary conditions for a potential spike to be an actual spike.
- detect_spikes: Detect unusual large jumps or ''spikes'' in the time series for `var`.

Intended Use
------------
Script functions for the spike QA/QC test, as a part of the QA/QC pipeline. 
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from log_config import logger

from qaqc_plot import *
from qaqc_utils import *


def qaqc_unusual_large_jumps(
    df: pd.DataFrame,
    iqr_thresh: int = 6,
    min_datapoints: int = 50,
    plot: bool = True,
) -> pd.DataFrame:
    """
    Test for unusual large jumps or spikes, given the statistics of the series. Analysis for each individual month in
    time series to account for seasonal cycles in different regions.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    iqr_thresh : int
        critical value (iqr_thresh*IQR) for spike detection (default=6)
    min_datapoints : int, optional
        minimum data points in each month to be valid for testing (default=50)
    plot : bool, optional
        if True, produces plot and uploads it to AWS

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)

    Notes
    -----
    Flag meaning : 23,qaqc_unusual_large_jumps,Unusual jump (spike) in variable
    """

    logger.info("Running: qaqc_unusual_large_jumps")
    INDEX = df.index
    df = df.copy()
    df.set_index(df["time"], inplace=True)
    df.drop(columns=["time"], inplace=True)

    # Define test variables and check if they are in the dataframe
    check_vars = [
        "pr",
        "pr_5min",
        "pr_15min",
        "pr_1h",
        "pr_24h",
        "pr_localmid",
        "tas",
        "tdps",
        "tdps_derived",
        "ps",
        "psl",
        "ps_altimeter",
        "ps_derived",
    ]
    variables = [var for var in check_vars if var in df.columns]

    logger.info(
        f"Running qaqc_unusual_large_jumps on {variables}"
    )

    # Loop through test variables
    for var in variables:
        try:
            logger.info(
                f"Running unusual large jumps check on: {var}"
            )
            new_df = grab_valid_obs(df, var)  # subset for valid obs

            # first scans suspect values using entire record
            if new_df[var].isna().all() == True:
                continue  # bypass to next variable if all obs are nans

            # Detect spikes
            new_df = detect_spikes(
                new_df, var=var, iqr_thresh=iqr_thresh, min_datapoints=min_datapoints
            )

            # Retrieve location of spikes
            ind = new_df.index[np.where(new_df[var + "_spikes"])[0]]

            # Flag _eraqc variable
            df.loc[ind, var + "_eraqc"] = 23  # see qaqc_flag_meanings.csv

            bad = df.loc[ind, ["year", "month"]]

            df_plot = df.copy()

            if plot:
                ## Plotting by month/year will reduce the number of plots
                keys = bad.groupby(["year", "month"]).groups.keys()
                logger.info(
                    f"Plotting {len(keys)} year/month cases"
                )
                for k in keys:
                    ind = np.logical_and(
                        df_plot["year"] == k[0], df_plot["month"] == k[1]
                    )
                    unusual_jumps_plot(df_plot.loc[ind, :], var, flagval=23)

        except Exception as e:
            logger.info(
                "qaqc_unusual_large_jumps failed on {} with Exception: {}".format(
                    var, e
                )
            )
            continue

    df["time"] = df.index.values
    df = df.set_index(INDEX)
    return df


def potential_spike_check(
    potential_spike: pd.Series, diff: pd.Series, crit: pd.Series, hours_diff: pd.Series
) -> pd.DataFrame:
    """
    Checks for neccessary conditions for a potential spike to be an actual spike.

    Parameters
    ----------
    potential_spike : pd.Series
        bool pd.Series with True on potential spike location
    diff : pd.Series
        float pd.Series with differences in the test variable
    crit : pd.Series
        float pd.Series with the critical value for the differences in the test variable
    hours_diff : pd.Series
        float pd.Series with the hour differences between data points in the test variable

    Returns
    -------
    spikes : pd.DataFrame
        input df with added `var`_spike column True where data matches the spike conditions

    Notes
    -----
    1. Spikes are considered 1-value spike up to 3-values spike
    2. Difference right before the spike should be lower than half the critical value
    3. Difference at the actual spike must be higher than the critical value
    4. Differences within the multi-value spike must lower than half the critical value
    5. Difference right after (spike exit) the spike should be higher than the critical value and
       of opposite sign of the actual spike
    """

    potential_spike = potential_spike.copy()

    ind = np.where(potential_spike)[0]
    spikes = pd.Series(
        np.zeros_like(potential_spike).astype("bool"), index=potential_spike.index
    )
    dates = pd.Series(potential_spike.index.values)

    for i in ind:
        try:
            # Ignore edges for now
            if i == 1 or i >= len(potential_spike) - 4:
                continue
            # Indices, critical values, and values before and after potential spike
            im1, i0, ip1, ip2, ip3, ip4 = [i - 1, i, i + 1, i + 2, i + 3, i + 4]
            tm1, t0, tp1, tp2, tp3, tp4 = diff.iloc[[im1, i0, ip1, ip2, ip3, ip4]]
            cm1, c0, cp1, cp2, cp3, cp4 = crit.iloc[[im1, i0, ip1, ip2, ip3, ip4]]
            # Three-values spike
            if (
                np.sign(t0) != np.sign(tp2)
                and np.abs(tm1) < 0.5 * cm1
                and np.abs(tp1) < 0.5 * cp1
                and np.abs(tp2) < 0.5 * cp2
                and np.abs(tp3) > cp3
                and np.abs(tp4) < 0.5 * cp4
            ):
                spikes.iloc[[i0, ip1, ip2]] = True
                # i += 3
                # continue

            # Two-values spike
            elif (
                np.sign(t0) != np.sign(tp2)
                and np.abs(tm1) < 0.5 * cm1
                and np.abs(tp1) < 0.5 * cp1
                and np.abs(tp2) > cp2
                and np.abs(tp3) < 0.5 * cp3
            ):
                spikes.iloc[[i0, ip1]] = True
                # i += 2
                # continue

            # One-value spike
            elif (
                np.sign(t0) != np.sign(tp1)
                and np.abs(tm1) < 1.0 * cm1
                and np.abs(tp1) > cp1
                and np.abs(tp2) < 1.0 * cp2
            ):
                spikes.iloc[i0] = True
                # i += 1
                # continue
        except Exception as e:
            logger.info("potential_spike_check failed on {}: {}".format(i, e))
            continue

    return spikes


def detect_spikes(
    df: pd.DataFrame, var: str, iqr_thresh: int = 6, min_datapoints: int = 50
) -> pd.DataFrame:
    """
    Detect  unusual large jumps or ''spikes'' in the time series for `var`.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    var : str
        variable to test
    iqr_thresh : int, optional
        critical value (iqr_thresh*IQR) for spike detection (default=6)
    min_datapoints : int, optional
        minimum data points in each month to be valid for testing (default=50)

    Returns
    -------
    df : pd.DataFrame
        input df with added columns for spike check

    Notes
    -----
    1. Find potential unusual large jumps or ''spikes'' by comparing the differences in `var` to each
    month's critical value (crit = iqr_thresh * IQR)
    2. `potential_spike_check` checks for neccessary conditions for a potential spike to be an actual spike
    """

    # Make a copy of the original dataframe
    df = df.copy()

    # Calculate difference in var values
    df[var + "_difference"] = df[var].diff().fillna(0)

    # Calculate dates
    df["date"] = df.index.values

    # Calculate time difference
    df["time_diff"] = df["date"].diff().fillna(pd.Timedelta(0))

    # Calculate time differece in hours
    df["hours_diff"] = df["time_diff"] / np.timedelta64(1, "h")
    df = df[np.logical_and(df["hours_diff"] > 0, df["hours_diff"] <= 12)]

    # Group by month to avoid strong seasonal cycle
    # grouped = df.groupby([pd.Grouper(freq='M'), df['hours_diff']])
    grouped = df.groupby(pd.Grouper(freq="M"))

    # Count number of data per month
    counts = grouped[var + "_difference"].transform("count")
    df[var + "_counts"] = counts
    # Keep only months with more than 50 values to be statistically valid
    df = df[df[var + "_counts"] > min_datapoints]

    # Define modified IQR
    # kwargs = {'rng':(20, 80),}
    kwargs = {}

    # Calculate iqr
    iqr = grouped[var + "_difference"].transform(stats.iqr, **kwargs)
    df[var + "_iqr"] = iqr

    # Calculate critical value as rounded-up 6 (or defined by argument) times IQR
    df[var + "_critical"] = np.ceil(iqr_thresh * df[var + "_iqr"])

    # Find potential spike values where var diff is higher than the critical value
    df[var + "_potential_spikes"] = (
        np.abs(df[var + "_difference"]) > df[var + "_critical"]
    )

    # Filter real spikes using `potential_spike_check` function
    spikes = potential_spike_check(
        df[var + "_potential_spikes"],
        df[var + "_difference"],
        df[var + "_critical"],
        df["hours_diff"],
    )
    df[var + "_spikes"] = spikes

    return df
