"""
This is a script where Stage 3: QA/QC function(s) on climatological outlier values in the data observations are flagged. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import boto3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.signal as signal

# New logger function
from log_config import logger

## Import plotting functions
try:
    from qaqc_plot import *
except:
    logger.debug("Error importing qaqc_plot.py")

try:
    from qaqc_unusual_gaps import *
except:
    logger.debug("Error importing qaqc_unusual_gaps.py")

try:
    from qaqc_utils import *
except Exception as e:
    logger.debug("Error importing qaqc_utils: {}".format(e))


# ----------------------------------------------------------------------
## climatological outlier check
def qaqc_climatological_outlier(
    df,
    winsorize=True,
    winz_limits=[0.05, 0.05],
    bin_size=0.25,
    plot=True,
    verbose=False,
    local=False,
):
    """Flags individual gross outliers from climatological distribution

    Parameters
    -----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    winsorize : bool
        if True, raw observations are winsorized to remove spurious outliers first
    winz_limits : list of floats
        if winsorize is True, values represent the low and high percentiles to standardize to
    bin_size : float
        size of distribution bins
    plot : bool, optional
        if True, produces plots of any flagged data and saved to AWS
    verbose : bool, optional
        if True, provides runtime output to local terminal
    local : bool, optional
        if True, retains local copy of figures

    Returns
    -------
    qaqc success:
        new_df : pd.DataFrame
            QAQC dataframe with flagged values (see below for flag meaning)
    qaqc failure:
        None
            This function does not return a value

    Notes
    ------
    Flag meaning : 26,qaqc_climatological_outlier,Value flagged as a climatological outlier

    References
    ----------
    [1] GHCN data description, "Global Historical Climatology Network daily (GHCNd)", URL: https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily

    """
    new_df = df.copy()

    vars_to_check = ["tas", "tdps", "tdps_derived"]
    pr_vars = ["pr_5min", "pr_15min", "pr_1h", "pr_24h", "pr_localmid"]
    vars_to_anom = [v for v in vars_to_check if v in df.columns]
    pr_vars_to_anom = [v for v in pr_vars if v in df.columns]

    try:
        logger.info(
            "Running {} on {}".format(
                "qaqc_climatological_outlier", vars_to_anom + pr_vars_to_anom
            )
        )

        # whole station bypass check first
        # df is already flagged by gaps function (careful if the order is modified)

        for var in vars_to_anom:
            # only work with non-flagged values
            logger.info("Checking for climatological outliers in: {}".format(var))
            df_valid = grab_valid_obs(
                new_df, var, kind="drop"
            )  # subset for valid obs, distribution drop yellow flags

            df_valid = df_valid.dropna(subset=var)  # Keep only useful columns
            df_valid = df_valid[[var, "year", "month", "day", "hour", "time"]]

            # Bypass if there are no valid observations remaining
            if df_valid[var].size == 0:
                logger.info(
                    "No valid (unflagged) observations for: {} to proceed through qaqc_climatological_outlier. Bypassing station.".format(
                        var
                    )
                )
                continue

            # Winsorize data and calculate climatology by month/hour with winsorized data
            if winsorize:
                clim = df_valid.groupby(["month", "hour"])[var].transform(
                    lambda row: stats.mstats.winsorize(row, limits=winz_limits)
                )
            clim = pd.DataFrame(
                data={var: clim, "hour": df_valid.hour, "month": df_valid.month},
                index=df_valid.index,
            )
            clim = clim.groupby(["month", "hour"])[var].transform(
                lambda row: np.nanmean(row)
            )
            clim = pd.DataFrame(
                data={var: clim, "hour": df_valid.hour, "month": df_valid.month},
                index=df_valid.index,
            )

            # Anomalize using winsorized month/hour climatologies
            anom = df_valid[var] - clim[var]
            anom = pd.DataFrame(
                data={var: anom, "hour": df_valid.hour, "month": df_valid.month},
                index=df_valid.index,
            )

            # Calculate IQR by month/hour
            iqr = anom.groupby(["month", "hour"])[var].transform(
                lambda row: max(
                    np.nanpercentile(row, 75) - np.nanpercentile(row, 25), 1.5
                )
            )
            iqr = pd.DataFrame(
                data={var: iqr, "hour": df_valid.hour, "month": df_valid.month},
                index=df_valid.index,
            )

            # Standardise anomalies by month/hour IQR
            std = anom[var] / iqr[var]
            std = pd.DataFrame(
                data={var: std, "hour": df_valid.hour, "month": df_valid.month},
                index=df_valid.index,
            )

            # Low-pass standardised data
            cut_freq = 1 / (3600 * 24 * 365 / 30)  # In Hz (cut_period : 1 month)
            data_freq = 1 / (
                df_valid["time"].diff().mode().values[0].astype("float") / 1e9
            )  # In Hz
            sos = signal.butter(1, cut_freq, "lp", output="sos", fs=data_freq)
            filtered = signal.sosfilt(sos, std[var].interpolate(method="linear"))
            df_valid["raw_" + var] = df_valid[var]
            df_valid[var] = filtered

            # Flag outliers
            df_valid["flag"] = df_valid.groupby(["month", "hour"])[var].transform(
                lambda row: flag_clim_outliers(row, bin_size=bin_size)
            )

            # Drop all non-flagged values
            df_valid = df_valid.dropna(subset=["flag"])
            if len(df_valid) != 0:
                logger.info(
                    "Flagging outliers in climatological outlier check for {0}".format(
                        var
                    )
                )  # only print statement if flags are set

            # Flag original data
            new_df.loc[new_df.time.isin(df_valid.time), var + "_eraqc"] = df_valid[
                "flag"
            ]

    except Exception as e:
        logger.info(
            "qaqc_climatological_outlier failed with Exception: {}".format(e),
        )
        return None

    try:
        # precip focused check
        for var in pr_vars_to_anom:
            new_df = qaqc_climatological_outlier_precip(new_df, var)

    except Exception as e:
        logger.info(
            "qaqc_climatological_outlier_precip failed with Exception: {}".format(e)
        )
        return None

    # Plot flagged values
    if plot:

        ## CURRENT PLOT IS FAILING -- NEED TO RE-EVALUATE
        ## CAN USE FLAGGED_TIMESERIES_PLOT FOR THE TIME BEING

        # # Save original df for plotting
        # df_plot = new_df.copy()

        # station = df["station"].unique()[0]
        # for var in vars_to_anom:
        #     print(var)
        #     if 26 in new_df[var + "_eraqc"].unique():  # only plot if flag is present
        #         # Extract only flagged values to loop over those months and hours
        #         df_plot = df_plot[
        #             ["year", "hour", "month", "time", var+"_eraqc", var]
        #         ].set_index(["month", "hour"])
        #         print(df_plot.columns)

        #         # Loop over flagged months/hours
        #         index = df_valid.set_index(["month", "hour"]).index.unique()
        #         print(index)
        #     for i, ind in enumerate(index):
        #         # Extract actual month/hour from index
        #         month, hour = ind

        #         # Plot distribution
        #         clim_outlier_plot(
        #             df_plot.loc[ind][var],
        #             month,
        #             hour,
        #             bin_size=bin_size,
        #             station=station,
        #             local=local,
        #         )
        for var in pr_vars_to_anom:
            if 32 in new_df[var + "_eraqc"].unique():  # only plot if flag is present
                climatological_precip_plot(new_df, var, flag=32)

    return new_df


# ----------------------------------------------------------------------
def flag_clim_outliers(series, bin_size=0.25):
    """Identifies climatological outliers to flag.

    Parameters
    ----------
    series : pd.DataFrame
        QAQC dataframe
    bin_size : float
        bin size for distribution

    Returns
    -------
    clim_outliers : pd.DataFrame
        QAQC dataframe with flags [?]
    """
    # If series is small (less than 5 years) skip to next month/hour
    if len(series) <= 5:
        return np.ones_like(series) * np.nan

    # Calculate frequency, normal fit, and boumdaries for clim outliers
    # freq, bins, p, left, right = fit_normal(series, bin_size=0.10, plot=True)
    freq, bins, p, left, right = fit_normal(series, bin_size=bin_size, plot=False)

    # Calculate bins for frequency checks
    freq_bins = np.concatenate(
        (bins[1 : int(len(bins) / 2)], [0, 0], bins[int(len(bins) / 2) + 1 : -1])
    )

    # Flag series given the distribution thresholds (left and right)
    flag = gap_search(freq, left, right)

    # -------------------------------------------------------------------
    # Red left side of the distribution
    left_bad_bins = freq_bins[np.logical_and(flag == -1, freq_bins < 0)]
    if len(left_bad_bins) > 0:
        red_left = series <= left_bad_bins.max()
    else:
        red_left = np.zeros_like(series).astype("bool")

    # Red right side of the distribution
    right_bad_bins = freq_bins[np.logical_and(flag == -1, freq_bins > 0)]
    if len(right_bad_bins) > 0:
        red_right = series >= right_bad_bins.max()
    else:
        red_right = np.zeros_like(series).astype("bool")

    # Red flags
    red = np.logical_or(red_left, red_right)

    # ----------------------------------------------------------------------
    # Yellow left side of the distribution
    left_probable_bins = freq_bins[np.logical_and(flag == 0, freq_bins < 0)]
    if len(left_probable_bins) > 0:
        yellow_left = np.logical_and(series <= left_probable_bins.max(), ~red_left)
    else:
        yellow_left = np.zeros_like(series).astype("bool")

    # Yellow right side of the distribution
    right_probable_bins = freq_bins[np.logical_and(flag == 0, freq_bins > 0)]
    if len(right_probable_bins) > 0:
        yellow_right = np.logical_and(series >= right_probable_bins.min(), ~red_right)
    else:
        yellow_right = np.zeros_like(series).astype("bool")

    # Yellow flags
    yellow = np.logical_or(yellow_left, yellow_right)

    # Create new array for the flags
    clim_outliers = np.ones_like(series) * np.nan

    #############################################################################
    # RED vs YELLOW flag with nearest neighbor stations for version 2
    # clim_outliers[red] = 30
    # clim_outliers[yellow] = 31

    # RED AND YELLOW FLAGGED
    clim_outliers[np.where(np.logical_or(red, yellow))[0]] = 26

    # ONLY RED FLAGGED
    # clim_outliers[np.where(red)[0]] = 26
    #############################################################################

    return clim_outliers


# ----------------------------------------------------------------------
def fit_normal(series, bin_size=0.25, plot=False):
    """Fits a guassian distribution to the series.

    Parameters
    ----------
    series : pd.DataFrame
        QAQC dataframe
    bin_size : float
        bin size for distribution
    plot : bool, optional
        whether to plot the data

    Returns
    -------
    freq : list of ?
        frequency [?]
    bins : list of ints/floats ?
        bins for distribution
    p : list of floats
        pdf of distribution
    left : int
        leftmost bin for distribution
    right : int
        rightmost bin for distribution
    """
    bins = create_bins(series, bin_size=bin_size)
    max_bin = np.abs(bins).max()
    bins = np.arange(-max_bin - bin_size, max_bin + 2 * bin_size, bin_size)

    freq, bins = np.histogram(
        series,
        bins=bins,
    )
    area = sum(np.diff(bins) * freq)

    # Fit a normal distribution to the data
    mu, std = stats.norm.fit(series)
    p = stats.norm.pdf(bins, mu, std) * area

    try:
        left = np.where(np.logical_and(np.gradient(p) > 0, p <= 0.1))[0][
            -1
        ]  # +1 # Manually shift the edge by one bin
    except:
        left = 1
    try:
        right = np.where(np.logical_and(np.gradient(p) < 0, p <= 0.1))[0][
            0
        ]  # -1 # Manually shift the edge by one bin
    except:
        right = len(bins) - 2

    if plot:  ## why is this here?
        # Plot the histogram of the series
        fig, ax = plt.subplots()
        ax.hist(series, bins=bins, density=False, alpha=0.35, label="Histogram")
        ax.hist(series, bins=bins, density=False, lw=1.5, color="C0", histtype="step")
        ax.plot(bins, p, "k", linewidth=2, label="Gaussian Fit")

        ax.set_yscale("log")
        ymin = min(0.08, freq[freq > 0].min())
        ymax = np.ceil(freq.max() / 100) * 100
        ax.set_ylim(ymin, ymax)
        ax.set_title("Histogram with Gaussian Fit")
        ax.set_xlabel("Value")
        ax.set_ylabel("Frequency")
        ax.legend()

        ax.axvline(bins[left], c="k", ls="--")
        ax.axvline(bins[right], c="k", ls="--")
        ax.axhline(0.1, c="k", ls=":")

    return freq, bins, p, left, right


# ----------------------------------------------------------------------
def gap_search(freq, left, right):
    """DOCUMENTATION NEEDED.

    Inputs
    ------
    freq : list of ?
        frequency [?]
    left : int
        leftmost bin for distribution
    right : int
        rightmost bin for distribution

    Returns
    -------
    flag : int [?]
        [?]

    """
    left_freq = freq[0:left]
    left_flag = np.zeros_like(
        left_freq
    )  # Yellow flag, all values beyond the threshold are flagged
    for i, f in zip(range(len(left_freq) - 1, -1, -1), left_freq[::-1]):
        if f < 0.1:
            left_flag[0 : i + 1] = (
                -1
            )  # Red flag, values and gap below 0.1 and beyond the threshold are flagged
            break

    right_freq = freq[right + 1 :]
    right_flag = np.zeros_like(
        right_freq
    )  # Yellow flag, all values beyond the threshold are flagged
    for i, f in zip(range(len(right_freq)), right_freq):
        if f < 0.1:
            right_flag[i:] = (
                -1
            )  # Red flag, values and gap below 0.1 and beyond the threshold are flagged
            break

    # Return flag of the size of freq
    flag = np.ones(len(freq) - len(left_freq) - len(right_freq))
    flag = np.concatenate((left_flag, flag, right_flag))
    return flag


# ----------------------------------------------------------------------
def qaqc_climatological_outlier_precip(df, var, factor=9):
    """Checks for daily precipitation totals that exceed the respective 29-day climatological 95th percentiles by at
    least a certain factor (9 when the day's mean temperature is above freezing, 5 when it is below freezing).
    This is a modification of a HadISD / GHCN-daily test, in which sub-daily data is aggregated to daily to identify flagged data,
    and flagged values are applied to all sub-daily observations within a flagged day.

    Parameters
    ----------
    df : pd.DataFrame
        QAQC dataframe
    var : str
        variable name
    factor : int, optional
        multiplication factor for severity of climatological exceedance, default 9

    Returns
    -------
    new_df : pd.DataFrame
        QAQC dataframe

    Notes
    -----
    Flag meaning : 32,qaqc_climatological_outlier_precip,Value flagged as a climatological outlier in daily precipitation check

    References
    ----------
    [1] GHCN data description, "Global Historical Climatology Network daily (GHCNd)", URL: https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily

    To Do
    ------
    1. Incorporate temperature check if temperature is present for station (V2)
    """

    logger.info("Checking for climatological outliers in: {}".format(var))

    new_df = df.copy()
    df_valid = grab_valid_obs(new_df, var)  # subset for valid obs

    # aggregate to daily
    df_sub = df_valid[["time", "month", "year", "day", var, var + "_eraqc"]]
    df_dy = (
        df_sub.resample("1D", on="time")
        .agg(
            {
                var: "sum",
                var + "_eraqc": "first",
                "month": "first",
                "year": "first",
                "day": "first",
            }
        )
        .reset_index()
    )

    # calculate the respective 29-day 95th percentile
    # v1: using the month each day is located in for "29-day"
    for mon in range(1, 13):
        df_mon = df_dy.loc[df_dy.month == mon]

        # subset for days with >0mm rain
        df_mon = df_mon.loc[df_mon[var] > 0]

        # calculate percentile
        p95 = df_mon[var].quantile(0.95)

        # identify where factor x percentile is exceeded and flag
        if p95 != 0:

            # handling for months with low number of non-zero precip days
            # some unusual spikes not being caught by other tests
            if p95 > 442.0:
                # largest recorded 1-day rainfall was 17.6 inches (442 mm) Feb 17 1986
                # https://cepsym.org/Sympro1994/Goodridge.pdf
                flagged_days = df.mon.loc[df_mon[var] > 442.0]
                new_df.loc[
                    (
                        new_df.year.isin(flagged_days.time.dt.year)
                        & new_df.month.isin(flagged_days.time.dt.month)
                        & new_df.day.isin(flagged_days.time.dt.day)
                    ),
                    var + "_eraqc",
                ] = 32
                if len(flagged_days) != 0:
                    logger.info(
                        "Flagging {} days in month {} for climatological outlier precip check for {}".format(
                            len(flagged_days), mon, var
                        )
                    )

            else:
                flagged_days = df_mon.loc[df_mon[var] > factor * p95]
                new_df.loc[
                    (
                        new_df.year.isin(flagged_days.time.dt.year)
                        & new_df.month.isin(flagged_days.time.dt.month)
                        & new_df.day.isin(flagged_days.time.dt.day)
                    ),
                    var + "_eraqc",
                ] = 32
                if len(flagged_days) != 0:
                    logger.info(
                        "Flagging {} days in month {} for climatological outlier precip check for {}".format(
                            len(flagged_days), mon, var
                        )
                    )

        elif p95 == 0:
            # p95 is zero in this case
            flagged_days = df_mon.loc[df_mon[var] > factor]
            new_df.loc[
                (
                    new_df.year.isin(flagged_days.time.dt.year)
                    & new_df.month.isin(flagged_days.time.dt.month)
                    & new_df.day.isin(flagged_days.time.dt.day)
                ),
                var + "_eraqc",
            ] = 32
            if len(flagged_days) != 0:
                logger.info(
                    "ZERO -- Flagging {} days in month {} for climatological outlier precip check for {}".format(
                        len(flagged_days), mon, var
                    )
                )

    return new_df
