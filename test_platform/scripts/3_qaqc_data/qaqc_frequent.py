"""
This is a script where Stage 3: QA/QC function(s) on unusually frequent values in the data observations are flagged. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import numpy as np
import pandas as pd
import datetime
import math

# New logger function
from log_config import logger

## Import plotting functions
try:
    from qaqc_plot import *
except:
    logger.debug("Error importing qaqc_plot.py")

try:
    from qaqc_utils import *
except Exception as e:
    logger.debug("Error importing qaqc_utils: {}".format(e))


## frequent values + helper functions
# -----------------------------------------------------------------------------
def qaqc_frequent_vals(df, rad_scheme, plots=True, verbose=False, local=False):
    """
    Test for unusually frequent values, run on temperatures and pressure. This check is performed in two phases.
    - Phase 1: Check is applied to all observations for a designated variable. If the current bin has >50% + >30 number of observations
    compared to +/- 3 surrounding bins, the current bin is highlighted for further check on the year-by-year basis. If the bin persists
    as unusually frequent, the bin is flagged.
    - Phase 2: Check is applied on a seasonal basis, for all observations within that season (mirroring phase 1). If a suspect bin is noted
    in the all observations stage, the check is performed on the year-by-year basis for that season.
    - This test is synergistically applied for air temperature and dew point temperature.

    For precipitation, the folloinwg test is performed:
    - Checks for clusters of 5-9 identical moderate to heavy daily totals in time series of non-zero precipitation observations.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    rad_scheme : str
        radiation handling for frequent occurence of valid zeros
    plots : bool, optional
        if True, produces plots of any flagged data and saved to AWS
    verbose : bool, optional
        if True, provides runtime output to local terminal
    local : bool, optional
        if True, saves plots and log files to local directory

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)
    If QAQC fails, returns None

    Notes
    -----
    Flag meaning : 24,qaqc_frequent_vals,Value flagged as unusually frequent in occurrence at the annual scale after assessing the entire observation record. Temperature and dew point temperature are synergistically flagged.
    Flag meaning : 25,qaqc_frequent_vals,Value flagged as unusually frequent in occurrence at the seasonal scale after assessing the entire observation record. Temperature and dew point temperature are synergistically flagged.

    References
    ----------
    [1] GHCN data description, "Global Historical Climatology Network daily (GHCNd)", URL: https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily
    """

    logger.info("Running: qaqc_frequent_vals")

    # list of var substrings to remove if present in var
    vars_to_remove = ["qc", "eraqc", "duration", "method", "flag", "depth", "accum"]

    non_pr_vars_to_run = [
        "tas",
        "tdps",
        "ps",
        "psl",
        "ps_altimeter",
        "ps_derived",
        "rsds",
    ]
    vars_to_check = [
        var
        for var in df.columns
        if any(True for item in non_pr_vars_to_run if item in var)
        and not any(True for item in vars_to_remove if item in var)
    ]

    pr_vars_to_run = [
        "pr_5min",
        "pr_15min",
        "pr_1h",
        "pr_24h",
        "pr_localmid",
        "pr",
    ]
    pr_vars_to_check = [
        var
        for var in df.columns
        if any(True for item in pr_vars_to_run if item in var)
        and not any(True for item in vars_to_remove if item in var)
    ]

    try:
        for var in vars_to_check:
            logger.info(
                "Running frequent values check on: {}".format(var),
            )
            df_valid = grab_valid_obs(df, var)  # subset for valid obs

            # first scans suspect values using entire record -- this is the per variable check?
            # all years
            if df_valid[var].isna().all() == True:
                continue  # bypass to next variable if all obs are nans

            df_valid = frequent_bincheck(
                df_valid, var, data_group="all", rad_scheme=rad_scheme, verbose=verbose
            )

            # if no values are flagged as suspect, end function, no need to proceed
            if len(df_valid.loc[df_valid[var + "_eraqc"] == 100]) == 0:
                logger.info(
                    "No unusually frequent values detected for entire {} observation record".format(
                        var
                    ),
                )
                # goes to seasonal check, no bypass

            else:
                # year by year
                # then scans for each value on a year-by-year basis to flag if they are a problem within that year
                # DECISION: the annual check uses the unfiltered data
                # previously flagged values are included here -- this would interfere with our entire workflow
                df_valid = frequent_bincheck(
                    df_valid,
                    var,
                    data_group="annual",
                    rad_scheme=rad_scheme,
                    verbose=verbose,
                )

            # seasonal scan (JF+D, MAM, JJA, SON)
            # each season is scanned over entire record to identify problem values
            # only flags applied on annual basis using the three months on their own
            # NOTE: HadISD approach is to use the current year's december, rather than the preceeding december

            # seasonal version because seasonal shift in distribution of temps/dewpoints can reveal hidden values
            # all years
            df_valid = frequent_bincheck(
                df_valid,
                var,
                data_group="seasonal_all",
                rad_scheme=rad_scheme,
                verbose=verbose,
            )  ## DECISION: December is from the current year
            if len(df_valid.loc[df_valid[var + "_eraqc"] == 100]) == 0:
                logger.info(
                    "No unusually frequent values detected for seasonal {} observation record".format(
                        var
                    ),
                )
                continue  # bypasses to next variable

            else:
                logger.info(
                    "Unusually frequent values detected in seasonal distribution, continuing to annual check",
                )
                # year by year --> December selection must be specific
                df_valid = frequent_bincheck(
                    df_valid,
                    var,
                    data_group="seasonal_annual",
                    rad_scheme=rad_scheme,
                    verbose=verbose,
                )

            # remove any lingering preliminary flags, data passed check
            df_valid.loc[df_valid[var + "_eraqc"] == 100, var + "_eraqc"] = np.nan

            # apply unique df_valid flags into full df
            isFlagged = df_valid.loc[df_valid[var + "_eraqc"].isnull() == False]
            for i in isFlagged.index:
                flag_to_place = isFlagged.loc[isFlagged.index == i][
                    var + "_eraqc"
                ].values[0]
                df.loc[isFlagged.index, var + "_eraqc"] = flag_to_place

        # synergistic flag on tas and tdps/tdps_derived
        # first establish at least tas and one tdps var present
        temp_vars = ["tas", "tdps", "tdps_derived"]
        num_temp_vars = [var for var in vars_to_check if var in temp_vars]
        if len(num_temp_vars) != 1 and "tas" in num_temp_vars:
            # proceed to synergistic check
            df = synergistic_flag(df, num_temp_vars)

    except Exception as e:
        logger.info(
            "qaqc_frequent_vals failed with Exception: {}".format(e),
        )
        return None

    try:
        # precip focused check
        for v in pr_vars_to_check:
            df = qaqc_frequent_precip(df, v)

    except Exception as e:
        logger.info(
            "qaqc_frequent_precip failed with Exception: {}".format(e),
        )
        return None

    # plots item
    if plots:
        for var in vars_to_check:
            if (
                24 in df[var + "_eraqc"].unique() or 25 in df[var + "_eraqc"].unique()
            ):  # only plot a figure if a value is flagged
                frequent_vals_plot(df, var, rad_scheme, local=local)

        for v in pr_vars_to_check:
            if 31 in df[v + "_eraqc"].unique():
                frequent_precip_plot(df, v, flag=31, local=local)

    return df


# -----------------------------------------------------------------------------
def frequent_bincheck(df, var, data_group, rad_scheme, verbose=False):
    """Identifies which bins should be flagged via the annual/seasonal frequent test.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    var : str
        variable to run check on
    data_group : str
        annual vs. seasonal handling, options: all, annual, seasonal_all, seasonal_annual
    rad_scheme : str
        radiation handling for frequent occurence of valid zeros
    verbose : bool, optional
        if True, provides runtime output to local terminal

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with flagged values (see below for flag meaning)

    Notes
    -----
    1. histograms created with 0.5 or 1.0 or hpa increments (depending on accuracy of instrument)
    2. each bin compared to the three on either side
    3. if this bin contains more than half the total population of the seven bins combined
    4. and more than 30 observations over the station record (20 for seasonal)
    5. then histogram bin is highlighted for further investigation
    6. minimum number limit imposted to avoid removing true tails of distribution
    """

    # seasons
    szns = [[3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 1, 2]]

    # radiation schemes for assessment
    if var == "rsds":
        if rad_scheme == "all_hours":
            # all valid observations included -- frequent flag will set on 0/nighttime hours
            logger.info(
                "Radiation frequent value check scheme: all_hours selected, will flag nighttime",
            )
            df_to_test = df

        elif rad_scheme == "day_hours":
            # only day hours -- 7am-8pm as "day"
            logger.info(
                "Radiation frequent value check scheme: day_hours selected, day set to 7am - 8pm",
            )
            # 6am PST ~ 1400 UTC, 8pm PST ~ 0400 UTC
            df_to_test = df.loc[(df["hour"] >= 14) | (df["hour"] <= 4)]

        elif rad_scheme == "remove_zeros":
            # remove all zeros -- may remove too many zeros, impact daytime cloudy conditions, regional (PNW)
            logger.info(
                "Radiation frequent value check scheme: remove_zeros selected, may remove valid daytime (cloudy) conditions",
            )
            df_to_test = df.loc[df[var] >= get_bin_size_by_var(var)]

    else:  # all other variables
        df_to_test = df

    # If df_to_test is empty, just skip the next part
    if len(df_to_test) == 0:
        return df

    # all data/annual checks
    if data_group == "all":
        bins = create_bins_frequent(df_to_test, var)
        bar_counts, bins = np.histogram(df_to_test[var], bins=bins)
        flagged_bins = bins_to_flag(bar_counts, bins)

        # flag values in that bin as suspect
        if len(flagged_bins) != 0:
            for sus_bin in flagged_bins:
                # indicate as suspect bins
                # DECISION: preliminary flag? and then remove if okay/reset to nan?
                df.loc[
                    (df[var] >= sus_bin) & (df[var] <= sus_bin + 1), var + "_eraqc"
                ] = 100  # highlight for further review flag, either overwritten with real flag or removed in next step

    # ============================================================================================================
    elif data_group == "annual":
        for yr in df_to_test.year.unique():
            df_yr = df_to_test.loc[df_to_test["year"] == yr]
            if df_yr[var].isna().all() == True:  # some vars will have nan years
                continue
            bins = create_bins_frequent(df_yr, var)  # using 1 degC/hPa bin width
            bar_counts, bins = np.histogram(df_yr[var], bins=bins)
            flagged_bins = bins_to_flag(
                bar_counts, bins, bin_main_thresh=20, secondary_bin_main_thresh=10
            )

            if len(flagged_bins) != 0:
                logger.info(
                    "Flagging bin: {0}".format(flagged_bins),
                )

                for sus_bin in flagged_bins:
                    df.loc[
                        (df["year"] == yr)
                        & (df[var] >= sus_bin)
                        & (df[var] <= sus_bin + 1),
                        var + "_eraqc",
                    ] = 24  # see era_qaqc_flag_meanings.csv

    # ============================================================================================================
    # seasonal checks require special handling
    elif data_group == "seasonal_all":
        for szn in szns:
            df_szn = df_to_test.loc[
                (df_to_test["month"] == szn[0])
                | (df_to_test["month"] == szn[1])
                | (df_to_test["month"] == szn[2])
            ]
            if df_szn[var].isna().all() == True:
                continue
            bins = create_bins_frequent(df_szn, var)  # using 1 degC/hPa bin width
            bar_counts, bins = np.histogram(df_szn[var], bins=bins)
            flagged_bins = bins_to_flag(
                bar_counts, bins, bin_main_thresh=20, secondary_bin_main_thresh=20
            )

            if len(flagged_bins) != 0:
                for sus_bin in flagged_bins:
                    df.loc[
                        (
                            (df["month"] == szn[0])
                            | (df["month"] == szn[1])
                            | (df["month"] == szn[2])
                        )
                        & (df[var] >= sus_bin)
                        & (df[var] <= sus_bin + 1),
                        var + "_eraqc",
                    ] = 100  # highlight for further review flag, either overwritten with real flag or removed in next step

    # ============================================================================================================
    elif data_group == "seasonal_annual":
        for yr in df_to_test.year.unique():
            for szn in szns:
                # all seasons except winter
                if szn != [12, 1, 2]:
                    df_szn = df_to_test.loc[
                        (df_to_test["year"] == yr)
                        & (
                            (df_to_test["month"] == szn[0])
                            | (df_to_test["month"] == szn[1])
                            | (df_to_test["month"] == szn[2])
                        )
                    ]

                    if (
                        df_szn[var].isna().all() == True
                    ):  # some vars will have nan years
                        continue

                    if yr == df_szn.loc[df_szn.index[-1], "year"]:
                        if len(df_szn) == 0:
                            break  # after last season in last year

                    bins = create_bins_frequent(
                        df_szn, var
                    )  # using 1 degC/hPa bin width
                    bar_counts, bins = np.histogram(df_szn[var], bins=bins)
                    flagged_bins = bins_to_flag(
                        bar_counts,
                        bins,
                        bin_main_thresh=15,
                        secondary_bin_main_thresh=10,
                    )

                    if len(flagged_bins) != 0:
                        logger.info(
                            "Flagging bins: {0}".format(flagged_bins),
                        )

                        for sus_bin in flagged_bins:
                            df.loc[
                                (df["year"] == yr)
                                & (
                                    (df["month"] == szn[0])
                                    | (df["month"] == szn[1])
                                    | (df["month"] == szn[2])
                                )
                                & (df[var] >= sus_bin)
                                & (df[var] <= sus_bin + 1),
                                var + "_eraqc",
                            ] = 25  # see era_qaqc_flag_meanings.csv

                # special handling for winter because of december
                else:
                    df_yr = df_to_test.loc[
                        df_to_test["year"] == yr
                    ]  # that year's jan, feb, and wrong dec
                    df_jf = df_yr.loc[
                        (df_yr["month"] == 1) | (df_yr["month"] == 2)
                    ]  # that specific year's jan and feb

                    df_d = df_to_test.loc[
                        (df_to_test["year"] == yr - 1) & (df_to_test["month"] == 12)
                    ]  # previous year's dec
                    if len(df_d) == 0:  # catching very first year instance
                        df_djf = df_jf
                        logger.info(
                            "Winter season: proceeding with just Jan/Feb, no previous Dec",
                        )  ## DECISION

                    else:
                        logger.info(
                            "Winter season: concatenating previous Dec",
                        )
                        df_djf = pd.concat([df_d, df_jf])

                    if (
                        df_djf[var].isna().all() == True
                    ):  # some vars will have nan years
                        continue

                    bins = create_bins_frequent(
                        df_djf, var
                    )  # using 1 degC/hPa bin width
                    bar_counts, bins = np.histogram(df_djf[var], bins=bins)
                    flagged_bins = bins_to_flag(
                        bar_counts,
                        bins,
                        bin_main_thresh=15,
                        secondary_bin_main_thresh=10,
                    )

                    if len(flagged_bins) != 0:
                        logger.info(
                            "Flagging frequent bins in: {0}".format(flagged_bins),
                        )

                        for sus_bin in flagged_bins:
                            # flag jan feb
                            df.loc[
                                (df["year"] == yr)
                                & ((df["month"] == szn[1]) | (df["month"] == szn[2]))
                                & ((df[var] >= sus_bin) & (df[var] <= sus_bin + 1)),
                                var + "_eraqc",
                            ] = 25  # see era_qaqc_flag_meanings.csv
                            # flag correct dec
                            df.loc[
                                ((df["year"] == yr - 1) & (df["month"] == szn[0]))
                                & ((df[var] >= sus_bin) & (df[var] <= sus_bin + 1)),
                                var + "_eraqc",
                            ] = 25  # see era_qaqc_flag_meanings.csv

    return df


# -----------------------------------------------------------------------------
def synergistic_flag(df, num_temp_vars):
    """
    In frequent values, if air temp is flagged, dew point is also flagged, and vice versa.
    Applies appropriate flag in corresponding vars

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    num_temp_vars : list
        list of temperature vars

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with flagged values (see below for flag meaning)
    """

    # need to identify which flag is placed
    # 24 for all obs/years check
    # 25 for all seasons/years check
    flags_to_set = [24, 25]

    # synergistically flag -- if tas is flagged, tdps needs to be flagged, and vice versa
    for flag_to_set in flags_to_set:
        if "tas" in num_temp_vars and "tdps" in num_temp_vars:
            df.loc[df["tas_eraqc"] == flag_to_set, "tdps_eraqc"] = flag_to_set
            df.loc[df["tdps_eraqc"] == flag_to_set, "tas_eraqc"] = flag_to_set

        if "tas" in num_temp_vars and "tdps_derived" in num_temp_vars:
            df.loc[df["tas_eraqc"] == flag_to_set, "tdps_derived_eraqc"] = flag_to_set
            df.loc[df["tdps_derived_eraqc"] == flag_to_set, "tas_eraqc"] = flag_to_set

    return df


# -----------------------------------------------------------------------------
def bins_to_flag(bar_counts, bins, bin_main_thresh=30, secondary_bin_main_thresh=30):
    """Returns the specific bins to flag as suspect.

    Parameters
    ----------
    bar_counts : list
        obs frequency per bin
    bins : list
        bin edges as determined by create_bins_frequent
    bin_main_thresh : int
        min num. of obs for all obs check to proceed to flagging
    secondary_bin_main_thresh : int
        min num of obs for annual/seasonal check to proceed to flagging

    Returns
    --------
    bins_to_flag : list of float
        list of bins to flag for frequent values check
    """

    bins_to_flag = []  # list of bins that will be flagged

    for i in range(0, len(bar_counts)):
        # identify main bin + 3 on either side
        bin_end = i + 4

        # need handling for first 3 blocks as there is no front
        if i < 3:
            bin_start = 0
        else:
            bin_start = i - 3

        bin_block_sum = bar_counts[
            bin_start:bin_end
        ].sum()  # num of obs in the 7-bin block
        bin_main_sum = bar_counts[i]  # num of obs in main bin

        # determine whether main bin is more than half sum in 7-block bin
        bin_block_50 = bin_block_sum * 0.5  # primary check at 50%
        bin_block_90 = bin_block_sum * 0.9  # secondary check at 90%

        if (bin_main_sum > bin_block_50) == True:
            # ensure that bin_main_sum is greater than bin_main_thresh
            if bin_main_sum > bin_main_thresh:
                bins_to_flag.append(math.floor(bins[i]))

                # annual/seasonal check
                if (bin_main_sum > bin_block_90) == True:
                    if bin_main_sum > secondary_bin_main_thresh:
                        bins_to_flag.append(math.floor(bins[i]))

            else:  # less than bin_main_thresh obs in bin_main_sum, do not indicate as suspect
                continue

    return bins_to_flag  # returns a list of values that are suspicious


# -----------------------------------------------------------------------------
# precipitation focused precip check
def qaqc_frequent_precip(df, var, moderate_thresh=18, day_thresh=5, verbose=False):
    """Checks for clusters of 5-9 identical moderate to heavy daily totals in time series of non-zero precipitation observations.
    This is a modification of a HadISD / GHCN-daily test, in which sub-hourly data is aggregated to daily to identify flagged data,
    and flagged values are applied to all subhourly observations within a flagged day.

    Parameters
    ----------
    df : pd.DataFrame
        QAQC dataframe to run through test
    var : str
        variable name
    moderate_thresh : int, optional
        moderate precipitation total to check, default 18mm (~0.7 inch for Santa Clara County, may be different for other regions)
    day_thresh : int, optional
        num. of min consecutive days to flag, default 5 days
    verbose : bool, optional
        whether to provide output to local env

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with test applied

    Notes
    -----
    Flag meaning : 31,qaqc_frequent_precip,Value flagged as unusually frequent in occurrence from daily precipitation above moderate rain total

    References
    ----------
    [1] GHCN data description, "Global Historical Climatology Network daily (GHCNd)", URL: https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily
    """

    logger.info("Running qaqc_frequent_precip on: {}".format(var))

    new_df = df.copy()
    df_valid = grab_valid_obs(new_df, var)  # subset for valid obs

    # add check in case valid_obs is now length 0
    if len(df_valid) == 0:
        logger.info("{} has 0 observations, moving to next variable.".format(var))
        return new_df

    # aggregate to daily, subset on time, var, and eraqc var
    df_sub = df_valid[["time", var, var + "_eraqc"]]
    df_dy = df_sub.resample("1D", on="time").sum().reset_index()

    # identify non-zero precip totals
    df_nozero = df_dy.loc[df_dy[var] > 0]
    df_nozero = df_nozero.copy()  # is this really necesssary?

    # creates new column to identify consecutive values
    df_nozero["consecutive"] = (
        df_nozero[var]
        .groupby((df_nozero[var] != df_nozero[var].shift()).cumsum())
        .transform("size")
    )

    # filter to get rows day_thresh min num. of consecutive days and obs above moderate threhsold
    flagged_days = df_nozero[
        (df_nozero["consecutive"] >= day_thresh) & (df_nozero[var] > moderate_thresh)
    ]

    # flag all values within flagged days
    if len(flagged_days) != 0:
        new_df.loc[
            (
                (
                    new_df["year"].isin(flagged_days.time.dt.year)
                    & new_df["month"].isin(flagged_days.time.dt.month)
                    & new_df["day"].isin(flagged_days.time.dt.day)
                )
            ),
            var + "_eraqc",
        ] = 31  # see flag meanings
        logger.info(
            "Flagging {} days for frequent value precip check for {}".format(
                len(flagged_days), var
            )
        )

    return new_df
