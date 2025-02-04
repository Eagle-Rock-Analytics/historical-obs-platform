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


# -----------------------------------------------------------------------------
## distributional gap (unusual gap) + helper functions
def qaqc_unusual_gaps(df, iqr_thresh=5, plots=True, verbose=False, local=False):
    """
    Runs all parts of the unusual gaps function, with a whole station bypass check first.

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    iqr_thresh : int, optional
        interquartile range year threshold, default set to 5
    plots : bool, optional
        if True, produces figures
    verbose : bool, optional
        if True, provides runtime output to the local terminal
    local : bool, optional
        if True, saves plots to local directory

    Returns
    -------
    If QAQC is successful, returns a dataframe with flagged values (see below for flag meaning)
    If QAQC fails, returns None

    Notes
    ------
    Flag meaning : 21,qaqc_unusual_gaps,Part 1: Monthly median value exceeds set threshold limits around monthly interquartile range for the monthly climatological median value
    Flag meaning : 22,qaqc_unusual_gaps,Part 2: Unusual gap in monthly distribution detected beyond PDF distribution
    """

    vars_for_gaps = [
        "tas",
        "tdps",
        "tdps_derived",
        "ps",
        "psl",
        "ps_altimeter",
        "ps_derived",
        "rsds",
    ]

    vars_for_pr = ['pr_5min', 'pr_15min', 'pr_1hr', 'pr_localmid']
    vars_to_check = [var for var in df.columns if var in vars_for_gaps]
    vars_to_pr = [var for var in df.columns if var in vars_for_pr] # precip vars

    try:
        logger.info(
            "Running {} on {}".format("qaqc_unusual_gaps", vars_to_check),
        )

        # whole station bypass check first
        df, stn_length = qaqc_dist_whole_stn_bypass_check(
            df, vars_to_check, min_num_months=iqr_thresh
        )

        # Calculate the number of years for each variable
        # It uses the month with the most (max) number of years (or should it be the min?)
        # TODO: Discuss with Victoria this threshold
        # df is already flagged, just bybass station?
        nYears = np.array([v.max() for k, v in stn_length.items()])
        if (
            nYears < 5
        ).all():  # If all variables have less than 5 years, bypass whole station
            return df
        else:
            df_part1 = qaqc_dist_gap_part1(
                df, vars_to_check, iqr_thresh, plots, verbose=verbose, local=local
            )
            df_part2 = qaqc_dist_gap_part2(
                df_part1, vars_to_check, plots, verbose=verbose, local=local
            )

    except Exception as e:
        logger.info(
            "qaqc_unusual_gaps failed with Exception: {}".format(e),
        )
        return None
    
    # precip flag
    try: 
        if len(vars_to_pr) != 0:
            logger.info("Running qaqc_unusual_gaps_precip on: {}".format(vars_to_pr))
            for var in vars_to_pr:
                df_part3 = qaqc_unusual_gaps_precip(df, var, threshold=200)
            
            return df_part3
    
    except Exception as e:
        logger.info("qaqc_unusual_gaps_precip failed with Exception: {}".format(e))
        return None



# -----------------------------------------------------------------------------
def qaqc_dist_gap_part1(
    df, vars_to_check, iqr_thresh, plot=True, verbose=False, local=False
):
    """Identifies suspect months and flags all obs within month

    Parameters
    ----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    vars_to_check : list of str
        list of variables to run test on
    iqr_thresh : int
        interquartile range year threshold, default set to 5
    plot : bool, optional
        if True, produces figures
    verbose : bool, optional
        if True, provides runtime output to the local terminal
    local : bool, optional
        if True, saves plots to local directory

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with flagged values (see below for flag meaning)

    Notes
    -----
    Part 1 / monthly check
        - compare anomalies of monthly median values
        - standardize against interquartile range
        - compare stepwise from the middle of the distribution outwards
        - asymmetries are identified and flagged if severe
    """
    network = df["station"].unique()[0].split("_")[0]

    for var in vars_to_check:
        logger.info(
            "Running unusual gaps check on: {}, qaqc_dist_gap_part1".format(var),
        )
        for month in range(1, 13):
            monthly_df = df.loc[df["month"] == month]

            # per variable bypass check (first yellow flag being set)
            monthly_df = qaqc_dist_var_bypass_check(monthly_df, var)  # flag here is 20

            # subset for valid obs, distribution drop yellow flags
            df_valid = grab_valid_obs(
                df, var, kind="drop"
            )  # drops data flagged with 20
            if len(df_valid) == 0 or df_valid[var].isnull().all():
                logger.info(
                    "No valid data present for {} in month {} -- skipping to next month".format(
                        var, month
                    ),
                )
                continue  # variable has no valid data

            # calculate monthly climatological median, and bounds
            mid, low, high = standardized_median_bounds(
                df_valid, var, iqr_thresh=iqr_thresh
            )

            # calculate monthly median per month / per year
            df_month = df_valid.groupby(["year"])[var].aggregate("median")

            # Index to flag finds where df_month is out of the distribution
            index_to_flag = (df_month < low) | (df_month > high)

            # Since grouping, the index of df_month is years
            years_to_flag = df_month[index_to_flag].index

            for year in years_to_flag:
                logger.info(
                    "Median {} value for {}-{} is beyond the {}*IQR limits -- flagging month".format(
                        var, month, int(year), iqr_thresh
                    ),
                )

            # flag all obs in that month
            bad = np.logical_and(df["month"] == month, df["year"].isin(years_to_flag))
            df.loc[bad, var + "_eraqc"] = 21  # see era_qaqc_flag_meanings.csv

            if plot:
                if (
                    21 in df[var + "_eraqc"].values
                ):  # don't plot a figure if nothing is flagged
                    dist_gap_part1_plot(
                        df,
                        month,
                        var,
                        flagval=21,
                        iqr_thresh=iqr_thresh,
                        network=network,
                        local=local,
                    )

    return df


# -----------------------------------------------------------------------------
def qaqc_dist_gap_part2(df, vars_to_check, plot=True, verbose=False, local=False):
    """Identifies individual suspect observations and flags the entire month

    Parameters
    -----------
    df : pd.DataFrame
        station dataset converted to dataframe through QAQC pipeline
    vars_to_check : list of str
        list of variables to run test on
    plot : bool, optional
        if True, produces figures
    verbose : bool, optional
        if True, provides runtime output to the local terminal
    local : bool, optional
        if True, saves plots to local directory

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe with flagged values (see below for flag meaning)

    Notes
    -----
    Part 2 / monthly check
        - compare all obs in a single month, all years
        - histogram created from all obs and gaussian distribution fitted
        - threshold values determined using positions where fitted freq falls below y=0.1
        - rounds outwards to next integer plus one
        - going outwards from center, distribution is scanned for gaps which occur outside threshold
        - obs beyond gap are flagged
    """

    for var in vars_to_check:
        logger.info(
            "Running unusual gaps check on: {}, qaqc_dist_gap_part2".format(var),
        )
        for month in range(1, 13):
            # Sel month data
            monthly_df = df.loc[df["month"] == month]

            # per variable bypass check
            monthly_df = qaqc_dist_var_bypass_check(monthly_df, var)  # flag here is 20
            # subset for valid obs, distribution drop yellow flags
            df_valid = grab_valid_obs(
                monthly_df, var, kind="drop"
            )  # drops data flagged with 20
            if len(df_valid) == 0 or df_valid[var].isnull().all() == True:
                logger.info(
                    "No valid data present for {} in month {} -- skipping to next month".format(
                        var, month
                    ),
                )
                continue  # variable has no valid data

            # from center of distribution, scan for gaps (where bin = 0)
            # when gap is found, and it is at least 2x bin width
            # any bins beyond end of gap + beyond threshold value are flagged

            # Bug due to repeated values giving iqr=0 and failing to calculate bins below
            if iqr_range(df_valid, var) == 0:
                logger.info(
                    "No valid data present for {} in month {} -- skipping to next month".format(
                        var, month
                    ),
                )
                continue

            # standardize against IQR range
            df_month_iqr = standardized_iqr(df_valid, var)

            # determine number of bins
            bins = create_bins(df_month_iqr)

            # pdf
            mu = np.nanmean(df_month_iqr)
            sigma = np.nanstd(df_month_iqr)
            y, left_bnd, right_bnd = pdf_bounds(df_month_iqr, mu, sigma, bins)

            # identify gaps as below y=0.1 from histogram, not pdf
            y_hist, bins = np.histogram(df_month_iqr, bins=bins, density=True)
            # identify climatology and iqr baselines in order to flag
            iqr_baseline = iqr_range(df_valid, var=var)
            clim = median_clim(df_valid, var=var)

            # gaps are only flagged for values beyond left_bnd, right_bnd, as long as gap is 2*bin_width (2*0.25)
            # considering that the # of bins for threshold is (4,7) from y=0.1
            # safe to assume that gap is present if values >0.1 outside of left_bnd, right_bnd

            # bins[1:] takes the right edge of bins, suitable for left_bnd
            bins_beyond_left_bnd = np.argwhere(bins[1:] <= left_bnd)
            if len(bins_beyond_left_bnd) != 0:
                for data in bins_beyond_left_bnd:
                    if y_hist[data] > 0.1:  # bins with data > 0.1 beyond left_bnd
                        # identify values beyond left bnd
                        vals_to_flag = clim + (
                            left_bnd * iqr_baseline
                        )  # left_bnd is negative
                        df.loc[df[var] <= vals_to_flag, var + "_eraqc"] = (
                            22  # see era_qaqc_flag_meanings.csv
                        )

            # bins[0:-1] takes the left edge of bins, suitable for left_bnd
            bins_beyond_right_bnd = np.argwhere(bins[:-1] >= right_bnd)
            if len(bins_beyond_right_bnd) != 0:
                for data in bins_beyond_right_bnd:
                    if y_hist[data] > 0.1:  # bins with data > 0.1 beyond right_bnd
                        # identify values beyond right bnd
                        vals_to_flag = clim + (
                            right_bnd * iqr_baseline
                        )  # upper limit threshold
                        # df.loc[df_valid[var] >= vals_to_flag, var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv
                        df.loc[df[var] >= vals_to_flag, var + "_eraqc"] = (
                            22  # see era_qaqc_flag_meanings.csv
                        )

            if plot:
                if (
                    22 in df[var + "_eraqc"].values
                ):  # don't plot a figure if nothing is flagged
                    dist_gap_part2_plot(
                        df,
                        month,
                        var,
                        network=df["station"].unique()[0].split("_")[0],
                        local=local,
                    )

    return df


# -----------------------------------------------------------------------------
def monthly_med(df):
    """Part 1: Calculates the monthly median.

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe

    Returns
    -------
    pd.DataFrame
        resampled dataframe of monthly medians
    """
    return df.resample("M", on="time").median(numeric_only=True)


# #-----------------------------------------------------------------------------
def iqr_range(df, var):
    """Part 1: Calculates the monthly interquartile range

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    var : str
        variable name

    Returns
    -------
    pd.DataFrame
        interquartile range
    """
    return df[var].quantile([0.25, 0.75]).diff().iloc[-1]


# -----------------------------------------------------------------------------
def standardized_iqr(df, var):
    """Part 2: Standardizes data against the interquartile range

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    var : str
        variable name

    Returns
    -------
    pd.DataFrame
        standardized interquartile range
    """
    return (df[var].values - df[var].median()) / iqr_range(df, var)


# -----------------------------------------------------------------------------
def median_clim(df, var):
    """Part 2: Calculate climatological median for a specific month and variable

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    var : str
        variable name

    Returns
    -------
    clim : pd.DataFrame
        climatological median
    """
    clim = df[var].median(numeric_only=True)
    return clim


# -----------------------------------------------------------------------------
def standardized_anom(df, month, var):
    """
    Part 1: Calculates the monthly anomalies standardized by IQR range

    Parameters
    ----------
    df : pd.DataFrame
        input QAQC dataframe
    month : int
        month
    var : str
        variable name

    Returns
    -------
    arr_std_anom : np.array
        array of monthly standardized anomalies for var
    """

    df_monthly_med = monthly_med(df)
    df_clim_med = median_clim(df)
    arr_anom = (
        df_monthly_med.loc[df_monthly_med["month"] == month][var].values
        - df_clim_med.loc[df_clim_med.index == month][var].values
    )
    arr_std_anom = arr_anom / iqr_range(df, month, var)

    return arr_std_anom

# -----------------------------------------------------------------------------

def check_differences(series, threshold=200):

    # Compute pairwise absolute differences, between values in the pandas series and all other values
    # for all values in the column
    diff_matrix = np.abs(series.values[:, None] - series.values)
    
    # Check for values exceeding threshold
    exceeds_threshold = diff_matrix > threshold 
 
    # Exclude self-comparison    
    np.fill_diagonal(exceeds_threshold, True)    
    
    # Identify Rows with Any Exceeding Differences
    rows_with_exceeding_diff = exceeds_threshold.all(axis=1)

    # row_has_diffs_above_threshold = pd.Series(
    return pd.Series(rows_with_exceeding_diff, name="exceeds_threshold", index=series.index)

# -----------------------------------------------------------------------------
def qaqc_unusual_gaps_precip(df, var, threshold=200, verbose=False):
    """
    Precipitation values that are at least threshold larger than all other precipitation totals for a given station and calendar month.
    This is a modification of a HadISD / GHCN-daily test, in which sub-hourly data is aggregated to daily to identify flagged data,
    and flagged values are applied to all subhourly observations within a flagged day.

    Parameters
    ----------
        df : pd.DataFrame
            QAQC dataframe to run through test
        var : str
            variable name
        threshold : int, optional 
            precipitation total to check, default 200 mm
        verbose : boolean, optional 
            whether to provide output to local env

    Returns
    -------
        new_df : pd.DataFrame
            QAQC dataframe with flagged values (see below for flag meaning)

    Notes
    ------
    1. PRELIMINARY: Thresholds/decisions may change with refinement.
    2. HadISD uses a threshold of 300 mm for global precip, we refined to 200 mm for California-specific but can be modified
    3. compare all precipitation obs in a single month, all years
    4. sums observations to daily timestep, then checks each daily sum to every other sum in that month
    5. flags days on which the sum is 200m more than any other daily observation in that month
    Flag Meaning : 33,qaqc_unusual_gaps_precip,Value flagged as an unusual gap in values in the daily precipitation check
    """
    ### Filter df to precipitation variables and sum daily observations

    logger.info("Running qaqc_unusual_gaps_precip on: {}".format(var))
    new_df = df.copy()
    df_valid = grab_valid_obs(new_df, var)

    # aggregate to daily, subset on time, var, and eraqc var
    df_sub = df_valid[["time", 'year','month', 'day', var, var+"_eraqc"]]
    df_dy = df_sub.resample("1D", on="time").agg({var: "sum",var+"_eraqc": "first", "year": "first", "month": "first", "day": "first"})#.reset_index()

    # returns a flag column with True or False
    output = df_dy.groupby("month")[var].transform(check_differences, threshold=threshold)

    # filter the boolean series with itself and set True entries to flag value "33"
    flagged = output[output]
    flagged_str = flagged.map({True: '33'})

    # backflag all observations in the input dataframe
    new_df[var+'_eraqc'] = new_df['time'].dt.date.map(flagged_str)

    # if plot:
    #             if (
    #                 33 in df[var + "_eraqc"].values
    #             ):  # don't plot a figure if nothing is flagged
    #                 output.astype("int").plot()
    
    return new_df
