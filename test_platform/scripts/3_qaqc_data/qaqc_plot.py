"""
This is a script where Stage 3: QA/QC related common plotting functions stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from io import BytesIO
import scipy.stats as stats

# New logger function
from log_config import logger
import time

# =======================================================================================================
# IGNORE PERFORMANCE WARNING FOR NOW
# Related to:
# PerformanceWarning: indexing past lexsort depth may impact performance.
# clim_outlier_plot(df_plot.loc[i][var], month, hour, bin_size=0.1, station=station, local=True)
# It's because we are using multiindex and is not sorted, for now leaving like that since
# it's not a big deal. Check for V2 if ordering de multiindex would preserve order for final/original df
import warnings

warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
# =======================================================================================================
try:
    from qaqc_utils import create_bins_frequent, create_bins
except Exception as e:
    logger.debug("Error importing qaqc_utils: {}".format(e))

from IPython.display import display


# ============================================================================================================
# All plots helper plotting function for labeling, units, min, maxes
def _plot_format_helper(var):
    """Helper function for unusual large jumps plots.

    Parameters
    ----------
    var : str
        variable name being plotted

    Returns
    -------
    ylab : str
        variable name for y-axis
    unit : str
        units of variable
    miny : float
        min var value for y axis
    maxy : float
        max var value for y axis
    """

    pr_vars = ["pr", "pr_5min", "pr_15min", "pr_1h", "pr_24h", "pr_localmid"]
    ps_vars = ["ps", "psl", "ps_derived", "ps_altimeter"]

    if var == "tas":
        ylab = "Air Temperature at 2m"
        unit = "K"

    elif var == "tdps" or var == "tdps_derived":
        ylab = "Dewpoint Temperature"
        unit = "K"

    elif var == "sfcWind":
        ylab = "Surface Wind Speed"
        unit = "$m s^{-1}$"

    elif var == "sfcWind_dir":
        ylab = "Surface Wind Direction"
        unit = "degrees"

    elif var == "rsds":
        ylab = "Surface Radiation"
        unit = "$W m^{-2}$"

    elif var == "hurs":
        ylab = "Humidity"
        unit = "%"

    elif var in pr_vars:
        ylab = "Precipitation"
        unit = "mm"

    elif var in ps_vars:
        ylab = "Pressure"
        unit = "Pa"

    elif var == "elevation":
        ylab = "Elevation"
        unit = "m"

    T_X = {"North_America": 329.92}  # K
    T_N = {"North_America": 210.15}  # K
    D_X = {"North_America": 329.85}  # K
    D_N = {"North_America": 173.15}  # K
    W_X = {"North_America": 113.2}  # m/s
    W_N = {"North_America": 0.0}  # m/s
    S_X = {"North_America": 108330}  # Pa
    S_N = {"North_America": 87000}  # Pa
    R_X = {"North_America": 1500}  # W/m2
    R_N = {"North_America": -5}  # W/m2

    # for other non-record variables (wind direction, precipitation)
    N_X = {"North_America": 360}  # degrees
    N_N = {"North_America": 0}  # degrees
    P_X = {"North_America": 1000}  # mm, arbitrarily set
    P_N = {"North_America": 0}  # mm
    H_X = {"North_America": 100}  # humidity max
    H_N = {"North_America": 0}  # humidity min
    E_X = {"North_America": 6210.0}  # m
    E_N = {"North_America": -100}  # m

    maxes = {
        "tas": T_X,
        "tdps": D_X,
        "tdps_derived": D_X,
        "sfcWind": W_X,
        "sfcWind_dir": N_X,
        "ps": S_X,
        "psl": S_X,
        "ps_altimeter": S_X,
        "ps_derived": S_X,
        "rsds": R_X,
        "pr": P_X,
        "pr_5min": P_X,
        "pr_15min": P_X,
        "pr_1h": P_X,
        "pr_24h": P_X,
        "pr_localmid": P_X,
        "hurs": H_X,
        "elevation": E_X,
    }
    mins = {
        "tas": T_N,
        "tdps": D_N,
        "tdps_derived": D_N,
        "sfcWind": W_N,
        "sfcWind_dir": N_N,
        "ps": S_N,
        "psl": S_N,
        "ps_altimeter": S_N,
        "ps_derived": S_N,
        "rsds": R_N,
        "pr": P_N,
        "pr_5min": P_N,
        "pr_15min": P_X,
        "pr_1h": P_N,
        "pr_24h": P_N,
        "pr_localmid": P_N,
        "hurs": H_N,
        "elevation": E_N,
    }
    miny = mins[var]["North_America"]
    maxy = maxes[var]["North_America"]

    return ylab, unit, miny, maxy


# ============================================================================================================
## flagged timeseries plot helper
def id_flag(flag_to_id):
    """Identifies flag based on numerical value assigned for plotting.

    Parameters
    ----------
    flag_to_id : int
        specific flag to identify

    Returns
    -------
    fn_name : str
        name of QA/QC flag
    """

    flag_df = pd.read_csv("era_qaqc_flag_meanings.csv")
    fn_name = flag_df.loc[flag_df["Flag_value"] == int(flag_to_id)][
        "QAQC_function"
    ].values[0]

    return fn_name


# ============================================================================================================
## flagged timeseries plot
def flagged_timeseries_plot(df, var, dpi=None, local=False, savefig=True):
    """Produces timeseries of variables that have flags placed.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    var : str
        variable name
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        if True, saves plot locally, else: only saves plot to AWS
    savefig : bool, optional
        if True, produces plot and saves either locally or to AWS

    Returns
    -------
    None
        This function does not return a value
    """

    # first check if var has flags, only produce plots of vars with flags
    if len(df[var + "_eraqc"].dropna().unique()) != 0:
        # create figure
        fig, ax = plt.subplots(figsize=(10, 3))

        # plot all observations
        df.plot(
            ax=ax,
            x="time",
            y=var,
            marker=".",
            ms=4,
            lw=1,
            color="k",
            alpha=0.5,
            label="Cleaned data",
        )

        # identify flagged data, can handle multiple flags
        for flag in df[var + "_eraqc"].dropna().unique():
            flag_name = id_flag(flag)
            flag_label = "{:.3f}% of data flagged by {}".format(
                100 * len(df.loc[df[var + "_eraqc"] == flag, var]) / len(df), flag_name
            )

            flagged_data = df[~df[var + "_eraqc"].isna()]
            flagged_data.plot(
                x="time",
                y=var,
                ax=ax,
                marker="o",
                ms=7,
                lw=0,
                mfc="none",
                color="C3",
                label=flag_label,
            )

            legend = ax.legend(loc=0, prop={"size": 8})

        # plot aesthetics
        ylab, units, miny, maxy = _plot_format_helper(var)
        plt.ylabel("{} [{}]".format(ylab, units))
        plt.xlabel("")
        plt.title(
            "Full station timeseries: {0}".format(df["station"].unique()[0]),
            fontsize=10,
        )

        # save to AWS
        if savefig:
            bucket_name = "wecc-historical-wx"
            directory = "3_qaqc_wx"
            img_data = BytesIO()
            plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
            img_data.seek(0)

            s3 = boto3.resource("s3")
            bucket = s3.Bucket(bucket_name)
            network = df["station"].unique()[0].split("_")[0]
            figname = "flagged_timeseries_{0}_{1}".format(
                df["station"].unique()[0], var
            )
            bucket.put_object(
                Body=img_data,
                ContentType="image/png",
                Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
            )

            # save locally if needed
            if local:
                fig.savefig(
                    "qaqc_figs/{}.png".format(figname),
                    format="png",
                    dpi=dpi,
                    bbox_inches="tight",
                )

            # close figure to save memory
            plt.close()

        # Useful completion statement
        logger.info("Flag summary plot produced on: {}".format(var))

        return None


# ============================================================================================================
## frequent values plotting functions
def frequent_plot_helper(df, var, bins, flag, yr, rad_scheme, dpi=None, local=False):
    """Plotting helper with common plotting elements for all 3 versions of this plot.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    var : str
        variable name
    bins : list
        histogram bin limits
    flag : int
        QA/QC flag to subset
    yr : int
        year of data subset, annotation
    rad_scheme : str
        radiation scheme, default is "remove_zeros"
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        if True, saves plot locally, else: only saves plot to AWS

    Returns
    -------
    None
        This function does not return a value
    """

    # plot all valid data within year/season
    _plot = df.plot.hist(column=var, bins=bins, color="k", legend=False, alpha=0.5)

    # plot flagged values
    # first identify which values are flagged
    vals_to_flag = df.loc[df[var + "_eraqc"] == flag][var].unique()

    bars_to_flag = []
    for i in vals_to_flag:
        if math.isnan(i) == False:
            bars_to_flag.append(math.floor(i))

    # flag bars if too frequent
    for bar in _plot.patches:
        x = bar.get_x()
        if x in bars_to_flag:  # right tail
            bar.set_color("r")

    # plot aesthetics
    xlab, units, miny, maxy = _plot_format_helper(var)
    plt.xlabel("{0} [{1}]".format(xlab, units))
    yr_formatted = str(yr).replace("_", " ")  # simple formatting for plot aesthetic
    plt.annotate(yr_formatted, xy=(0.02, 0.95), xycoords="axes fraction", fontsize=10)
    plt.title("Frequent value check: {}".format(df["station"].unique()[0]), fontsize=10)
    plt.legend(("Cleaned data", "Flagged"), loc="upper right")
    ax = plt.gca()
    leg = ax.get_legend()
    leg.legendHandles[0].set_color("k")  # set valid to black
    leg.legendHandles[-1].set_color("r")  # set flagged bar to red

    if var == "rsds":
        plt.annotate(
            "Sfc. radiation option: \n{}".format(rad_scheme),
            xy=(0.02, 0.85),
            xycoords="axes fraction",
            fontsize=10,
        )

    elif var == "tdps" or var == "tdps_derived":
        plt.annotate(
            "Dewpoint temperature\nand air temperature are\nsynergistically flagged.",
            xy=(0.02, 0.85),
            xycoords="axes fraction",
            fontsize=8,
        )

    # save figure to AWS
    network = df["station"].unique()[0].split("_")[0]

    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    img_data = BytesIO()
    plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    figname = "qaqc_frequent_{0}_{1}_{2}".format(df["station"].unique()[0], var, yr)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
    )

    # save locally if needed
    if local:
        fig.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figure to save memory
    plt.close()

    return None


# -----------------------------------------------------------------------------------------
def frequent_vals_plot(df, var, rad_scheme, local=False):
    """
    Produces a histogram of the diagnostic histogram per variable,
    and any bin that is indicated as "too frequent" by the qaqc_frequent_vals test
    is visually flagged.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    var : str
        variable name
    rad_scheme : str
        radiation scheme, default is "remove_zeros"
    local : bool, optional
        whether to produce local plots for review

    Returns
    -------
    None
        This function does not return a value
    """

    bins = create_bins_frequent(df, var)

    # first identify which values are flagged and "where"
    # Year-by-year flag (24): plot all data for that year
    flag_df = df.loc[df[var + "_eraqc"] == 24]
    if len(flag_df) != 0:
        # identify year(s) with flagged data
        plot_yrs = flag_df["year"].unique()

        for y in plot_yrs:
            df_to_plot = df.loc[df["year"] == y]
            _plot = frequent_plot_helper(
                df_to_plot, var, bins, flag=24, yr=y, rad_scheme=rad_scheme
            )

    # Seasonal flag (25): plot all data for that year and season + specific handling for winter
    flag_df = df.loc[df[var + "_eraqc"] == 25]
    if len(flag_df) != 0:
        # identify unique years with flagged seasonal data
        plot_yrs = flag_df["year"].unique()

        for y in plot_yrs:
            df_year = df.loc[df["year"] == y]  # grab the entire year

            flagged_szns = df_year.loc[df_year[var + "_eraqc"] == 25][
                "month"
            ].unique()  # identify flagged months in that year

            if (
                3 in flagged_szns or 4 in flagged_szns or 5 in flagged_szns
            ):  # Spring - MAM
                df_to_plot = df_year.loc[
                    (df_year["month"] == 3)
                    | (df_year["month"] == 4)
                    | (df_year["month"] == 5)
                ]
                _plot = frequent_plot_helper(
                    df_to_plot,
                    var,
                    bins,
                    flag=25,
                    yr=str(y) + "_spring",
                    rad_scheme=rad_scheme,
                )

            if (
                6 in flagged_szns or 7 in flagged_szns or 8 in flagged_szns
            ):  # Summer - JJA
                df_to_plot = df_year.loc[
                    (df_year["month"] == 6)
                    | (df_year["month"] == 7)
                    | (df_year["month"] == 8)
                ]
                _plot = frequent_plot_helper(
                    df_to_plot,
                    var,
                    bins,
                    flag=25,
                    yr=str(y) + "_summer",
                    rad_scheme=rad_scheme,
                )

            if (
                9 in flagged_szns or 10 in flagged_szns or 11 in flagged_szns
            ):  # Autumn - SON
                df_to_plot = df_year.loc[
                    (df_year["month"] == 9)
                    | (df_year["month"] == 10)
                    | (df_year["month"] == 11)
                ]
                _plot = frequent_plot_helper(
                    df_to_plot,
                    var,
                    bins,
                    flag=25,
                    yr=str(y) + "_autumn",
                    rad_scheme=rad_scheme,
                )

            if 12 in flagged_szns:  # Winter - current year D + next year JF
                # special handling as follows
                # if the next year has flagged jan/feb, this will overwrite, but will be identical figure
                # some years will not have current year december and next year jan/feb so need this edge case
                df_d = df_year.loc[df_year["month"] == 12]  # current year dec
                df_jf = df.loc[
                    (df["year"] == y + 1) & ((df["month"] == 1) | (df["month"] == 2))
                ]  # next year jan+feb
                df_to_plot = pd.concat([df_d, df_jf])
                _plot = frequent_plot_helper(
                    df_to_plot,
                    var,
                    bins,
                    flag=25,
                    yr=str(y + 1) + "_winter",
                    rad_scheme=rad_scheme,
                )

            if (
                1 in flagged_szns or 2 in flagged_szns
            ):  # Winter - previous year D + current year JF
                # special handling as follows
                # if the previous year has flagged december, this will overwrite, but will be identical figure
                # some years will not have previous year december and current jan/feb so need this edge case
                df_d = df.loc[
                    (df["year"] == y - 1) & (df["month"] == 12)
                ]  # previous year dec
                df_jf = df_year[
                    (df_year["month"] == 1) | (df_year["month"] == 2)
                ]  # current year jan+feb
                df_to_plot = pd.concat([df_d, df_jf])
                _plot = frequent_plot_helper(
                    df_to_plot,
                    var,
                    bins,
                    flag=25,
                    yr=str(y) + "_winter",
                    rad_scheme=rad_scheme,
                )

    return None


# -----------------------------------------------------------------------------------------
def frequent_precip_plot(df, var, flag, dpi=None, local=False):
    """Plot frequent values for precipitation.

    Parameters
    ----------
    df : pd.DataFrame
        input QA/QC dataframe to produce plot on
    var : str
        variable name, precipitation vars only
    flag : int
        qaqc_precip_check flag (31)
    dpi : int, optional
        resolution of figure
    local : bool, optional
        whether to save figure locally in addition to AWS

    Returns
    -------
    None
        This function does not return a value
    """
    # valid precipitation variables

    fig, ax = plt.subplots()

    # plot all cleaned data
    df.plot(
        ax=ax,
        x="time",
        y=var,
        marker=".",
        ms=4,
        lw=1,
        color="k",
        alpha=0.5,
        label="Cleaned data",
    )

    # plot all flagged data
    flagged_df = df.loc[df[var + "_eraqc"] == flag]
    flagged_df.plot(
        ax=ax,
        x="time",
        y=var,
        marker="o",
        ms=7,
        lw=0,
        mfc="none",
        color="C3",
        label="Flagged data",
    )

    # plot aesthetics
    plt.legend(loc="best")
    ylab, units, miny, maxy = _plot_format_helper(var)
    plt.ylabel("{} [{}]".format(ylab, units))
    plt.xlabel("")
    plt.title(
        "Frequent values -- precipitation: {0}".format(df["station"].unique()[0]),
        fontsize=10,
    )

    # save figure to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    network = df["station"].unique()[0].split("_")[0]
    img_data = BytesIO()
    plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    figname = "qaqc_frequent_value_check_{0}_{1}".format(df["station"].unique()[0], var)
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
    )

    # save locally if needed
    if local:
        plt.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figure to save memory
    plt.close()

    return None


# ============================================================================================================
## distribution gap plotting functions
def dist_gap_part1_plot(
    df, month, var, flagval, iqr_thresh, network, dpi=None, local=False
):
    """Produces a timeseries plots of specific months and variables for part 1 of the unusual gaps function.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    month : int
        month to subset
    var : str
        variable name
    flagval : int
        QA/QC flag to subset
    iqr_thresh : int
        min threshold multiplier for IQR to flag
    network : str
        name of network
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        whether to produce local plots for review

    Returns
    -------
    None
        This function does not return a value
    """

    # grab data by months
    df = df.loc[df["month"] == month]

    # Skip if all values are NaN
    if df[var].isnull().all():
        return

    # grab flagged data
    flag_vals = df.loc[df[var + "_eraqc"] == flagval]

    # plot valid data
    ax = df.plot.scatter(x="time", y=var, color="k", label="Cleaned data")

    # plot flagged data
    flag_name = id_flag(flagval)
    flag_label = "{:.3f}% of data flagged by {}".format(
        100 * len(flag_vals) / len(df), flag_name
    )

    flag_vals.plot(
        x="time",
        y=var,
        ax=ax,
        marker="o",
        ms=7,
        lw=0,
        mfc="none",
        color="C3",
        label=flag_label,
    )

    # plot climatological median and threshold * IQR range
    mid, low_bnd, high_bnd = standardized_median_bounds(df, var, iqr_thresh)

    plt.axhline(y=mid, color="k", lw=0.5, label="Climatological monthly median")
    plt.fill_between(
        x=df["time"].values,
        y1=low_bnd,
        y2=high_bnd,
        alpha=0.25,
        color="0.75",
        label="{} * IQR".format(iqr_thresh),
    )

    # plot aesthetics
    plt.legend(loc="best")
    ylab, units, miny, maxy = _plot_format_helper(var)
    plt.ylabel("{} [{}]".format(ylab, units))
    plt.xlabel("")
    plt.title(
        "Distribution gap check pt 1: {0} / month: {1}".format(
            df["station"].unique()[0], month
        ),
        fontsize=10,
    )

    # save figure to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    img_data = BytesIO()
    plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    figname = "qaqc_dist_gap_check_part1_{0}_{1}_{2}".format(
        df["station"].unique()[0], var, month
    )
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
    )

    # save locally if needed
    if local:
        plt.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figure to save memory
    plt.close()

    return None


# -----------------------------------------------------------------------------------------
def dist_gap_part2_plot(df, month, var, network, dpi=None, local=False):
    """Produces a histogram of the monthly standardized distribution
    with PDF overlay and threshold lines where pdf falls below y=0.1.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    month : int
        month to subset
    var : str
        variable name
    network : str
        name of network
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        whether to produce local plots for review

    Returns
    -------
    None
        This function does not return a value
    """

    # select month
    df = df.loc[df["month"] == month]

    # standardize against IQR range
    df_month_iqr = standardized_iqr(df, var)

    # determine number of bins
    bins = create_bins(df_month_iqr)

    # plot histogram
    ax = plt.hist(
        df_month_iqr, bins=bins, log=False, density=True, color="k", alpha=0.3
    )
    xmin, xmax = plt.xlim()
    plt.ylim(ymin=0.1)

    # plot pdf
    mu = np.nanmean(df_month_iqr)
    sigma = np.nanstd(df_month_iqr)
    y = stats.norm.pdf(bins, mu, sigma)
    l = plt.plot(bins, y, "k--", linewidth=1)

    # add vertical lines to indicate thresholds where pdf y=0.1
    try:
        pdf_bounds = np.argwhere(y > 0.1).squeeze()
        # v2 refine
        if len(pdf_bounds) != 0:
            # find first index
            left_bnd = np.floor(pdf_bounds[0])
            right_bnd = np.ceil(pdf_bounds[-1])
            thresholds = (left_bnd - 1, right_bnd + 1)

            plt.axvline(thresholds[1], color="r")  # right tail
            plt.axvline(thresholds[0], color="r")  # left tail

            # flag (visually) obs that are beyond threshold
            for bar in ax[2].patches:
                x = bar.get_x() + 0.5 * bar.get_width()
                if x > thresholds[1]:  # right tail
                    bar.set_color("r")
                elif x < thresholds[0]:  # left tail
                    bar.set_color("r")
    except:
        logger.info(
            "dist_gap_part2_plot: PDF boundaries issue -- skipping left and right tails"
        )

    # title and useful annotations
    plt.title(
        "Distribution gap check, {0}: {1}".format(df["station"].unique()[0], var),
        fontsize=10,
    )
    plt.annotate(
        "Month: {}".format(month),
        xy=(0.025, 0.95),
        xycoords="axes fraction",
        fontsize=8,
    )
    plt.annotate(
        "Mean: {}".format(round(mu, 3)),
        xy=(0.025, 0.9),
        xycoords="axes fraction",
        fontsize=8,
    )
    plt.annotate(
        "Std.Dev: {}".format(round(sigma, 3)),
        xy=(0.025, 0.85),
        xycoords="axes fraction",
        fontsize=8,
    )
    plt.ylabel("Frequency (obs)")

    # save figure to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    img_data = BytesIO()
    plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    figname = "qaqc_dist_gap_check_part2_{0}_{1}_{2}".format(
        df["station"].unique()[0], var, month
    )
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
    )

    # Save locally if needed
    if local:
        plt.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figure to save memory
    plt.close()

    return None


# ============================================================================================================
## unusual large jumps plotting functions
def unusual_jumps_plot(df, var, flagval=23, dpi=None, local=False):
    """Plots unusual large jumps qaqc.

    Parameters
    ----------
    df : pd.Dataframe
        station pd.DataFrame from qaqc pipeline
    var : str
        variable name
    flagval : int, optional
        flag value to plot (23 for unusual large jumps)
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        if True, saves plot locally, else: only saves plot to AWS

    Returns
    -------
    None
        This function does not return a value
    """

    fig, ax = plt.subplots(figsize=(10, 3))

    # Plot variable and flagged data
    df[var].plot(
        ax=ax, marker=".", ms=4, lw=1, color="k", alpha=0.5, label="Cleaned data"
    )

    flag_label = "{:.4f}% of data flagged".format(
        100 * len(df.loc[df[var + "_eraqc"] == flagval, var]) / len(df)
    )
    df.loc[df[var + "_eraqc"] == flagval, var].plot(
        ax=ax, marker="o", ms=7, lw=0, mfc="none", color="C3", label=flag_label
    )

    # plot other flags
    other_flags = np.logical_and(
        ~df[var + "_eraqc"].isnull() == True, df[var + "_eraqc"] != flagval
    )
    if other_flags.any():
        df.loc[other_flags, var].plot(
            ax=ax, marker="o", ms=7, lw=0, mfc="none", color="C4", label="other flags"
        )

    legend = ax.legend(loc=0, prop={"size": 8})

    station = df["station"].unique()[0]
    network = station.split("_")[0]

    # Plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper(var)
    ylab = "{0} [{1}]".format(ylab, units)
    ax.set_ylabel(ylab)
    ax.set_xlabel("")

    # Set time for fig namne
    month = df.month.unique()[0]
    year = df.year.unique()[0]

    title = "Unusual large jumps check: {0}".format(station)
    ax.set_title(title, fontsize=10)

    # save to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    figname = "qaqc_figs/qaqc_unusual_large_jumps_{0}_{1}_{2}-{3}".format(
        station, var, year, month
    )

    key = "{0}/{1}/{2}.png".format(directory, network, figname)
    img_data = BytesIO()
    fig.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(Body=img_data, ContentType="image/png", Key=key)

    if local:
        fig.savefig(
            "{}.png".format(figname), format="png", dpi=dpi, bbox_inches="tight"
        )

    # close figure to save memory
    plt.close("all")

    return None


# ============================================================================================================
def clim_outlier_plot(
    series, month, hour, bin_size=0.1, station=None, dpi=None, local=False
):
    """Produces a histogram of monthly standardized distribution
    with PDF overlay and threshold lines where pdf falls below y=0.1.
    Differs from dist_gap_part2_plot for the climatological outlier
    as IQR standardization does not occur within plotting.

    Parameters
    ----------
    series : pd.Dataframe
        station pd.DataFrame from qaqc pipeline
    month : int
        month to subset
    hour : int
        hour to subset
    bin_size : int
        width of bins
    station : str
        station name
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        if True, saves plot locally, else: only saves plot to AWS

    Returns
    -------
    None
        This function does not return a value
    """

    var = series._name

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

    index = np.arange(len(bins))
    good = np.where(np.logical_and(index >= left, index <= right))[0]
    good_freq = good[:-1]

    # Plot the histogram of the series
    fig, ax = plt.subplots()
    ax.stairs(
        freq, bins, alpha=1, color="C3", label="Clim outliers".format(month, hour)
    )

    ax.stairs(
        freq[good_freq],
        bins[good],
        alpha=0.35,
        fill=True,
        color="C0",
        label="Distribution",
    )
    ax.stairs(freq[good_freq], bins[good], alpha=1, color="C0")
    ax.set_yscale("log")
    ymin = min(0.05, freq[freq > 0].min())
    ymax = np.ceil(freq.max() / 100) * 100
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="upper right")

    ax.axvline(bins[left], c="k", ls=":", alpha=0.8)
    ax.axvline(bins[right], c="k", ls=":", alpha=0.8)
    ax.axhline(0.1, c="k", ls=":", alpha=0.8)

    # title and useful annotations
    box = dict(facecolor="white", edgecolor="white", alpha=0.85)
    plt.title(
        "Climatological outlier check, {0}: {1}".format(station, var),
        fontsize=10,
    )
    plt.annotate(
        "Month: {}".format(month),
        xy=(0.025, 0.93),
        xycoords="axes fraction",
        fontsize=10,
        bbox=box,
    )
    plt.annotate(
        "Hour: {}".format(hour),
        xy=(0.025, 0.87),
        xycoords="axes fraction",
        fontsize=10,
        bbox=box,
    )
    plt.ylabel("Frequency (obs)")

    # save figure to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    img_data = BytesIO()
    plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)

    network = station.split("_")[0]
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    figname = "qaqc_climatological_outlier_{0}_{1}_{2}_{3}".format(
        station, var, month, hour
    )
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
    )

    if local:
        plt.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figures to save memory
    plt.close(fig)

    return None


# ============================================================================================================
def climatological_precip_plot(df, var, flag, dpi=None, local=False):
    """Plot frequent values for precipitation.

    Parameters
    -----------
    df : pd.DataFrame
        input QA/QC dataframe to produce plot on
    var : str
        variable name, precipitation vars only
    flag : int
        qaqc_precip_check flag (31)
    dpi : int, optional
        resolution of figure
    local : bool, optional
        if True, whether to save figure locally in addition to AWS

    Returns
    -------
    None
        This function does not return a value
    """
    # valid precipitation variables

    fig, ax = plt.subplots(figsize=(10, 3))

    # plot all cleaned data
    df.plot(
        ax=ax,
        x="time",
        y=var,
        marker=".",
        ms=4,
        lw=1,
        color="k",
        alpha=0.5,
        label="Cleaned data",
    )

    # plot all flagged data
    flagged_df = df.loc[df[var + "_eraqc"] == flag]
    flagged_df.plot(
        ax=ax,
        x="time",
        y=var,
        marker="o",
        ms=7,
        lw=0,
        mfc="none",
        color="C3",
        label="Flagged data",
    )

    # plot aesthetics
    plt.legend(loc="best")
    ylab, units, miny, maxy = _plot_format_helper(var)
    plt.ylabel("{} [{}]".format(ylab, units))
    plt.xlabel("")
    plt.title(
        "Climatological outliers -- precipitation: {0}".format(
            df["station"].unique()[0]
        ),
        fontsize=10,
    )

    # save figure to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    network = df["station"].unique()[0].split("_")[0]
    img_data = BytesIO()
    plt.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)

    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    figname = "qaqc_climatological_outlier_{0}_{1}".format(
        df["station"].unique()[0], var
    )
    bucket.put_object(
        Body=img_data,
        ContentType="image/png",
        Key="{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname),
    )

    # save locally if needed
    if local:
        plt.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figure to save memory
    plt.close()

    return None


# ============================================================================================================
def unusual_streaks_plot(
    df, var, flagvals=(27, 28, 29), station=None, dpi=None, local=False
):
    """Plots unusual large jumps qaqc result.

    Parameters
    ----------
    df : pd.Dataframe
        station data from qaqc pipeline
    var : str
        variable name
    flagval : int, optional
        flag value to plot (27, 28, 29 for unusual streaks)
    dpi : int, optional
        resolution for png plots
    local : bool, optional
        if True, saves plot locally, else: only saves plot to AWS

    Returns
    -------
    None
        This function does not return a value
    """

    fig, ax = plt.subplots(figsize=(10, 3))

    # Plot variable and flagged data
    df.plot(
        ax=ax,
        x="time",
        y=var,
        marker=".",
        ms=4,
        lw=1,
        color="k",
        alpha=0.5,
        label="Cleaned data",
    )

    # grab flagged data
    flag_vals_0 = df.copy()[["time", var]]
    flag_vals_0.loc[:, var] = np.nan
    flag_vals_0.loc[df[var + "_eraqc"] == 27, var] = df.loc[
        df[var + "_eraqc"] == 27, var
    ]

    flag_vals_1 = df.copy()[["time", var]]
    flag_vals_1.loc[:, var] = np.nan
    flag_vals_1.loc[df[var + "_eraqc"] == 28, var] = df.loc[
        df[var + "_eraqc"] == 28, var
    ]

    flag_vals_2 = df.copy()[["time", var]]
    flag_vals_2.loc[:, var] = np.nan
    flag_vals_2.loc[df[var + "_eraqc"] == 29, var] = df.loc[
        df[var + "_eraqc"] == 29, var
    ]

    # Amount of data flagged
    flag_label_0 = "Same hour replication"
    flag_label_1 = "Consecutive replication"
    flag_label_2 = "Whole-day replication"

    # if no flags are present, it messes with time x axis labels
    if len(flag_vals_0) != 0:
        flag_vals_0.plot(
            ax=ax,
            x="time",
            y=var,
            marker="s",
            ms=7,
            lw=0,
            mfc="none",
            color="C3",
            label=flag_label_0,
        )

    if len(flag_vals_1) != 0:
        flag_vals_1.plot(
            ax=ax,
            x="time",
            y=var,
            marker="x",
            ms=7,
            lw=0,
            mfc="none",
            color="C4",
            label=flag_label_1,
        )

    if len(flag_vals_2) != 0:
        flag_vals_2.plot(
            ax=ax,
            x="time",
            y=var,
            marker="o",
            ms=7,
            lw=0,
            mfc="none",
            color="C2",
            label=flag_label_2,
        )

    legend = ax.legend(loc=0, prop={"size": 8})
    network = station.split("_")[0]

    # Plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper(var)
    ylab = "{} [{}]".format(ylab, units)

    ax.set_ylabel(ylab)
    ax.set_xlabel("")

    title = "Unusual repeated streaks check: {0}".format(station)
    ax.set_title(title, fontsize=10)
    month = df.month.unique()[0]
    year = df.year.unique()[0]

    # save to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    figname = "qaqc_unusual_repeated_streaks_{0}_{1}_{2}-{3}".format(
        station, var, year, month
    )
    key = "{0}/{1}/qaqc_figs/{2}.png".format(directory, network, figname)
    img_data = BytesIO()
    fig.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(Body=img_data, ContentType="image/png", Key=key)

    if local:
        fig.savefig(
            "qaqc_figs/{}.png".format(figname),
            format="png",
            dpi=dpi,
            bbox_inches="tight",
        )

    # close figure to save memory
    plt.close()

    return None


# ============================================================================================================
## precip de-accumulation plot
def precip_deaccumulation_plot(df, flags, local=False, dpi=300):
    """
    Generate and save a precipitation de-accumulation plot with flagged data points.

    This function visualizes the de-accumulation of precipitation by plotting 
    the original accumulated precipitation and the de-accumulated values. 
    It highlights flagged oscillating or ringing values and saves the figure 
    either locally or to an AWS S3 bucket.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing the accumulated precipitation (`accum_pr`), 
        de-accumulated precipitation (`pr`), timestamps (`time`), and station information (`station`).
    flags : pandas.Series (bool)
        A boolean Series indicating flagged data points that exhibit oscillating or ringing behavior.
    local : bool, optional
        If `True`, the plot is saved locally. If `False` (default), the plot is uploaded to AWS S3.
    dpi : int, optional
        Resolution of the saved figure in dots per inch (default is 300).

    Returns
    -------
    None
        The function saves the plot but does not return any values.

    Notes
    -----
    - The top subplot shows the original accumulated precipitation (`accum_pr`).
    - The bottom subplot shows the de-accumulated precipitation (`pr`).
    - Flagged ringing values are marked in red on the accumulated precipitation plot.
    - The function automatically adjusts the y-axis limits to mitigate the effect of outliers.
    - The plot is saved either locally or to AWS S3 in the "wecc-historical-wx" bucket.
    - The `_plot_format_helper("pr")` function is used to determine y-axis labels and units.
    """
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 7))

    # Plot variable and flagged data
    df.plot(
        x="time",
        y="accum_pr",
        ax=ax0,
        marker=".",
        ms=3,
        lw=1,
        color="k",
        alpha=0.5,
        label="Original accumulated precip",
    )
    ax0.set_xticklabels([])
    df.plot(
        x="time",
        y="pr",
        ax=ax1,
        marker=".",
        ms=3,
        lw=1,
        color="C0",
        alpha=0.5,
        label="De-accumulated precip",
    )

    # Set ylims in a way we can avoid big ranges due to outliers/spikes
    mean = np.mean(df["pr"])
    std = np.std(df["pr"])
    z_scores = (df["pr"] - mean) / std
    ylim0 = -0.5
    ylim1 = np.max(df["pr"][z_scores <= 4])
    ax1.set_ylim(ylim0, ylim1)

    station = df["station"].unique()[0]
    network = station.split("_")[0]

    # Plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper("pr")
    ylab = "{0} [{1}]".format(ylab, units)
    title = "Precipitation deaccumulation: {0}".format(station)
    ax0.set_title(title, fontsize=10)

    # Plot oscillating/ringing flags
    df[flags].plot(
        x="time",
        y="accum_pr",
        ax=ax0,
        marker="o",
        mfc="none",
        lw=0,
        ms=4,
        color="red",
        alpha=0.5,
        label="Bad osscillating/ringing values",
    )

    for ax in (ax0, ax1):
        ax.set_ylabel(ylab)
        ax.set_xlabel("")
        legend = ax.legend(loc=0, prop={"size": 8})

    # save to AWS
    bucket_name = "wecc-historical-wx"
    directory = "3_qaqc_wx"
    figname = "qaqc_figs/qaqc_recip_deaccumulation_{0}".format(station)

    key = "{0}/{1}/{2}.png".format(directory, network, figname)
    img_data = BytesIO()
    fig.savefig(img_data, format="png", dpi=dpi, bbox_inches="tight")
    img_data.seek(0)
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(Body=img_data, ContentType="image/png", Key=key)

    if local:
        fig.savefig(
            "{}.png".format(figname), format="png", dpi=dpi, bbox_inches="tight"
        )

    # close figure to save memory
    plt.close("all")

    return None


# ============================================================================================================
## V2 research: these should live in qaqc_unusual_gaps.py but running into some circular import issues
def standardized_median_bounds(df, var, iqr_thresh):
    """Part 1: Calculates the standardized median.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    var : str
        variable name
    iqr_thresh : int
        min threshold multiplier for IQR to flag

    Returns
    -------
    std_med : float
        climatological median of specific month passed by df
    lower_bnd : float
        lower bound determined by IQR range and std_med
    upper_bnd : float
        upper bound determined by IQR range and std_med

    Notes
    -----
    1. v2: for some reason plotting only works if this function is in plot and not unusual gaps
    """

    std_med = df[var].median()  # climatological median for that month
    iqr = iqr_range(df, var)
    lower_bnd = std_med - (iqr_thresh * iqr)
    upper_bnd = std_med + (iqr_thresh * iqr)

    return (std_med, lower_bnd, upper_bnd)


def iqr_range(df, var):
    """Part 1: Calculates the monthly interquartile range.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    var : str
        variable name

    Returns
    -------
    range_to_return : float
        interquartile range
    """
    range_to_return = df[var].quantile([0.25, 0.75]).diff().iloc[-1]
    return range_to_return


def standardized_iqr(df, var):
    """Part 2: Standardizes data against the interquartile range.

    Parameters
    ----------
    df : pd.DataFrame
        QA/QC dataframe to produce plot on
    var : str
        variable name

    Returns
    -------
    std_data : pd.DataFrame
        standardized data based on IQR range
    """
    std_data = (df[var].values - df[var].median()) / iqr_range(df, var)
    return std_data
