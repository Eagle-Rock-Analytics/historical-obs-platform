"""
This is a script where Stage 3: QA/QC function(s) on unusually frequent values in the data observations are flagged. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import boto3
import geopandas as gp
import numpy as np
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

try:
    from qaqc_utils import *
except Exception as e:
    print("Error importing qaqc_utils: {}".format(e))


def open_log_file_logic(file):
    global log_file
    log_file = file


# #####################################
# #FOR DEBUG
# #UNCOMMENT FOR NOTEBOOK DEBUGGING
# global log_file
# log_file = open("logtest.log","w")
# verbose=True
# #####################################


# -----------------------------------------------------------------------------
## logic check: dew point must not exceed air temperature
def qaqc_crossvar_logic_tdps_to_tas_supersat(df, verbose=False):
    """
    Checks that dewpoint temperature does not exceed air temperature.
    If fails, only dewpoint temperature is flagged.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        12,qaqc_crossvar_logic_tdps_to_tas_supersat,Cross-variable logic check failure: dewpoint temperature exceeds air temperature
    """

    printf(
        "Running: qaqc_crossvar_logic_tdps_to_tas_supersat",
        log_file=log_file,
        verbose=verbose,
    )

    try:
        # first check that tdps and/or tdps_derived are provided
        dew_vars = [col for col in df.columns if "tdps" in col]
        all_dew_vars = [
            var for var in dew_vars if "qc" not in var
        ]  # remove all qc variables so they do not also run through: raw, eraqc

        # dew point is not present
        if not all_dew_vars:
            printf(
                "Station does not report dew point temperature - bypassing temperature cross-variable logic check",
                log_file=log_file,
                verbose=verbose,
            )

        # dew point is present
        else:
            for dew_var in all_dew_vars:
                # only use valid obs for both dewpoint and air temp
                df_valid = grab_valid_obs(df, var="tas", var2=dew_var)
                isBad = df_valid.loc[df_valid[dew_var] > df_valid["tas"]]
                df.loc[isBad.index, dew_var + "_eraqc"] = (
                    12  # see qaqc_flag_meanings.csv
                )

        return df

    except Exception as e:
        printf(
            "qaqc_crossvar_logic_tdps_to_tas_supersat failed with Exception: {}".format(
                e
            ),
            log_file=log_file,
            verbose=verbose,
        )
        return None


# ----------------------------------------------------------------------
def qaqc_crossvar_logic_tdps_to_tas_wetbulb(df, verbose=False):
    """
    Checks for extended periods of a dewpoint depression of 0°C.
    Extended period is defined as a 24-hour period
    If fails, only dewpoint temperature is flagged.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        13,qaqc_crossvar_logic_tdps_to_tas_wetbulb,Cross-variable logic check failure: extended streak of a zero dewpoint depression (indicative of instrument failure)
    """

    printf(
        "Running: qaqc_crossvar_logic_tdps_to_tas_wetbulb",
        log_file=log_file,
        verbose=verbose,
    )

    try:
        df_dpt = df.copy(deep=True)
        # first check that tdps and/or tdps_derived are provided
        dew_vars = [col for col in df_dpt.columns if "tdps" in col]
        all_dew_vars = [
            var for var in dew_vars if "qc" not in var
        ]  # remove all qc variables so they do not also run through: raw, eraqc

        # dew point is not present
        if not all_dew_vars:
            printf(
                "Station does not report dew point temperature - bypassing temperature cross-variable logic check",
                log_file=log_file,
                verbose=verbose,
            )

        # dew point is present
        else:
            for dew_var in all_dew_vars:
                # only use valid obs for both dewpoint and air temp
                df_valid = grab_valid_obs(df, var="tas", var2=dew_var)
                df_valid = df_valid.assign(
                    dew_depression=df_valid["tas"] - df_valid[dew_var]
                )
                df_to_check = df_valid.loc[df_valid["dew_depression"] == 0]

                # identify and flag long streak of dew point depression values = 0
                for t in df_to_check.time:
                    dpd_to_check = df_valid.loc[
                        (df_valid.time >= t)
                        & (df_valid.time <= (t + datetime.timedelta(days=1)))
                    ]["dew_depression"]

                    if all(v == 0 for v in dpd_to_check):
                        df_dpt.loc[
                            (df_dpt.time >= t)
                            & (df_dpt.time <= (t + datetime.timedelta(days=1))),
                            dew_var + "_eraqc",
                        ] = 13  # see qaqc_flag_meanings.csv

                # only print warning flag once
                if 13 in df_dpt[dew_var + "_eraqc"].unique():
                    printf(
                        "Flagging extended streak in dewpoint depression",
                        log_file=log_file,
                        verbose=verbose,
                    )

        return df_dpt

    except Exception as e:
        printf(
            "qaqc_crossvar_logic_tdps_to_tas_wetbulb failed with Exception: {}".format(
                e
            ),
            log_file=log_file,
            flush=True,
        )
        return None


# ----------------------------------------------------------------------
## logic check: precip does not have any negative values
def qaqc_precip_logic_nonegvals(df, verbose=False):
    """
    Ensures that precipitation values are positive. Negative values are flagged as impossible.
    Provides handling for the multiple precipitation variables presently in the cleaned data.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        10,qaqc_precip_logic_nonegvals,Precipitation value reported below 0 (negative value)
    """

    printf("Running: qaqc_precip_logic_nonegvals", log_file=log_file, verbose=verbose)

    df_neg_pr = df.copy(deep=True)

    # identify which precipitation vars are reported by a station
    vars_to_remove = ["qc", "duration", "method", "depth"]
    all_pr_vars = [
        var for var in df_neg_pr.columns if "pr" in var
    ]  # can be variable length depending if there is a raw qc var
    pr_vars = [
        var
        for var in all_pr_vars
        if not any(True for item in vars_to_remove if item in var)
    ]  # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    printf(
        "Running qaqc_precip_logic_nonegvals on: {}".format(pr_vars),
        log_file=log_file,
        verbose=verbose,
    )

    try:
        if not pr_vars:  # precipitation variable(s) is not present
            printf(
                "Station does not report precipitation - bypassing precip logic nonnegvals check",
                log_file=log_file,
                verbose=verbose,
            )
        else:
            for item in pr_vars:
                df_valid = grab_valid_obs(df_neg_pr, item)  # subset for valid obs
                df_valid.loc[df_valid[item] < 0, item + "_eraqc"] = (
                    10  # see era_qaqc_flag_meanings.csv
                )

        return df_valid

    except Exception as e:
        printf(
            "qaqc_precip_logic_nonegvals failed with Exception: {0}".format(e),
            log_file=log_file,
            verbose=verbose,
        )
        return None


# ----------------------------------------------------------------------
## logic check: precip accumulation amounts balance for time period
def qaqc_precip_logic_accum_amounts(df, verbose=False):
    """
    Ensures that precipitation accumulation amounts are consistent with reporting time frame.
    Only needs to be applied when 2 or more precipitation duration specific variables are present (pr_5min, pr_1h, pr_24h)
    For example: pr_5min should not be larger than pr_1h

    Rules:
    ------
        1) pr_5min < pr_1h < pr_24h
        2) pr_localmid should never exceed pr_24h

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        16,qaqc_precip_logic_accum_amounts,Cross-variable logic check failure: accumulated precipitation value in shorter window is larger than in longer window (e.g. pr_5min > pr_1h)
        17,qaqc_precip_logic_accum_amounts,Cross-variable logic check failure: accumulated precipitation value in longer window is smaller than in shorter window (e.g. pr_24h < pr_1h)
        18,qaqc_precip_logic_accum_amounts,Cross-variable logic check failure: accumulated precipitation in a 24h period is too low compared to accumulated precipitation since local midnight
    """
    printf(
        "Running: qaqc_precip_logic_accum_amounts", log_file=log_file, verbose=verbose
    )

    # identify which precipitation vars are reported by a station
    vars_to_remove = ["qc", "duration", "method", "depth"]
    all_pr_vars = [
        var for var in df.columns if "pr" in var
    ]  # can be variable length depending if there is a raw qc var
    pr_vars = [
        var
        for var in all_pr_vars
        if not any(True for item in vars_to_remove if item in var)
    ]  # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process

    try:
        # if station does not report any precipitation values, or only one, bypass
        if len(pr_vars) == 0 or len(pr_vars) == 1:
            printf(
                "Station does not report multiple precipitation variables - bypassing precip logic accum check",
                log_file=log_file,
                verbose=verbose,
            )
            return df

        # checks accumulated precip vars against each other
        # noting that these flags are essentially identical in operation
        # flag assignment is logically dependent on the first var to determine
        # which flag is placed (i.e. to determine if too larges/small)

        # checks accumulated precip vars against each other
        # noting that these flags are essentially identical in operation
        # flag assignment is logically dependent on the first var to determine which flag is placed (i.e. to determine if too larges/small)

        if "pr_5min" in pr_vars:
            if "pr_1h" in pr_vars:
                # only use valid obs for precip vars
                df_valid = grab_valid_obs(df, var="pr_5min", var2="pr_1h")
                ind = (df_valid["pr_5min"] > df_valid["pr_1h"]).index
                df.loc[ind, "pr_5min_eraqc"] = 16  # see era_qaqc_flag_meanings.csv

            if "pr_24h" in pr_vars:
                df_valid = grab_valid_obs(df, var="pr_5min", var2="pr_24h")
                ind = (df_valid["pr_5min"] > df_valid["pr_24h"]).index
                df.loc[ind, "pr_5min_eraqc"] = 16  # see era_qaqc_flag_meanings.csv

        if "pr_1h" in pr_vars:
            if "pr_5min" in pr_vars:
                df_valid = grab_valid_obs(df, var="pr_5min", var2="pr_1h")
                ind = (df_valid["pr_1h"] < df_valid["pr_5min"]).index
                df.loc[ind, "pr_1h_eraqc"] = 17  # see era_qaqc_flag_meanings.csv

            if "pr_24h" in pr_vars:
                df_valid = grab_valid_obs(df, var="pr_1h", var2="pr_24h")
                ind = (df_valid["pr_1h"] > df_valid["pr_24h"]).index
                df.loc[ind, "pr_1h_eraqc"] = 17  # see era_qaqc_flag_meanings.csv

        if "pr_24h" in pr_vars:
            if "pr_5min" in pr_vars:
                df_valid = grab_valid_obs(df, var="pr_5min", var2="pr_24h")
                ind = (df_valid["pr_24h"] < df_valid["pr_5min"]).index
                df.loc[ind, "pr_24h_eraqc"] = 17  # see era_qaqc_flag_meanings.csv

            if "pr_1h" in pr_vars:
                df_valid = grab_valid_obs(df, var="pr_1h", var2="pr_24h")
                ind = (df_valid["pr_24h"] < df_valid["pr_1h"]).index
                df.loc[ind, "pr_24h_eraqc"] = 17  # see era_qaqc_flag_meanings.csv

            if "pr_localmid" in pr_vars:
                df_valid = grab_valid_obs(df, var="pr_localmid", var2="pr_24h")
                ind = (df_valid["pr_24h"] < df_valid["pr_localmid"]).index
                df.loc[ind, "pr_24h_eraqc"] = 18  # see era_qaqc_flag_meanings.csv

        return df

    except Exception as e:
        printf(
            "qaqc_precip_logic_accum_amounts failed with Exception: {0}".format(e),
            log_file=log_file,
            verbose=verbose,
        )
        return None


# ----------------------------------------------------------------------
## logic check: wind direction must be 0 if wind speed is 0
def qaqc_crossvar_logic_calm_wind_dir(df, verbose=False):
    """
    Checks that wind direction is zero when wind speed is also zero.
    If fails, wind direction is flagged.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        14,qaqc_crossvar_logic_calm_wind_dir,Cross-variable logic check failure: wind direction is not zero when wind speed is zero
        15,qaqc_crossvar_logic_calm_wind_dir,Cross-variable logic check failure: wind direction manually reset to 360 to represent true northerly winds
    """

    # import pdb; pdb.set_trace()
    printf(
        "Running: qaqc_crossvar_logic_calm_wind_dir", log_file=log_file, verbose=verbose
    )

    try:
        # Noting that a wind direction value of 0 is a valid value
        # Only a problem when wind speed is also 0, where 0 now means no winds for there to be a direction

        # check that wind direction is provided
        if "sfcWind_dir" not in df.columns:
            printf(
                "Station does not report wind direction - bypassing wind cross-variable logic check",
                log_file=log_file,
                verbose=verbose,
            )
            return df
        elif "sfcWind_dir" in df.columns and "sfcWind" not in df.columns:
            printf(
                "Station does reports wind direction, but not wind speed - bypassing wind cross-variable logic check",
                log_file=log_file,
                verbose=verbose,
            )
            return df

        # use only valid observations
        df_valid = grab_valid_obs(df, var="sfcWind", var2="sfcWind_dir")

        # identify calm winds but with incorrect wind directions
        isBad = df_valid.loc[
            (df_valid["sfcWind"] == 0)
            & (df_valid["sfcWind_dir"] != 0)
            & ~(df_valid["sfcWind_dir"]).isnull()
            == True
        ]
        df.loc[isBad.index, "sfcWind_dir_eraqc"] = 14  # see qaqc_flag_meanings.csv

        # identify non-zero winds but with incorrect wind directions
        # non-zero northerly winds should be coded as 360 deg, not 0 deg
        isBad = df_valid.loc[
            (df_valid["sfcWind"] != 0) & (df_valid["sfcWind_dir"] == 0)
        ]
        df.loc[isBad.index, "sfcWind_dir"] = 360
        df.loc[isBad.index, "sfcWind_dir_eraqc"] = 15  # see qaqc_flag_meanings.csv

        return df

    except Exception as e:
        printf(
            "qaqc_crossvar_logic_calm_wind_dir failed with Exception: {}".format(e),
            log_file=log_file,
            verbose=verbose,
        )
        return None


# -----------------------------------------------------------------------------
## temporary fix on pressure variables being in the wrong unit
## fn to be removed from pipeline on next full cleaning update
def qaqc_pressure_units_fix(df, verbose=False):
    """
    Ensures that stations consistently report pressure vars in Pa units. This largely impacts ASOSAWOS stations,
    where the pressure unit conversion did not take.

    This is a temporary fix; in the next cleaning update, unit conversions will be applied and checked.
    No flag is placed in this fix, if variable fails it will be caught by the world records check.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Notes:
    ------
    grab_valid_data is not applied here as this is a temporary fix, applied uniformly
    """

    printf("Running: qaqc_pressure_units_fix", log_file=log_file, verbose=verbose)

    try:
        # identify pressure variables to check conversion on
        ps_vars = ["ps", "psl", "ps_altimeter", "ps_derived"]

        for var in ps_vars:
            if var in df.columns:
                if df[var].mean() < 10000:
                    df[var] = df[var] * 100.0
                    printf(
                        "Pressure units on {} updated to be Pa".format(var),
                        log_file=log_file,
                        verbose=verbose,
                    )
        return df

    except Exception as e:
        printf(
            "qaqc_pressure_units_fix failed with Exception: {}".format(e),
            log_file=log_file,
            verbose=verbose,
        )
        return None
