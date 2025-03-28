"""
This is a script where Stage 3: QA/QC function(s) on spurious buoy issues within the NDBC and MARITIME networks are flagged. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import numpy as np
import datetime

# New logger function
from log_config import logger

# Import utils
try:
    from qaqc_utils import *
except Exception as e:
    logger.debug("Error importing qaqc_utils: {}".format(e))


## NDBC and MARITIME only
# -----------------------------------------------------------------------------
def spurious_buoy_check(df, qc_vars, verbose=False):
    """
    Checks the end date on specific buoys to confirm disestablishment/drifting dates of coverage.
    If station reports data past disestablishment date, data records are flagged as suspect.

    Parameters
    ----------
    df : pd.DataFrame
        QAQC dataframe to check
    qc_vars : list of str
        QC'd variables to avoid
    verbose : bool, optional
        If True, provides runtime output to local terminal

    Returns
    -------
    df : pd.DataFrame
        QAQC dataframe post-check

    Notes
    -----
    Flag meaning : 1,spurious_buoy_check,Suspect observation (i.e. buoy reports wind during day with known proximity issue)
    Flag meaning : 2,spurious_buoy_check,Out of range for station official data temporal record (i.e. buoy reports data past disestablishment date)
    """
    logger.info("Running: spurious_buoy_check")

    # buoys with known issues, or potential issues
    known_issues = [
        "NDBC_46023",
        "NDBC_46045",
        "NDBC_46051",
        "NDBC_46044",
        "MARITIME_PTAC1",
        "MARITIME_PTWW1",
        "MARITIME_MTYC1",
        "MARITIME_MEYC1",
        "MARITIME_SMOC1",
        "MARITIME_ICAC1",
    ]
    potential_issues = [
        "NDBC_46290",
        "NDBC_46404",
        "NDBC_46212",
        "NDBC_46216",
        "NDBC_46220",
        "NDBC_46226",
        "NDBC_46227",
        "NDBC_46228",
        "NDBC_46230",
        "NDBC_46234",
        "NDBC_46245",
        "NDBC_46250",
    ]

    # extract station name
    station = df["station"].unique()[0]

    if station in known_issues:
        logger.info(
            "{0} has a known issue, checking data coverage.".format(station),
        )

        # buoys with "data" past their disestablishment dates
        if station == "NDBC_46023":  # disestablished 9/8/2010
            isBad = df.loc[df["time"] >= np.datetime64("2010-09-09")]
            for new_var in qc_vars:
                if new_var != "elevation_qaqc":
                    df.loc[df.time.isin(isBad.time), new_var] = (
                        2  # see era_qaqc_flag_meanings.csv
                    )

        elif station == "NDBC_46045":  # disestablished 11/1997
            isBad = df.loc[df["time"] >= np.datetime64("1997-12-01")]
            for new_var in qc_vars:
                if new_var != "elevation_qaqc":
                    df.loc[df.time.isin(isBad.time), new_var] = (
                        2  # see era_qaqc_flag_meanings.csv
                    )

        elif (
            station == "NDBC_46051"
        ):  # disestablished 4/1996, and out of range of DEM (past coastal range) but reports nan elevation
            isBad = df.loc[df["time"] >= np.datetime64("1996-05-01")]
            for new_var in qc_vars:
                if new_var != "elevation_qaqc":
                    df.loc[df.time.isin(isBad.time), new_var] = (
                        2  # see era_qaqc_flag_meanings.csv
                    )

        elif (
            station == "MARITIME_PTAC1"
        ):  # data currently available 1984-2012, but disestablished 2/9/2022
            # only flag if new data is added after 2022 in a new data pull
            isBad = df.loc[df["time"] >= np.datetime64("2022-02-09")]
            for new_var in qc_vars:
                if new_var != "elevation_qaqc":
                    df.loc[df.time.isin(isBad.time), new_var] = (
                        2  # see era_qaqc_flag_meanings.csv
                    )

        # adrift buoy that reports valid data during adrift period (5/2/2015 1040Z to 5/3/2015 1600Z)
        elif station == "NDBC_46044":
            isBad = df.loc[
                (df["time"] >= np.datetime64("2015-05-02 10:40:00"))
                & (df["time"] <= np.datetime64("2015-05-03 15:50:00"))
            ]
            for new_var in qc_vars:
                if new_var != "elevation_qaqc":
                    df.loc[df.time.isin(isBad.time), new_var] = (
                        2  # see era_qaqc_flag_meanings.csv
                    )

        # other known issues
        elif (
            station == "MARITIME_PTWW1"
        ):  # wind data obstructed by ferries docking at pier during day hours
            # only wind vars need flag during "day" hours, currently set for 6am to 8pm every day
            df_valid = grab_valid_obs(df, var="sfcWind", var2="sfcWind_dir")
            isBad = df_valid.loc[(df_valid["hour"] >= 6) & (df_valid["hour"] <= 20)]
            df.loc[df.time.isin(isBad.time), "sfcWind_eraqc"] = (
                1  # see era_qaqc_flag_meanings.csv
            )
            df.loc[df.time.isin(isBad.time), "sfcWind_dir_eraqc"] = (
                1  # see era_qaqc_flag_meanings.csv
            )

        # elif station == "MARITIME_MTYC1" or station == "MARITIME_MEYC1": # buoy was renamed, no relocation; MTYC1 2005-2016, MEYC1 2016-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests

        # elif station == "MARITIME_SMOC1" or station == "MARITIME_ICAC1": # buoy was renamed, small relocation (see notes); SMOC1 2005-2010, ICAC1 2010-2021
        #     # modify attribute/naming with note
        #     # this will get flagged in station proximity tests

    elif station in potential_issues:
        # other stations have partial coverage of their full data records as well as disestablishment dates
        # if new data is added in the future, needs a manual check and added to known issue list if requires handling
        # most of these should be caught by not having a cleaned data file to begin with, so if this print statement occurs it means new raw data was cleaned and added to 2_clean_wx/
        if verbose:
            logger.info(
                "{0} has a reported disestablishment date, requires manual confirmation of dates of coverage.".format(
                    station
                )
            )

        for new_var in qc_vars:
            if new_var != "elevation_qaqc":
                df.loc[:, new_var] = 2  # see era_qaqc_flag_meanings.csv

    return df
