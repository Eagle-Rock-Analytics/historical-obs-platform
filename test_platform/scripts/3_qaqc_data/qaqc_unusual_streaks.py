"""
This is a script where Stage 3: QA/QC related common functions, conversions, and operations is stored for ease of use
for the Historical Observations Platform.
"""

## Import Libraries
import boto3
import numpy as np
import pandas as pd
import requests
import urllib
import xarray as xr
import datetime
import matplotlib.pyplot as plt
from io import BytesIO, StringIO

try:
    from qaqc_plot import *
except:
    print("Error importing qaqc_plot.py")
## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client("s3")  # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = (
    "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
)
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

try:
    from qaqc_utils import *
except Exception as e:
    print("Error importing qaqc_utils: {}".format(e))


def open_log_file_streaks(file):
    global log_file
    log_file = file


# #####################################
# #FOR DEBUG
# #UNCOMMENT FOR NOTEBOOK DEBUGGING
# global log_file
# log_file = open("logtest.log","w")
# verbose=True
# #####################################


# ----------------------------------------------------------------------
def infere_freq(df):
    """ """
    # Calculate time differences in data
    time_diff = pd.Series(df["time"]).diff()

    # Count resolution occurences and use only those that occure in 5% of series or more
    frequencies = time_diff.value_counts()
    frequencies = frequencies[frequencies > len(df) // 5]

    # Create resolutions dictionary
    frequencies = {
        np.round(n / len(df), decimals=5): float(f / 1e9)
        for n, f in zip(frequencies.values, frequencies.index.values)
    }
    return frequencies


# ----------------------------------------------------------------------
def infere_res_var(df, var):
    """ """
    # If var is tdps_derived, use temp measurements
    # if var == "tdps_derived":
    #     var = "tas"

    # Extract var data
    data = df.copy()[var]

    # TODO: use function from calc_clean.py
    # Convert from Pa to hPa -- only to determine mode and resolution of data, data not returned in hPa
    if data.mean() > 10000 and (
        var == "ps" or var == "psl" or var == "ps_altimeter" or var == "ps_derived"
    ):
        data = data / 100

    data_diff = data.sort_values().diff().clip(lower=0.01, upper=1).dropna()
    if len(data_diff) <= 10:
        return 0.5
    else:
        # Calculate modified mode from avg mode and median
        # mode0 = data.sort_values().diff().replace(0, pd.NaT).mode().values[0]
        # mode1 = data.sort_values().diff().replace(0, pd.NaT).median()
        mode0 = data.sort_values().diff().mode().values[0]
        mode1 = data.sort_values().diff().median()
        mode = (mode0 + mode1) / 2

        # Round to the nearest 0.5
        multiplied = round(mode * 2)
        rounded_to_whole = multiplied / 2

        # If mode is 0.25 or less, round to 0.1
        if rounded_to_whole <= 0.25:
            mode = 0.1
        else:
            mode = rounded_to_whole

        if mode <= 1:
            return mode
        else:
            return 1.0


# ----------------------------------------------------------------------
def infere_res(df, verbose=False):
    """ """
    check_vars = [
        "tas",
        "tdps",
        "tdps_derived",
        "ps",
        "psl",
        "ps_derived",
        "ps_altimeter",
        "sfcWind",
    ]
    variables = [var for var in check_vars if var in df.columns]
    # variables = [var for var in df.columns if any(True for item in check_vars if item in var)]
    # printf(variables, log_file=log_file, verbose=verbose)

    resolutions = {}
    for var in variables:
        # if var in df.columns:
        if df[var].isnull().all() or len(np.where(df[var].isnull())[0]) < 100:
            resolutions[var] = 0.1
        else:
            resolutions[var] = infere_res_var(df, var)

    return resolutions


# ----------------------------------------------------------------------
# Straight repeat streak criteria
straight_repeat_criteria = {
    "tas": {
        1: [40, 14],  # 40 values of 14 days
        0.5: [30, 10],  # 30 values or 10 days
        0.1: [24, 7],  # 24 values or 7 days
    },
    "tdps": {1: [80, 14], 0.5: [60, 10], 0.1: [48, 7],},  # of  # or  # or
    "psl": {1: [120, 28], 0.5: [100, 21], 0.1: [72, 14]},  # of  # or  # or
    "sfcWind": {1: [40, 14], 0.5: [30, 10], 0.1: [24, 7],},  # of  # or  # or
}
straight_repeat_criteria["tdps_derived"] = straight_repeat_criteria["tdps"]
straight_repeat_criteria["ps"] = straight_repeat_criteria["psl"]
straight_repeat_criteria["ps_derived"] = straight_repeat_criteria["psl"]
straight_repeat_criteria["ps_altimeter"] = straight_repeat_criteria["psl"]

# ----------------------------------------------------------------------
# Hour repeat streak criteria
hour_repeat_criteria = {
    "tas": {1: 25, 0.5: 20, 0.1: 15,}  # 40 days  # 20 days  # 15 days
}
# All variables have the same hourly criteria
hour_repeat_criteria["tdps"] = hour_repeat_criteria["tas"]
hour_repeat_criteria["tdps_derived"] = hour_repeat_criteria["tas"]
hour_repeat_criteria["psl"] = hour_repeat_criteria["tas"]
hour_repeat_criteria["ps"] = hour_repeat_criteria["tas"]
hour_repeat_criteria["ps_altimeter"] = hour_repeat_criteria["tas"]
hour_repeat_criteria["ps_derived"] = hour_repeat_criteria["tas"]
hour_repeat_criteria["sfcWind"] = hour_repeat_criteria["tas"]

# ----------------------------------------------------------------------
# Day repeat streak criteria
day_repeat_criteria = {"tas": {1: 10, 0.5: 7, 0.1: 5,}}  # 10 days  #  7 days  #  5 days
# All variables have the same daily criteria
day_repeat_criteria["tdps"] = day_repeat_criteria["tas"]
day_repeat_criteria["tdps_derived"] = day_repeat_criteria["tas"]
day_repeat_criteria["psl"] = day_repeat_criteria["tas"]
day_repeat_criteria["ps"] = day_repeat_criteria["tas"]
day_repeat_criteria["ps_altimeter"] = day_repeat_criteria["tas"]
day_repeat_criteria["ps_derived"] = day_repeat_criteria["tas"]
day_repeat_criteria["sfcWind"] = day_repeat_criteria["tas"]

# ----------------------------------------------------------------------
# Min wind value for straight repeat test
# TODO: HadISD thresholds: does it make sense to change them in future versions?
# More analysis needs to be done to ensure what is a good threshold for calm wind conditions for this test
WIND_MIN_VALUE = {1: 1.0, 0.5: 0.5, 0.1: 0.5}


# ---------------------------------------------------------------------------------------------------
# Function to create a new column for consecutive months
def consecutive_months(series):

    indices = np.where(np.diff(series.values) > 1)[0] + 1
    clusters = np.split(series.values, indices)
    isin = [series.isin(c) for c in clusters]
    groups = np.zeros_like(series.values, dtype="int")

    for i, ind in enumerate(isin):
        groups[ind.values] = int(i)
    return pd.Series(groups, index=series.index)


# ---------------------------------------------------------------------------------------------------
def qaqc_unusual_repeated_streaks(
    df, plot=True, local=False, verbose=False, min_sequence_length=10
):
    """
    Test for repeated streaks/unusual spell frequenc.
    Three test are conducted here:
       - Consecutive observation replication
       - Same hour observation replication over a number of days
         (either using a threshold of a certain number of observations,
          or for sparser records, a number of days during which all the
          observations have the same value)
       - Whole day replication for a streak of days

    This test is done for ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter" "sfcWind"]

    Input:
    -----
        df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
        local [bool] : if True, saves the plot to local directory
        plot [bool] : if True, produces plot and uploads it to AWS
    Output:
    ------
        if qaqc succeded:
            df [pandas dataframe] : QAQC dataframe with flagged values (see below for flag meaning).
        else if qaqc failed:
            None
    Flag meaninig:
    -------------
        27,qaqc_unusual_repeated_streaks,Same hour observation replication over a number of days
        28,qaqc_unusual_repeated_streaks,Same hour observation replication over a number of days or number of observations
        29,qaqc_unusual_repeated_streaks,Whole day replication for a streak of days

    NOTES:
    Threshold for different variables/resolutions are noted on:
    https://doi.org/10.5194/cp-8-1649-2012 : Table 4

    TODO:
    min_sequence_length can be tweaked, althought HadISD uses 10
    """
    printf("Running: qaqc_unusual_repeated_streaks", log_file=log_file, verbose=verbose)

    station = df["station"].dropna().unique()[0]

    # import pdb; pdb.set_trace()
    try:
        # Infere resolution from data
        resolutions = infere_res(df)

        # Save original df multiindex and create station column
        new_df = df.copy()

        # Define test variables and check if they are in the dataframe
        check_vars = [
            "tas",
            "tdps",
            "tdps_derived",
            "ps",
            "psl",
            "ps_derived",
            "ps_altimeter",
            "sfcWind",
        ]
        variables = [var for var in check_vars if var in new_df.columns]
        printf(
            "Running {} on {}".format("qaqc_unusual_repeated_streaks", variables),
            verbose=verbose,
            log_file=log_file,
            flush=True,
        )

        # Loop through test variables
        for var in variables:
            printf(
                "Running unusual streaks check on: {}".format(var),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            # Create a copy of the original dataframe and drop NaNs in the testing variable
            test_df = new_df.copy().dropna(subset=var)

            # Use only values that have not been flagged by previous QAQC tests
            test_df = grab_valid_obs(test_df, var)  # subset for valid obs

            # first scans suspect values using entire record
            if test_df[var].isna().all() == True:
                printf(
                    "All values for {} are flagged, bypassing qaqc_unusual_repeated_streaks".format(
                        var
                    ),
                    verbose=verbose,
                    log_file=log_file,
                    flush=True,
                )
                continue  # bypass to next variable if all obs are nans

            # Choose resolution
            res = resolutions[var]

            ##########################################################################################
            ## NOTE for V2:
            ## tdps_derived (and probably other derived quantities) have a much higher resolution than
            ## what the measured variable would have, since it's the combination of two or more variables
            ## This makes the computation of streaks much slower, since the grouping by var value
            ## will create much more groups. This can be avoided by rounding to a decimal approximation
            ## of the variable. From my tests (Héctor) this is 1 decimal place, the number of groups
            ## decreases, the computation time decreases, and the number of streaks does not increase,
            ## so it does not create artifact streaks
            ##
            ## test = stn_to_qaqc.copy()
            ## test.loc[:,'tdps_derived'] = test['tdps_derived'].round(decimals=1)
            ##
            ##########################################################################################
            # ------------------------------------------------------------------------------------------------
            # Hour repeat streak criteria
            printf(
                "Running hourly repeats on {}".format(var),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            tt00 = time.time()
            threshold = hour_repeat_criteria[var][res]
            bad_hourly = hourly_repeats(
                test_df, var=var, threshold=threshold
            )  # Bad hourly returns a pd.Series of time stamps
            new_df.loc[
                new_df["time"].isin(bad_hourly), var + "_eraqc"
            ] = 27  # Flag _eraqc variable
            printf(
                "Hourly repeats flagged for {}. Ellapsed time: {:.2f}".format(
                    var, time.time() - tt00
                ),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            # ------------------------------------------------------------------------------------------------
            # Straight repeat streak criteria
            printf(
                "Running straight repeats on {}".format(var),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            tt00 = time.time()
            threshold = straight_repeat_criteria[var][res]
            if var == "sfcWind":
                wind_min_value = WIND_MIN_VALUE[res]
            else:
                wind_min_value = None
            bad_straight = consecutive_repeats(
                test_df,
                var,
                threshold,
                wind_min_value,
                min_sequence_length=min_sequence_length,
            )  # Bad straight returns a pd.Series of time stamps
            new_df.loc[
                new_df["time"].isin(bad_straight), var + "_eraqc"
            ] = 28  # Flag _eraqc variable
            printf(
                "Straight repeats flagged for {}. Ellapsed time: {:.2f}".format(
                    var, time.time() - tt00
                ),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            # ------------------------------------------------------------------------------------------------
            # Whole day replication for a streak of days
            printf(
                "Running whole day repeats on {}".format(var),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            tt00 = time.time()
            threshold = day_repeat_criteria[var][res]
            bad_whole = consecutive_fullDay_repeats(
                test_df, var, threshold
            )  # Bad whole returns a pd.Series of time stamps
            new_df.loc[
                new_df["time"].isin(bad_whole), var + "_eraqc"
            ] = 29  # Flag _eraqc variable
            printf(
                "Whole day repeats flagged for {}. Ellapsed time: {:.2f}".format(
                    var, time.time() - tt00
                ),
                verbose=verbose,
                log_file=log_file,
                flush=True,
            )
            # ------------------------------------------------------------------------------------------------
            #
            bad_keys = np.concatenate(
                [
                    len(bad_hourly) * [27],
                    len(bad_straight) * [28],
                    len(bad_whole) * [29],
                ]
            )
            bad_times = np.concatenate(
                [bad_hourly.values, bad_straight.values, bad_whole.values]
            )
            bad = pd.DataFrame(data={"time": bad_times, "flag": bad_keys})
            bad["year"] = bad.time.dt.year
            bad["month"] = bad.time.dt.month

            # # Groups for zoom plots
            # bad_times = df.groupby([pd.Grouper(key="time", freq='W')])['time'].min() # start of week
            # bad_counts = df.groupby([pd.Grouper(key="time", freq='W')])[var + '_eraqc'].count() # how many non-flagged counts

            # # trim out any 0 counts for only flagged weeks
            # bad = pd.DataFrame({"dt":bad_times, "count":bad_counts})
            # bad_to_run = bad.loc[(bad['count']!=0)]

            # --------------------------------------------------------
            if plot:
                ## Plotting by month/year will reduce the number of plots
                keys = bad.groupby(["year", "month"]).groups.keys()
                for k in keys:
                    ind = np.logical_and(
                        new_df["year"] == k[0], new_df["month"] == k[1]
                    )
                    unusual_streaks_plot(new_df[ind], var, station="test", local=local)
                printf(
                    "{} subset plots produced for flagged obs in {}".format(
                        len(keys), var
                    ),
                    verbose=verbose,
                    log_file=log_file,
                    flush=True,
                )

        return new_df
    except Exception as e:
        printf(
            "qaqc_unusual_repeated_streaks failed with Exception: {}".format(e),
            verbose=verbose,
            log_file=log_file,
            flush=True,
        )
        return None


# ---------------------------------------------------------------------------------------------------
def find_date_clusters(dates, threshold):
    """ """
    # Ensure the dates are sorted
    dates = pd.Series(dates).sort_values().reset_index(drop=True)

    # Calculate the difference between consecutive dates
    if dates.size > 0:
        diff = dates.diff().dt.days.fillna(1)

        # Identify clusters by assigning a cluster number to each date
        cluster_id = (diff > 1).cumsum()

        # Group the dates by their cluster id and filter by cluster size
        clusters = dates.groupby(cluster_id).filter(lambda x: len(x) > threshold)

        # Create a list of clusters
        cluster_list = [group.tolist() for _, group in clusters.groupby(cluster_id)]

        if len(cluster_list) > 0:
            return np.concatenate(cluster_list)
        else:
            return np.nan
    else:
        return np.nan


# ---------------------------------------------------------------------------------------------------
def hourly_repeats(df, var, threshold):
    """ """
    ##########################################################################################
    ## NOTE for V2:
    ## when selectin original data for a specific hour, it is possible that that series
    ## will have small (or large, which are not a concern here, although how large is large)
    ## gaps in the time. For now, this hourly test will consider same hour consecutive
    ## streaks as data separated by 1 day. IF there is a gap (for example:
    ## [1-1-1980 12:00, 1-1-1980 13:00, 1-1-1980 15:00, ...], the gap between 13 and 15 hrs
    ## would prevent the function to flag this. Or if the streak is sufficiently large before
    ## and/or after the gap, it would flag before or after the gap
    ## This could be addressed by regularization to hourly data before running this test
    ## or by interpolating "small" gaps within this test
    ##########################################################################################
    hourly_streaks = []
    values = []
    for hour in range(24):
        da = df[df["hour"] == hour]
        streaks = (
            da.groupby(var, group_keys=True)["time"]
            .apply(find_date_clusters, threshold=15)
            .dropna()
        )
        if streaks.size > 0:
            for ind in streaks.index:
                streaks_dates = list(streaks.loc[ind])
                hourly_streaks.extend(streaks_dates)
                values.extend(len(streaks_dates) * [ind])
    return pd.Series(hourly_streaks, index=values, dtype="datetime64[ns]").sort_values()


# ---------------------------------------------------------------------------------------------------
def consecutive_repeats(
    df, var, threshold, wind_min_value=None, min_sequence_length=10
):
    """
    Consecutive observation replication
         (either using a threshold of a certain number of observations,
          or for sparser records, a number of days during which all the
          observations have the same value)

    This test is done for ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter", "sfcWind"]

    Input:
    -----
        df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
        var [str] : variable to test
        threshold (int,int) : comes from straight_repeat_criteria[var][res]
    Output:
    ------
        bad [numpy array] : dates that mark the flagged values (from df.index)

    NOTES (TODO:)
    """

    da = df.copy()[[var, "time"]]

    # If variable is wind, only use values above min wind value
    if wind_min_value is not None:
        da = da[da[var] > wind_min_value]
    # Identify sequences of similar values
    da.loc[:, "group"] = (da[var] != da[var].shift()).cumsum()
    da.loc[:, "start_date"] = da["time"].values
    da.loc[:, "end_date"] = da["time"].values

    start_date = da.copy().groupby([var, "group"]).min().sort_values(by="group").copy()
    end_date = da.copy().groupby([var, "group"]).max().sort_values(by="group").copy()

    # Calculate the length of each sequence
    sequence_lengths = (
        da.copy()
        .groupby([var, "group"])
        .size()
        .reset_index(name="sequence_length")
        .sort_values(by="group")
    )
    sequence_lengths.loc[:, "start_date"] = start_date.loc[:, "start_date"].values
    sequence_lengths.loc[:, "end_date"] = end_date.loc[:, "end_date"].values

    # Filter sequences with a minimum length
    min_sequence_length = 10
    filtered_sequences = sequence_lengths.copy()[
        sequence_lengths.loc[:, "sequence_length"] >= min_sequence_length
    ]
    filtered_sequences.loc[:, "dt"] = (
        filtered_sequences.loc[:, "end_date"] - filtered_sequences.loc[:, "start_date"]
    ).values.astype("timedelta64[D]")

    # Get limits
    nvalues, ndays = threshold

    # TODO: fix this commented line, 2x sequence_length to test, it should be dt
    condition = np.logical_or(
        filtered_sequences.loc[:, "sequence_length"] > nvalues,
        filtered_sequences["dt"].astype("int") > ndays,
    )

    # Find bad groups and index in the original dataset
    bad_groups = filtered_sequences[condition].group.values
    bad = da[da["group"].isin(bad_groups)].copy()

    if len(bad) > 0:
        bad["month"] = pd.to_datetime(bad.loc[:, "time"]).dt.month.values
        bad["year"] = pd.to_datetime(bad.loc[:, "time"]).dt.year.values
    else:
        bad = pd.DataFrame(
            {
                "time": np.array([], dtype="datetime64[ns]"),
                "month": np.array([], dtype=np.int64),
                "year": np.array([], dtype=np.int64),
            }
        )
    if bad.size > 0:
        return pd.Series(bad["time"].values, index=bad[var]).sort_values()
    else:
        return pd.Series([], dtype="datetime64[ns]")
    # return bad[['time', 'month', 'year']]


# ---------------------------------------------------------------------------------------------------
def full_day_compare(series0, series1):

    ind = []
    groups = []
    g = 0
    for a, b in zip(series0.values, series1.values):
        if type(a) == np.ndarray and type(b) == np.ndarray and len(a) == len(b):
            if (a == b).all():
                groups.append(g)
            else:
                groups.append(-1)
                g += 1
        else:
            groups.append(-1)
            g += 1
    if len(np.unique(groups)) > 1:
        return np.array(groups) - np.max(groups)
    else:
        return groups


# ---------------------------------------------------------------------------------------------------
def consecutive_fullDay_repeats(df, var, threshold):
    """
    Consecutive observation replication
         (either using a threshold of a certain number of observations,
          or for sparser records, a number of days during which all the
          observations have the same value)

    This test is done for ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter", "sfcWind"]

    Input:
    -----
        df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
        var [str] : variable to test
        threshold (int,int) : comes from straight_repeat_criteria[var][res]
    Output:
    ------
        bad [numpy array] : dates that mark the flagged values (from df.index)
                  (A==B).all()
    NOTES (TODO:)
    """

    # Temporary dataframe to work with
    da = df.copy()[[var, "time", "year", "month", "day", "hour"]].dropna()
    datavar = da.groupby(by=["year", "month", "day", "hour"])[var].apply(np.nanmean)
    othervars = (
        da.drop(columns=var).groupby(by=["year", "month", "day", "hour"]).first()
    )
    data = {
        "hour": othervars.index.get_level_values(-1),
        var: datavar.values,
        "time": othervars.time.values,
    }
    da = pd.DataFrame(data)
    da["date"] = pd.to_datetime(da["time"]).dt.date.values

    # Whole days to analysis
    whole_days = da.groupby(by=["date"])[var].apply(
        lambda x: np.round(x.values, decimals=1)
    )
    whole_days = pd.DataFrame({var: whole_days, "date": whole_days.index.values})
    whole_days["group"] = full_day_compare(whole_days[var], whole_days[var].shift())

    sequence_lengths = (
        whole_days.groupby(["group"])
        .size()
        .reset_index(name="length")
        .sort_values(by="group")
    )
    sequence_lengths = sequence_lengths[sequence_lengths["group"] >= 0]

    # for g,l in zip(sequence_lengths['group'].values, sequence_lengths['length'].values):
    #     print(g,l)
    # import pdb; pdb.set_trace()
    bad_groups = np.array(
        [
            g
            for g, l in zip(
                sequence_lengths["group"].values, sequence_lengths["length"].values
            )
            if l > threshold
        ]
    )

    bad_dates = whole_days[whole_days["group"].isin(bad_groups)][["date", "group"]]

    df["date"] = pd.to_datetime(df["time"]).dt.date.values
    bad = df.copy()[df["date"].isin(bad_dates["date"])]
    bad["group"] = [bad_dates["group"].loc[d] for d in bad["date"]]

    if len(bad) > 0:
        bad["month"] = pd.to_datetime(bad.loc[:, "time"]).dt.month.values
        bad["year"] = pd.to_datetime(bad.loc[:, "time"]).dt.year.values
    else:
        bad = pd.DataFrame(
            {
                "time": np.array([], dtype="datetime64[ns]"),
                "month": np.array([], dtype=np.int64),
                "year": np.array([], dtype=np.int64),
            }
        )
    if bad.size > 0:
        return pd.Series(bad["time"].values, index=bad[var]).sort_values()
    else:
        return pd.Series([], dtype="datetime64[ns]")
    # return bad[['time', 'month', 'year']]
