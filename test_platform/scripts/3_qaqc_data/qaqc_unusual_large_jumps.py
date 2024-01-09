"""
This is a script where Stage 3: QA/QC function(s) on unusual large jumps / spikes with data observations are flagged. 
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

## Import plotting functions
try:
    from qaqc_plot import *
except:
    print("Error importing qaqc_plot.py")

#-----------------------------------------------------------------------------
## unusual large jumps (spike) + helper functions
def qaqc_unusual_large_jumps(df, iqr_thresh=6, min_datapoints=50, plot=True, local=False, verbose=True):
    """
    Test for unusual large jumps or ''spikes'', given the statistics of the series. Analysis for each individual month in 
    time series to account for seasonal cycles in different regions.
    
    This test is done for ["tas", "tdps", "ps", "psl", "ps_altimeter"]
    Should it be done for more vars?
    
    Input:
    -----
        df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
        iqr_thresh [int] : critical value (iqr_thresh*IQR) for spike detection (default=6)
        min_datapoints [int] : minimum data points in each month to be valid for testing (default=50)
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
        23,qaqc_unusual_large_jumps,Unusual jump (spike) in variable

    Notes:
    ------
    To Do:
    - iqr_thresh is something can be modified of tweaked down the line (6 is what HadISD uses)
    - min_datapoints is the minimum data points in a group for threshold calculation (month/hours between data points)
    - HadISD uses 100, this can be modified and tweaked in future development
    """
    df = df.copy(deep=True)
    df.set_index(df['time'], inplace=True)
    df.drop(columns=['time'], inplace=True)

    try:
    # if True:

        # Drop station index
        # df = df.droplevel(level="station")

        # Define test variables and check if they are in the dataframe
        check_vars = ["tas", "tdps", "tdps_derived", 'ps', 'psl', 'ps_altimeter', 'ps_derived']
        variables = [var for var in check_vars if var in df.columns]

        if verbose:
            print("Running {} on {}".format("qaqc_unusual_large_jumps", variables))

        # Loop through test variables
        for var in variables:
            print('Running unusual large jumps check on: {}'.format(var))

            # Create a copy of the original dataframe and drop NaNs in the testing variable
            # new_df = df.copy(deep=True)
            # new_df = new_df.dropna(subset=var)#.drop(columns=["lat","lon","elevation"])

            # Use only values that have not been flagged by previous QAQC tests
            # valid = np.where(np.isnan(df[var+"_eraqc"]))[0]
            # new_df = df.iloc[valid]

            new_df = df.loc[df[var+'_eraqc'].isnull() == True]

           # first scans suspect values using entire record
            if new_df[var].isna().all() == True:
                continue # bypass to next variable if all obs are nans 
            
            # Detect spikes
            new_df = detect_spikes(new_df, var=var, iqr_thresh=iqr_thresh, min_datapoints=min_datapoints)

            # Retrieve location of spikes
            ind = new_df.index[np.where(new_df[var+"_spikes"])[0]]

            # Flag _eraqc variable
            for i in ind:
                df.loc[df.index == i, var+"_eraqc"] = 23 # see qaqc_flag_meanings.csv
            
            if plot:
                unusual_jumps_plot(df, var, flagval=23, local=local)
                for i in ind:
                    try:
                    # if True:
                        subset = df.loc[(df.index >= i - datetime.timedelta(hours=48)) & 
                                        (df.index <= i + datetime.timedelta(hours=48))]
                        # subset = np.logical_and(df.index >= i - np.timedelta64(48,'h'), 
                        #                     df.index <= i + np.timedelta64(48,'h'))
                        unusual_jumps_plot(subset, var, flagval=23, date=i, local=local)
                    except:
                    # else:
                        print('Unable to plot {0} detailed unusual jumps figure for {1}'.format(i, var))
                        continue

        return df.reset_index()
        return df

    except Exception as e:
        print("qaqc_unusual_large_jumps failed with Exception: {}".format(e))
        return None

#-----------------------------------------------------------------------------
def potential_spike_check(potential_spike, diff, crit, hours_diff):
    """
    Checks for neccessary conditions for a potential spike to be an actual spike:
     - Spikes are considered 1-value spike up to 3-values spike
     - Difference right before the spike should be lower than half the critical value
     - Difference at the actual spike must be higher than the critical value
     - Differences within the multi-value spike must lower than half the critical value
     - Difference right after (spike exit) the spike should be higher than the critical value and
       of opposite sign of the actual spike

    Input:
    -----
        potential_spike [pandas series] : bool pd.Series with True on potential spike location
        diff [pandas series] : float pd.Series with differences in the test variable
        crit [pandas series] : float pd.Series with the critical value for the differences in the test variable
        crit [pandas series] : float pd.Series with the hour differences between data points in the test variable

    Output:
    ------
        df [pandas dataframe] : input df with added `var`_spike column True where data matches the spike conditions
    """

    potential_spike = potential_spike.copy(deep=True)
    
    ind = np.where(potential_spike)[0]
    spikes = pd.Series(np.zeros_like(potential_spike).astype("bool"), index=potential_spike.index)
    dates = pd.Series(potential_spike.index.values)
    # print(len(ind))
    # print(dates[ind])
    
    for i in ind:
        
        #Ignore edges for now
        if i==1 or i>=len(potential_spike)-4:
            continue
        # Indices, critical values, and values before and after potential spike
        im1, i0, ip1, ip2, ip3, ip4 = [i-1, i, i+1, i+2, i+3, i+4]
        tm1, t0, tp1, tp2, tp3, tp4 = diff.iloc[[im1, i0, ip1, ip2, ip3, ip4]]
        cm1, c0, cp1, cp2, cp3, cp4 = crit.iloc[[im1, i0, ip1, ip2, ip3, ip4]]
        # print(dates.iloc[[i0,ip1]])
        # print(" sgn(t0)[{}] != sgn(tm1)[{}] \n abs(tm1)[{}] < cm1[{}] \n abs(tp1)[{}] > cp1[{}] \n abs(tp2)[{}] < cp2[{}]".format(np.sign(t0),np.sign(tp1), np.abs(tm1), cm1, np.abs(tp1), cp1, np.abs(tp2), cp2))
        # print()
        # Three-values spike
        if (
            np.sign(t0) != np.sign(tp2) and 
            np.abs(tm1) < 0.5*cm1 and 
            np.abs(tp1) < 0.5*cp1 and 
            np.abs(tp2) < 0.5*cp2 and 
            np.abs(tp3) > cp3 and 
            np.abs(tp4) < 0.5*cp4
        ):
            # print(dates.iloc[i0])
            spikes.iloc[[i0,ip1,ip2]] = True
            # i += 3
            # continue
            
        # Two-values spike
        elif (
            np.sign(t0) != np.sign(tp2) and 
            np.abs(tm1) < 0.5*cm1 and 
            np.abs(tp1) < 0.5*cp1 and 
            np.abs(tp2) > cp2 and 
            np.abs(tp3) < 0.5*cp3
        ):
            # print(dates.iloc[i0])
            spikes.iloc[[i0,ip1]] = True
            # i += 2
            # continue
        
        # One-value spike
        elif( 
            np.sign(t0) != np.sign(tp1) and 
            np.abs(tm1) < 1.0*cm1 and 
            np.abs(tp1) > cp1 and 
            np.abs(tp2) < 1.0*cp2
        ):
            # print(dates.iloc[i0])
            spikes.iloc[i0] = True
            # i += 1
            # continue
        
    return spikes

#-----------------------------------------------------------------------------
def detect_spikes(df, var, iqr_thresh=6, min_datapoints=50):
    """
    Detect  unusual large jumps or ''spikes'' in the time series for `var`:
      1- Find potential unusual large jumps or ''spikes'' by comparing the differences in `var` to each 
         month's critical value (crit = iqr_thresh * IQR)
      2- `potential_spike_check` checks for neccessary conditions for a potential spike to be an actual
          spike
    
    This test is done for ["tas", "tdps", "tdps_derived", 'ps', 'psl', 'ps_altimeter', 'ps_derived']
    Should it be done for more vars?
    
    Input:
    -----
        df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
        var [str] : variable to test
        iqr_thresh [int] : critical value (iqr_thresh*IQR) for spike detection (default=6)
        min_datapoints [int] : minimum data points in each month to be valid for testing (default=50)

    Output:
    ------
        df [pandas dataframe] : input df with added columns for spike check
                 
    Notes:
    ------
    To Do:
    - iqr_thresh is something can me modified of tweaked down the line (6 is what HadISD uses)
    - min_datapoints is the minimum data points in a group for threshold calculation (month/hours between data points)
    - HadISD uses 100, this can be modified and twaked in future development
    """
    
    # Make a copy of the original dataframe
    df = df.copy(deep=True)
    
    # Calculate difference in var values
    df[var+'_difference'] = df[var].diff().fillna(0)
    
    # Calculate dates
    df['date'] = df.index.values
    
    # Calculate time difference
    df['time_diff'] = df['date'].diff().fillna(pd.Timedelta(0))
    
    # Calculate time differece in hours
    df['hours_diff'] = df['time_diff']/np.timedelta64(1, 'h')
    df = df[np.logical_and(df['hours_diff'] > 0, df['hours_diff'] <= 12)]
    
    # Group by month to avoid strong seasonal cycle
    # grouped = df.groupby([pd.Grouper(freq='M'), df['hours_diff']]) 
    grouped = df.groupby(pd.Grouper(freq='M')) 

    # Count number of data per month
    counts = grouped[var+'_difference'].transform("count")
    df[var+'_counts'] = counts
    # Keep only months with more than 50 values to be statistically valid
    df = df[df[var+'_counts']>min_datapoints]

    # Define modified IQR 
    # kwargs = {'rng':(20, 80),}
    kwargs = {}
    
    # Calculate iqr
    iqr = grouped[var+'_difference'].transform(stats.iqr, **kwargs)
    df[var+'_iqr'] = iqr

    # Calculate critical value as rounded-up 6 (or defined by argument) times IQR
    df[var+'_critical'] = np.ceil(iqr_thresh*df[var+'_iqr'])
    
    # Find potential spike values where var diff is higher than the critical value
    df[var+'_potential_spikes'] = np.abs(df[var+'_difference'])>df[var+'_critical']

    # Filter real spikes using `potential_spike_check` function
    spikes = potential_spike_check(df[var+'_potential_spikes'], df[var+'_difference'], df[var+'_critical'], df['hours_diff'])
    df[var+'_spikes'] = spikes
    
    return df
