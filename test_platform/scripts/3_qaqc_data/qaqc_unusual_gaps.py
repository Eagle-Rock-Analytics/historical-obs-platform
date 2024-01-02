"""
This is a script where Stage 3: QA/QC function(s) on unusual gaps / gaps within the monthly distribution with data observations are flagged. 
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
## distributional gap (unusual gap) + helper functions
def qaqc_unusual_gaps(df, iqr_thresh=5, plots=True):
    '''
    Runs all parts of the unusual gaps function, with a whole station bypass check first.

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        iqr_thresh [int]: interquartile range year threshold, default set to 5

    Output:
    -------
        if QAQC success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        if failure:
            None

    Flag meaning:
    -------------
        21,qaqc_unusual_gaps,Part 1: Monthly median value exceeds set threshold limits around monthly interquartile range for the monthly climatological median value
        22,qaqc_unusual_gaps,Part 2: Unusual gap in monthly distribution detected beyond PDF distribution
    
    Notes:
    ------
    PRELIMINARY: This function has not been fully evaluated or finalized in full qaqc process. Thresholds/decisions may change with refinement.
        - iqr_thresh preliminarily set to 5 years, pending revision
    '''

    # bypass check
    vars_to_remove = ['index','station','qc','duration','method',
                        'anemometer_height_m','thermometer_height_m',
                        'lat','lon','elevation','time','month','year',
                        'sfcWind_dir','hurs'] # list of var substrings to exclude if present in var
    vars_to_check = [var for var in df.columns if not any(True for item in vars_to_remove if item in var)] # remove all non-primary variables
    
    try:
        # whole station bypass check first
        df, pass_flag = qaqc_dist_whole_stn_bypass_check(df, vars_to_check)
        
        if pass_flag == 'fail':
            return df
        else:
            df_part1 = qaqc_dist_gap_part1(df, vars_to_check, iqr_thresh, plots)
            df_part2 = qaqc_dist_gap_part2(df_part1, vars_to_check, plots)

            if plots == True:
                for var in vars_to_check:
                    if (21 in df[var+'_eraqc'].values or 22 in df[var+'_eraqc'].values): # don't plot a figure if it's all nans/not enough months
                        flagged_timeseries_plot(df_part2, vars_to_check, flag_to_viz = [21, 22])
        
        return df_part2
    
    except Exception as e:
        print("qaqc_unusual_gaps failed with Exception: {}".format(e))
        return None


#-----------------------------------------------------------------------------
def qaqc_dist_whole_stn_bypass_check(df, vars_to_check, min_num_months=5):
    """
    Part 1: Checks the number of valid observation months in order to proceed through monthly distribution checks. 
    Identifies whether a station record has too few months and produces a fail pass flag. 

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        vars_to_check [list]: list of variables to run whole station bypass check on
        min_num_months [int]: minimum number of months required to pass check, default is 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        pass_flag [str]: pass flag indicating whether whole station pass/fail number of minimum months

    Flag meaning:
    -------------
        19,qaqc_unusual_gaps,Warning: Whole station has too few monthly observations to proceed through the monthly distribution gap check. Observations were not assessed for quality
    """

    # in order to grab the time information more easily -- would prefer not to do this
    # df['month'] = pd.to_datetime(df['time']).dt.month # sets month to new variable
    # df['year'] = pd.to_datetime(df['time']).dt.year # sets year to new variable
             
    # set up a "pass_flag" to determine if station proceeds through distribution function
    pass_flag = 'pass'
    
    for var in vars_to_check:
        for month in range(1,13):

            # first check num of months in order to continue
            month_to_check = df.loc[df['time'].dt.month == month]

            # check for number of obs years
            if (len(month_to_check.year.unique()) < 5):
                df[var+'_eraqc'] = 19 # see era_qaqc_flag_meanings.csv
                pass_flag = 'fail'

    err_statement = '{} has too short of an observation record to proceed through the monthly distribution qa/qc checks -- bypassing station'.format(df['station'].unique()[0])
    
    if pass_flag == 'fail':
        print(err_statement)
                
    return (df, pass_flag) 

#-----------------------------------------------------------------------------
def qaqc_dist_var_bypass_check(df, vars_to_check, min_num_months=5):
    """
    Part 1: Checks the number of valid observation months per variable to proceed through monthly distribution checks.
    Primarily assesses whether if null values persist for a month

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        vars_to_check [list]: list of variables to run whole station bypass check on
        min_num_months [int]: minimum number of months required to pass check, default is 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)

    Flag meaning:
    -------------
        20,qaqc_unusual_gaps,Warning: Variable has too few monthly observations to proceed through the monthly distribution gap check. Observations were not assessed for quality
    """
        
    for var in vars_to_check:
        for month in range(1,13):
            monthly_df = df.loc[df['time']dt.month == month]
            
            # if all values are null for that month across years
            if monthly_df[var].isnull().all() == True:
                df[var+'_eraqc'] = 20 # see era_qaqc_flag_meanings.csv
            
            # if not all months have nans, need to assess how many years do
            elif monthly_med(df).loc[monthly_med(df)['time']dt.month == month][var].isna().sum() > min_num_months:                
                df[var+'_eraqc'] = 20 # see era_qaqc_flag_meanings.csv
        
    return df

#-----------------------------------------------------------------------------
def qaqc_dist_gap_part1(df, vars_to_check, iqr_thresh, plot=True):
    """
    Part 1 / monthly check
        - compare anomalies of monthly median values
        - standardize against interquartile range
        - compare stepwise from the middle of the distribution outwards
        - asymmetries are identified and flagged if severe
    Goal: identifies suspect months and flags all obs within month

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        vars_to_check [list]: list of variables to run test on
        iqr_thresh [int]: interquartile range year threshold, default set to 5

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)

    Notes:
    ------
    PRELIMINARY: This function has not been fully evaluated or finalized in full qaqc process. Thresholds/decisions may change with refinement.
        - iqr_thresh preliminarily set to 5 years, pending revision
    """
        
    for var in vars_to_check:
        for month in range(1,13): 

            # per variable bypass check
            df = qaqc_dist_var_bypass_check(df, vars_to_check) # flag here is 20
            if 20 in df[var+'_eraqc']:
                continue # skip variable

            # station has above min_num_months number of valid observations, proceed with dist gap check
            else:
                # valid obs only
                df_valid = df.loc[df[var+'_eraqc'].isnull() == True]

                # calculate monthly climatological median, and bounds
                mid, low, high = standardized_median_bounds(df_valid, month, var, iqr_thresh=iqr_thresh)

                # calculate monthly median per month
                df_month = monthly_med(df_valid)

                for i in df_month.loc[df_month['time'].dt.month == month][var]:
                    if (i < low) or (i > high):
                        year_to_flag = (df_month.loc[(df_month[var]==i) & 
                                           (df_month['time'].dt.month == month)].dt.year.values[0])
                        print('Median {} value for {}-{} is beyond the {}*IQR limits -- flagging month'.format(
                            var,
                            month, 
                            int(year_to_flag),
                            iqr_thresh)
                        )

                        # flag all obs in that month
                        df.loc[(df_valid['time'].dt.month == month) & 
                               (df_valid['time'].dt.year == year_to_flag), var+'_eraqc'] = 21 # see era_qaqc_flag_meanings.csv

        if plot==True:
            for month in range(1,13):
                for var in vars_to_check:
                    if 21 in df[var+'_eraqc'].values: # don't plot a figure if nothing is flagged
                        dist_gap_part1_plot(df, month, var, flagval=21, iqr_thresh=iqr_thresh,
                                            network=df['station'].unique()[0].split('_')[0])
                
    return df

#-----------------------------------------------------------------------------
def qaqc_dist_gap_part2(df, vars_to_check, plot=True):
    """
    Part 2 / monthly check
        - compare all obs in a single month, all years
        - histogram created from all obs and gaussian distribution fitted
        - threshold values determined using positions where fitted freq falls below y=0.1
        - rounds outwards to next integer plus one
        - going outwards from center, distribution is scanned for gaps which occur outside threshold
        - obs beyond gap are flagged
    Goal: identifies individual suspect observations and flags the entire month 

    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        vars_to_check [list]: list of variables to run test on

    Output:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)

    Notes:
    ------
    PRELIMINARY: This function has not been fully evaluated or finalized in full qaqc process. Thresholds/decisions may change with refinement.
        - iqr_thresh preliminarily set to 5 years, pending revision 
    """

    # whole station bypass check first
    df, pass_flag = qaqc_dist_whole_stn_bypass_check(df, vars_to_check)
    
    if pass_flag != 'fail':
        
        for var in vars_to_check:
            for month in range(1,13):
                
                # per variable bypass check
                df = qaqc_dist_var_bypass_check(df, vars_to_check) # flag here is 20
                if 20 in df[var+'_eraqc']:
                    continue # skip variable 
                
                # station has above min_num_months number of valid observations, proceed with dist gap check
                else:
                    # valid obs only
                    df_valid = df.loc[df[var+'_eraqc'].isnull() == True]

                    # from center of distribution, scan for gaps (where bin = 0)
                    # when gap is found, and it is at least 2x bin width
                    # any bins beyond end of gap + beyond threshold value are flagged
                    
                    # subset by month
                    df_m = df_valid.loc[df_valid['time'].dt.month == month]
                    
                    # standardize against IQR range
                    df_month_iqr = standardized_iqr(df_m, var)

                    # determine number of bins
                    bins = create_bins(df_month_iqr)
                    
                    # pdf
                    mu = np.nanmean(df_month_iqr)
                    sigma = np.nanstd(df_month_iqr)

                    y, left_bnd, right_bnd = pdf_bounds(df_month_iqr, mu, sigma, bins)
                    
                    # identify gaps as below y=0.1 from histogram, not pdf                    
                    y_hist, bins = np.histogram(df_iqr, bins=bins, density=True)
                    
                    # identify climatology and iqr baselines in order to flag
                    iqr_baseline = iqr_range(df_valid, month=month, var=var)
                    clim = median_clim(df_valid, month=month, var=var)
                                        
                    # gaps are only flagged for values beyond left_bnd, right_bnd, as long as gap is 2*bin_width (2*0.25)
                    # considering that the # of bins for threshold is (4,7) from y=0.1
                    # safe to assume that gap is present if values >0.1 outside of left_bnd, right_bnd
                    bins_beyond_left_bnd = np.argwhere(bins <= left_bnd)
                    if len(bins_beyond_left_bnd) != 0: 
                        for data in bins_beyond_left_bnd:
                            if y_hist[data] > 0.1: # bins with data > 0.1 beyond left_bnd
                                
                                # identify values beyond left bnd
                                vals_to_flag = clim + (left_bnd * iqr_baseline) # left_bnd is negative
                                df.loc[df_valid[var] <= vals_to_flag[0], var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv

                    bins_beyond_right_bnd = np.argwhere(bins >= right_bnd)
                    if len(bins_beyond_right_bnd) != 0:
                        for data in bins_beyond_right_bnd:
                            if y_hist[data] > 0.1: # bins with data > 0.1 beyond right_bnd
                                
                                # identify values beyond right bnd
                                vals_to_flag = clim + (right_bnd * iqr_baseline) # upper limit threshold
                                df.loc[df_valid[var] >= vals_to_flag[0], var+'_eraqc'] = 22 # see era_qaqc_flag_meanings.csv
                    
    if plot==True:
        for month in range(1,13):
            for var in vars_to_check:
                if 20 not in df[var+'_eraqc'].values: # don't plot a figure if it's all nans/not enough months
                    if 22 in df[var+'_eraqc'].values: # don't plot a figure if nothing is flagged
                        dist_gap_part2_plot(df, month, var,
                                            network=df['station'].unique()[0].split('_')[0])
    
    return df  

#-----------------------------------------------------------------------------
def monthly_med(df):
    """Part 1: Calculates the monthly median"""
    return df.resample('M', on='time').median(numeric_only=True)

#-----------------------------------------------------------------------------
def iqr_range(df, month, var):
    """Part 1: Calculates the monthly interquartile range"""
    q1 = df.groupby('month').quantile(0.25, numeric_only=True)
    q3 = df.groupby('month').quantile(0.75, numeric_only=True)
    iqr_df = q3 - q1
    
    iqr_val = iqr_df.loc[iqr_df.index == month]
    
    # # inflated to 4°C or 4 hPa for months with very small IQR
    # var_check = ['tas', 'tdps', 'tdps_derived', 'ps', 'psl', 'psl_altimeter']
    # if iqr_val[var].values < 4:
    #     if var in var_check:
    #         iqr_val[var].values = 4
    
    return iqr_val[var].values

#-----------------------------------------------------------------------------
def standardized_iqr(df, var):
    """Part 2: Standardizes data against the interquartile range"""
    q1 = df[var].quantile(0.25)
    q3 = df[var].quantile(0.75)
    iqr = q3 - q1

    return (df[var].values - df[var].median()) / iqr
    
#-----------------------------------------------------------------------------
def median_clim(df, month, var):
    '''Part 2: Calculate climatological median for a specific month and variable'''
    clim = df[var].median(numeric_only=True)

    return clim

#-----------------------------------------------------------------------------
def standardized_anom(df, month, var):
    """
    Part 1: Calculates the monthly anomalies standardized by IQR range
    
    Output:
    -------
        arr_std_anom: array of monthly standardized anomalies for var
    """
    
    df_monthly_med = monthly_med(df)
    df_clim_med = median_clim(df)
    
    arr_anom = (df_monthly_med.loc[df_monthly_med['time'].dt.month == month][var].values -
                df_clim_med.loc[df_clim_med.index == month][var].values)
        
    arr_std_anom = arr_anom / iqr_range(df, month, var)
    
    return arr_std_anom
    
#-----------------------------------------------------------------------------
def standardized_median_bounds(df, month, var, iqr_thresh):
    """Part 1: Calculates the standardized median"""
    std_med = df.loc[df['time'].dt.month == month][var].median() # climatological median for that month
    
    lower_bnd = std_med - (iqr_thresh * iqr_range(df, month, var))
    upper_bnd = std_med + (iqr_thresh * iqr_range(df, month, var))
    
    return (std_med, lower_bnd[0], upper_bnd[0])

#-----------------------------------------------------------------------------
def create_bins(data, bin_size=0.25):
    '''Create bins from data covering entire data range'''

    # set up bins
    b_min = np.floor(np.nanmin(data))
    b_max = np.ceil(np.nanmax(data))
    # bins = np.arange(b_min - bin_size, b_max + (3. * bin_size), bin_size)
    bins = np.arange(b_min, b_max, bin_size)

    return bins

#-----------------------------------------------------------------------------
def pdf_bounds(df, mu, sigma, bins):
    '''Calculate pdf distribution, return pdf and threshold bounds'''

    y = stats.norm.pdf(bins, mu, sigma)
    
    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1)
    
    # find first index
    left_bnd = round(bins[pdf_bounds[0][0] -1])
    right_bnd = round(bins[pdf_bounds[-1][0] + 1])
    thresholds = (left_bnd - 1, right_bnd + 1)
    
    return (y, left_bnd - 1, right_bnd + 1)