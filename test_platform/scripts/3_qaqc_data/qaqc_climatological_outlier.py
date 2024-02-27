"""
This is a script where Stage 3: QA/QC function(s) on climatological outlier values in the data observations are flagged. 
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

try:
    from qaqc_unusual_gaps import *
except:
    print("Error importing qaqc_unusual_gaps.py")

try:
    from qaqc_utils import *
except Exception as e:
    print("Error importing qaqc_utils: {}".format(e))
    
def open_log_file_clim(file):
    global log_file
    log_file = file
#----------------------------------------------------------------------
## climatological outlier check
def qaqc_climatological_outlier(df, winsorize=True, winz_limits=[0.05,0.05], plot=True, verbose=False, local=False):
    '''
    Flags individual gross outliers from climatological distribution.
    Only applied to air temperature and dew point temperature
    
    Input:
    ------
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        plots [bool]: if True, produces plots of any flagged data and saved to AWS
        winsorize [bool]: if True, raw observations are winsorized to remove spurious outliers first
        winz_limits [list]: if winsorize is True, values represent the low and high percentiles to standardize to
            
    Returns:
    --------
        qaqc success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        qaqc failure:
            None
            
    Flag meaning:
    -------------
        26,qaqc_climatological_outlier,Value flagged as a climatological outlier
    '''

    printf("Running: qaqc_climatological_outlier", log_file=log_file, verbose=verbose)
    
    vars_to_check = ['tas', 'tdps', 'tdps_derived']
    vars_to_anom = [v for v in vars_to_check if v in df.columns]

    try:
        # whole station bypass check
        df, pass_flag = qaqc_dist_whole_stn_bypass_check(df, vars_to_anom)
        
        if pass_flag == 'fail':
            return df
        else:
            for var in vars_to_anom:      
                # only work with non-flagged values
                printf('Checking for climatological outliers in: {}'.format(var), log_file=log_file, verbose=verbose)
                df_valid = df.loc[df[var+'_eraqc'].isnull() == True]

                # winsorize data by percentiles
                if winsorize == True:
                    df_std = winsorize_temps(df_valid, vars_to_anom, winz_limits)
                else:
                    df_std = df_valid

                # standardize data by monthly climatological anomalies by hour
                df_std = clim_standardized_anom(df_std, vars_to_anom)

                # apply low pass filter
                df_std = low_pass_filter(df_std, vars_to_anom)

                # gaussian is fitted to the histogram of anomalies for each month
                # similar to distributional gap check
                # FUTURE DEV (v2 of product):
                    # HadISD: obs that fall between critical threshold value and gap or critical threshold and end of distribution are tentatively flagged
                    # May be later reinstated on comparison with good data from neighboring stations

                for month in range(1,13):

                    df_m = df_std.loc[df_valid.time.dt.month == month]

                    # some months will be missing data
                    if len(df_m) == 0:
                        continue

                    # determine number of bins
                    bins = create_bins(df_m[var])

                    # pdf
                    mu = np.nanmean(df_m[var])
                    sigma = np.nanstd(df_m[var])
                    y, left_bnd, right_bnd = pdf_bounds_clim(df_m[var], mu, sigma, bins)

                    # identify gaps as below y=0.1 from histogram, not pdf
                    y_hist, bins = np.histogram(df_m[var], bins=bins, density=True)

                    # identify bin indices outside of thresholds and check if bin is above 0.1
                    bins_to_check = [i for i, n in enumerate(bins) if n <= left_bnd or n >= right_bnd][:-1] # remove last item due to # of bins exceeding hist by 1
                    if len(bins_to_check) != 0:
                        for b in bins_to_check:
                            if y_hist[b] > 0.1:
                                printf('Flagging outliers in {0}, month {1}'.format(var, month), log_file=log_file, verbose=verbose)
                                # list of index of full df to flag, not standardized df
                                idx_to_flag = [i for i in df_m.loc[(df_m[var] >= bins[b]) & (df_m[var] < bins[b+1])].index]
                                for i in idx_to_flag:
                                    df.loc[df.index == i, var+'_eraqc'] = 26 # see era_qaqc_flag_meanings.csv 

            if plot:
                for var in vars_to_anom:
                    for month in range(1,13):
                        if 26 in df[var+'_eraqc'].values: # only plot a figure if flag is present
                            clim_outlier_plot(df, var, month, network=df['station'].unique()[0].split('_')[0], local=local) 
        # Drop month,year vars used for calculations
        df = df.drop(columns=['month','year'])
        return df
    
    except Exception as e:
        printf("qaqc_climatological_outlier failed with Exception: {}".format(e), log_file=log_file, verbose=verbose)
        return None

#----------------------------------------------------------------------
def clim_mon_mean_hourly(df, var, month, hour):
    '''Calculate the monthly mean climatology for each of the day'''
    
    df_m_h = df.loc[(df.time.dt.month == month) & (df.time.dt.hour == hour)]
    clim_value = df_m_h[var].mean(numeric_only = True)
    
    # special handling if value is nan? 
    
    return clim_value

#----------------------------------------------------------------------
def iqr_range_monhour(df, var, month, hour):
    '''Calculates the monthly interquartile range per hour'''
    
    q1 = df.loc[(df.time.dt.month == month) & (df.time.dt.hour == hour)].quantile(0.25, numeric_only=True)
    q3 = df.loc[(df.time.dt.month == month) & (df.time.dt.hour == hour)].quantile(0.75, numeric_only=True)
    
    iqr_df = q3 - q1
    iqr_df_val = iqr_df[var]
    
    # iqr cannot be less than 1.5°C in order to preserve low variance stations
    if iqr_df_val < 1.5:
        iqr_df_val = 1.5
    else:
        iqr_df_val = iqr_df_val
            
    return iqr_df_val

#----------------------------------------------------------------------
def clim_standardized_anom(df, vars_to_anom):
    '''
    First anomalizes data by monthly climatology for each hour, then
    standardizes by the monthly climatological anomaly IQR for each hour
    '''
    
    df2 = df.copy()
    
    for var in vars_to_anom:
        for m in range(1,13,1):
            for h in range(0,24,1):
                # each hour in each month
                anom_value = clim_mon_mean_hourly(df, var, month=m, hour=h)
                iqr_value = iqr_range_monhour(df, var, month=m, hour=h)
                
                # locate obs within specific month/hour
                df_m_h = df.loc[(df.time.dt.month == m) & (df.time.dt.hour == h)]
                
                # calculate the monthly climatological anomaly by hour and standardize by iqr
                df2.loc[(df.time.dt.month == m) & 
                        (df.time.dt.hour == h), 
                        var] = (df_m_h[var] - anom_value) / iqr_value
                
    return df2

#----------------------------------------------------------------------
def winsorize_temps(df, vars_to_anom, winz_limits):
    '''
    Replaces potential spurious outliers by limiting the extreme values
    using the winz_limits set (default is 5% and 95% percentiles)
    '''
    
    df2 = df.copy()
    
    for var in vars_to_anom:
        for m in range(1,13,1):
            for h in range(0,24,1):
                if h not in df.loc[df.time.dt.hour == h]:
                    continue # some stations only report some hours
                else:
                    df_m_h = df.loc[(df.time.dt.month == m) & (df.time.dt.hour == h)]

                    # winsorize only vars in vars_to_anom
                    df_w = stats.mstats.winsorize(df_m_h[var], limits=winz_limits, nan_policy='omit')

                    df2.loc[(df.time.dt.month == m) & (df.time.dt.hour == h),
                           var] = df_w
                
    return df2

#----------------------------------------------------------------------
def median_yr_anom(df, var):
    '''Get median anomaly per year'''
    
    monthly_anoms = []
    
    # identify years in data
    years = df.time.dt.year.unique()
    
    for yr in years:
        df_yr = df.loc[df.time.dt.year == yr]

        ann_anom = df_yr[var].median()
        monthly_anoms.append(ann_anom)
        
    return monthly_anoms

#----------------------------------------------------------------------
def low_pass_filter_weights(median_anoms, month_low, month_high, filter_low, filter_high):
    '''Calculates weights for low pass filter'''
    
    filter_wgts = [1, 2, 3, 2, 1]
    
    if np.sum(filter_wgts[filter_low:filter_high] * 
              np.ceil(median_anoms[month_low:month_high] - 
                      np.floor(median_anoms[month_low:month_high]))) == 0:
        weight = 0
    
    else:
        weight = (
            np.sum(filter_wgts[filter_low:filter_high] * np.ceil(median_anoms[month_low:month_high])) / 
            np.sum(filter_wgts[filter_low:filter_high] * np.ceil(median_anoms[month_low:month_high] - 
                                                                 np.floor(median_anoms[month_low:month_high])))
        )
        
    return weight

#----------------------------------------------------------------------
def low_pass_filter(df, vars_to_anom):
    '''
    Low pass filtering on observations to remove any climate change signal 
    causing overzealous removal at ends of time series
    '''
    # identify years in data
    years = df.time.dt.year.unique()
    
    for var in vars_to_anom:
        
        median_anoms = median_yr_anom(df, var)
    
        for yr in range(len(years)):
            if yr == 0:
                month_low, month_high = 0, 3
                filter_low, filter_high = 2, 5
                
            elif yr == 1:
                month_low, month_high = 0, 4
                filter_low, filter_high = 1, 5
                
            elif yr == len(years)-2:
                month_low, month_high = -4, -1
                filter_low, filter_high = 0, 3

            elif yr == len(years)-1:
                month_low, month_high = -3, -1
                filter_low, filter_high = 0, 2

            else:
                month_low, month_high = yr-2, yr+3
                filter_low, filter_high = 0, 5
                            
            if np.sum(np.abs(median_anoms[month_low:month_high])) != 0:
                weights = low_pass_filter_weights(median_anoms, month_low, month_high, filter_low, filter_high)
                      
            # want to return specific year of data at a specific variable, the variable minus weight value
            df.loc[(df.time.dt.year == years[yr]), var] = df.loc[df.time.dt.year == years[yr]][var] - weights
            
    return df

#----------------------------------------------------------------------
def pdf_bounds_clim(df, mu, sigma, bins):
    '''Calculate pdf distribution, return pdf and threshold bounds'''

    y = stats.norm.pdf(bins, mu, sigma)
    
    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1)
    
    # find first index
    left_bnd = round(bins[pdf_bounds[0][0] -1])
    right_bnd = round(bins[pdf_bounds[-1][0] + 1])
    
    return (y, left_bnd - 1, right_bnd)