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
import scipy.signal as signal

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
######################################
# #FOR DEBUG
# #UNCOMMENT FOR NOTEBOOK DEBUGGING
# global log_file
# log_file = open("logtest.log","w")
# verbose=True
######################################

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

    new_df = df.copy()
    
    vars_to_check = ['tas', 'tdps', 'tdps_derived']
    # vars_to_check = ['tas']
    vars_to_anom = [v for v in vars_to_check if v in df.columns]

    try:
    # if True:
        printf("Running {} on {}".format("qaqc_climatological_outlier", vars_to_anom), verbose=verbose, log_file=log_file)

        # whole station bypass check first
        # df is already flagged by gaps function (careful if the order is modified)
        # df,stn_length = qaqc_dist_whole_stn_bypass_check(df, vars_to_check, min_num_months=iqr_thresh)
    
        for var in vars_to_anom:      
            # only work with non-flagged values
            printf('Checking for climatological outliers in: {}'.format(var), log_file=log_file, verbose=verbose)
            df_valid = grab_valid_obs(new_df, var, kind='drop') # subset for valid obs, distribution drop yellow flags
            # Keep only useful columns
            df_valid = df_valid[[var,"year","month","day","hour","time"]]
            
            # Bypass if there are not valid observations
            if df_valid[var].size ==0:
                printf('Not valid observations for: {}'.format(var), log_file=log_file, verbose=verbose)
                continue
            
            # Winsorize data and calculate climatology by month/hour with winsorized data
            clim = df_valid.groupby(["month","hour"])[var].transform(lambda row:stats.mstats.winsorize(row, limits=[0.05,0.05]))
            clim = pd.DataFrame(data={var:clim, "hour":df_valid.hour, "month":df_valid.month}, index=df_valid.index)
            clim = clim.groupby(["month","hour"])[var].transform(lambda row: np.nanmean(row))
            clim = pd.DataFrame(data={var:clim, "hour":df_valid.hour, "month":df_valid.month}, index=df_valid.index)

            # Anomalize using winsorized month/hour climatologies
            anom = df_valid[var]-clim[var]
            anom = pd.DataFrame(data={var:anom, "hour":df_valid.hour, "month":df_valid.month}, index=df_valid.index)

            # Calculate IQR by month/hour
            iqr = anom.groupby(["month","hour"])[var].transform(lambda row: max(np.nanpercentile(row, 75)-np.nanpercentile(row, 25), 1.5))
            iqr = pd.DataFrame(data={var:iqr, "hour":df_valid.hour, "month":df_valid.month}, index=df_valid.index)
            
            # Standardise anomalies by month/hour IQR
            std = anom[var]/iqr[var]
            std = pd.DataFrame(data={var:std, "hour":df_valid.hour, "month":df_valid.month}, index=df_valid.index)

            # Low-pass standardised data
            cut_freq = 1/(3600*24*365/30) # In Hz (cut_period : 1 month)
            data_freq = 1/(df_valid['time'].diff().mode().values[0].astype("float")/1e9) # In Hz
            sos = signal.butter(1, cut_freq, 'lp', output='sos', fs=data_freq)
            filtered = signal.sosfilt(sos, std[var].interpolate(method="linear"))
            # data = std[var].interpolate(method="linear")-filtered
            df_valid['raw_'+var] = df_valid[var] 
            df_valid[var] = filtered

            # Flag outliers
            printf('Flagging outliers in {0}'.format(var), log_file=log_file, verbose=verbose)
            df_valid['flag'] = df_valid.groupby(["month","hour"])[var].transform(lambda row: flag_clim_outliers(row))

            # Save original for plotting
            df_plot = df_valid.copy()

            # Drop all non-flagged values
            df_valid = df_valid.dropna(subset=['flag'])

            # # Debug
            # if "tdps" in var:
            #     display(df_valid)
        
            # Flag original data
            new_df.loc[new_df.time.isin(df_valid.time), var+'_eraqc'] = df_valid['flag']

            # Plot flagged values
            if plot:
                # Extract only flagged values to loop over those months and hours
                df_plot = df_plot[['year','hour','month','time','flag',var]].set_index(["month","hour"])

                # Loop over flagged months/hours
                for i in df_valid.set_index(["month","hour"]).index.unique():

                    # Extract actual month/hour from index
                    month, hour = i
                    # Extrac station name
                    station=df['station'].unique()[0]
                    # Plot distribution
                    clim_outlier_plot(df_plot.loc[i][var], month, hour, bin_size=0.1, station=station, local=True)
                        
        return new_df

    except Exception as e:
        printf("qaqc_climatological_outlier failed with Exception: {}".format(e), log_file=log_file, verbose=verbose)
        return None
        
#----------------------------------------------------------------------
def flag_clim_outliers(series):
    """
    """
    # Calculate frequency, normal fit, and boumdaries for clim outliers
    freq, bins, p, left, right = fit_normal(series, bin_size=0.10, plot=False)

    # Calculate bins for frequency checks
    freq_bins = np.concatenate((bins[1:int(len(bins)/2)], [0,0], bins[int(len(bins)/2)+1:-1]))

    # Flag series given the distribution thresholds (left and right)
    flag = gap_search(freq, left, right)
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Red left side of the distribution
    left_bad_bins = freq_bins[np.logical_and(flag==-1,freq_bins<0)]
    if len(left_bad_bins)>0:
        red_left = series <= left_bad_bins.max()
    else:
        red_left = np.zeros_like(series).astype("bool")
    
    # Red right side of the distribution
    right_bad_bins = freq_bins[np.logical_and(flag==-1,freq_bins>0)]
    if len(right_bad_bins)>0:
        red_right = series >= right_bad_bins.max()
    else:
        red_right = np.zeros_like(series).astype("bool")

    # Red flags 
    red = np.logical_or(red_left, red_right)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Yellow left side of the distribution
    left_probable_bins = freq_bins[np.logical_and(flag==0, freq_bins<0)]
    if len(left_probable_bins)>0:
        yellow_left = np.logical_and(series<=left_probable_bins.max(), ~red_left)
    else:
        yellow_left = np.zeros_like(series).astype("bool")
        
    # Yellow right side of the distribution
    right_probable_bins = freq_bins[np.logical_and(flag==0, freq_bins>0)]
    if len(right_probable_bins)>0:
        yellow_right = np.logical_and(series>=right_probable_bins.min(), ~red_right)
    else:
        yellow_right = np.zeros_like(series).astype("bool")
        
    # Yellow flags
    yellow = np.logical_or(yellow_left, yellow_right)

    # Create new array for the flags
    clim_outliers = np.ones_like(series)*np.nan
    
    #############################################################################
    # RED vs YELLOW flag with nearest neighbor stations for version 2
    # clim_outliers[red] = 30
    # clim_outliers[yellow] = 31
    
    # RED AND YELLOW FLAGGED
    clim_outliers[np.where(np.logical_or(red,yellow))[0]] = 26
    
    # ONLY RED FLAGGED
    # clim_outliers[np.where(red)[0]] = 26
    #############################################################################

    return clim_outliers

#----------------------------------------------------------------------
def fit_normal(series, bin_size=0.1, plot=False):
    """
    """
    bins = create_bins(series, bin_size=bin_size)
    max_bin = np.abs(bins).max()
    bins = np.arange(-max_bin-bin_size, max_bin+2*bin_size, bin_size)

    freq, bins = np.histogram(series, bins=bins, )
    area = sum(np.diff(bins) * freq)

    # Fit a normal distribution to the data
    mu, std = stats.norm.fit(series)
    p = stats.norm.pdf(bins, mu, std)*area

    try:
        left = np.where(np.logical_and(np.gradient(p)>0, p<=0.1))[0][-1]   # +1 # Manually shift the edge by one bin
    except:
        left = 1
    try:
        right = np.where(np.logical_and(np.gradient(p)<0, p<=0.1))[0][0]   # -1 # Manually shift the edge by one bin
    except:
        right = len(bins)-2
            
    if plot:
        # print(len(bins))
        # print(left, right)
        # Plot the histogram of the series
        fig,ax = plt.subplots()
        ax.hist(series, bins=bins, density=False, alpha=0.35, label='Histogram')
        ax.hist(series, bins=bins, density=False, lw=1.5, color="C0", histtype=u'step')
        ax.plot(bins, p, 'k', linewidth=2, label='Gaussian Fit')
        
        ax.set_yscale("log")
        ymin = min(0.08, freq[freq>0].min())
        # print(ymin, freq[freq>0].min())
        ymax = np.ceil(freq.max()/100)*100
        # ax.set_ylim(ymin,ymax) 
        # ax.set_xlim(xmin, xmax)
        ax.set_title('Histogram with Gaussian Fit')
        ax.set_xlabel('Value')
        ax.set_ylabel('Frequency')
        ax.legend()
        
        ax.axvline(bins[left],c="k",ls='--')
        ax.axvline(bins[right],c="k",ls='--')
        ax.axhline(0.1, c="k",ls=':')
        plt.show()
    return freq, bins, p, left, right
    
#----------------------------------------------------------------------
def gap_search(freq, left, right):

    """
    """
    left_freq = freq[0:left]
    left_flag = np.zeros_like(left_freq) # Yellow flag, all values beyond the threshold are flagged
    for i,f in zip(range(len(left_freq)-1,-1,-1), left_freq[::-1]):
        if f<0.1:
            left_flag[0:i+1] = -1  # Red flag, values and gap below 0.1 and beyond the threshold are flagged
            break

    right_freq = freq[right+1:]
    right_flag = np.zeros_like(right_freq) # Yellow flag, all values beyond the threshold are flagged
    for i,f in zip(range(len(right_freq)),right_freq):
        if f<0.1:
            right_flag[i:] = -1  # Red flag, values and gap below 0.1 and beyond the threshold are flagged
            break

    # Return flag of the size of freq
    flag = np.ones(len(freq)-len(left_freq)-len(right_freq))
    flag = np.concatenate((left_flag, flag, right_flag))
    return flag
    
#----------------------------------------------------------------------
def clim_mon_mean_hourly(df, var, month, hour):
    '''Calculate the monthly mean climatology for each of the day'''
    
    df_m_h = df.loc[(df['month'] == month) & (df['hour'] == hour)]
    clim_value = df_m_h[var].mean(numeric_only = True)
    
    # special handling if value is nan? 
    
    return clim_value

#----------------------------------------------------------------------
def iqr_range_monhour(df, var, month, hour):
    '''Calculates the monthly interquartile range per hour'''
    
    q1 = df.loc[(df['month'] == month) & (df['hour'] == hour)].quantile(0.25, numeric_only=True)
    q3 = df.loc[(df['month'] == month) & (df['hour'] == hour)].quantile(0.75, numeric_only=True)
    
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
                df_m_h = df.loc[(df['month'] == m) & (df['hour'] == h)]
                
                # calculate the monthly climatological anomaly by hour and standardize by iqr
                df2.loc[(df['month'] == m) & 
                        (df['hour'] == h), 
                        var] = (df_m_h[var] - anom_value) / iqr_value
                
    return df2

#----------------------------------------------------------------------
# def winsorize_temps(df, vars_to_anom, winz_limits):
def winsorize_temps(df, var, winz_limits):
    '''
    Replaces potential spurious outliers by limiting the extreme values
    using the winz_limits set (default is 5% and 95% percentiles)
    '''
    df2 = df.copy()
    df2.loc[:,[var]] = df2.groupby(["month","hour"])[var].transform(lambda row: stats.mstats.winsorize(row, limits=[0.05,0.05]))
    
    return df2
    
    # df2 = df.copy()
    
    # for var in vars_to_anom:
    #     for m in range(1,13,1):
    #         for h in range(0,24,1):
    #             if h not in df.loc[df['hour'] == h]:
    #                 continue # some stations only report some hours
    #             else:
    #                 df_m_h = df.loc[(df['month'] == m) & (df['hour'] == h)]

    #                 # winsorize only vars in vars_to_anom
    #                 df_w = stats.mstats.winsorize(df_m_h[var], limits=winz_limits, nan_policy='omit')

    #                 df2.loc[(df['month'] == m) & (df['hour'] == h),
    #                        var] = df_w
                
    # return df2

#----------------------------------------------------------------------
def median_yr_anom(df, var):
    '''Get median anomaly per year'''
    
    monthly_anoms = []
    
    # identify years in data
    years = df['year'].unique()
    
    for yr in years:
        df_yr = df.loc[df['year'] == yr]

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
def get_filter_indices(years):

    month_low, month_high, filter_low, filter_high = [],[],[],[]
    
    for yr in range(len(years)):
        if yr == 0:
            month_low.append(0)
            month_high.append(3)
            filter_low.append(2)
            filter_high.append(5)
        elif yr == 1:
            month_low.append(0)
            month_high.append(4)
            filter_low.append(1)
            filter_high.append(5)
        elif yr == len(years)-2:
            month_low.append(-4)
            month_high.append(-1)
            filter_low.append(0)
            filter_high.append(3)
        elif yr == len(years)-1:
            month_low.append(-3)
            month_high.append(-1)
            filter_low.append(0)
            filter_high.append(2)
        else:
            month_low.append(yr-2)
            month_high.append(yr+3)
            filter_low.append(0)
            filter_high.append(5)
    return month_low, month_high, filter_low, filter_high
    
#----------------------------------------------------------------------
# def low_pass_filter(df, vars_to_anom):
def low_pass_filter(df, var):
    '''
    Low pass filtering on observations to remove any climate change signal 
    causing overzealous removal at ends of time series
    '''

    median_anoms = df.groupby(["year"])[var].median().reset_index()
    month_low, month_high, filter_low, filter_high =  get_filter_indices(df['year'].unique())

    filter_weights = [low_pass_filter_weights(median_anoms[var], ml, mh, fl, fh) 
                      for ml, mh, fl, fh 
                      in zip(month_low, month_high, filter_low, filter_high)]
    median_anoms['filter_weights'] = filter_weights

    df = df.reset_index().\
            set_index("year").\
            assign(filter_weights=lambda x: median_anoms.set_index("year")['filter_weights'])
    df["year"] = df.index.values
    df = df.set_index("index")

    df.loc[:, [var]] = df[var] - df["filter_weights"] 
    
    # for yr in range(len(years)):
    
    # # identify years in data
    # years = df['year'].unique()
    
    # for var in vars_to_anom:
        
    #     median_anoms = median_yr_anom(df, var)
    
    #     for yr in range(len(years)):
    #         if yr == 0:
    #             month_low, month_high = 0, 3
    #             filter_low, filter_high = 2, 5
                
    #         elif yr == 1:
    #             month_low, month_high = 0, 4
    #             filter_low, filter_high = 1, 5
                
    #         elif yr == len(years)-2:
    #             month_low, month_high = -4, -1
    #             filter_low, filter_high = 0, 3

    #         elif yr == len(years)-1:
    #             month_low, month_high = -3, -1
    #             filter_low, filter_high = 0, 2

    #         else:
    #             month_low, month_high = yr-2, yr+3
    #             filter_low, filter_high = 0, 5
                            
    #         if np.sum(np.abs(median_anoms[month_low:month_high])) != 0:
    #             weights = low_pass_filter_weights(median_anoms, month_low, month_high, filter_low, filter_high)
                      
    #         # want to return specific year of data at a specific variable, the variable minus weight value
    #         df.loc[(df['year'] == years[yr]), var] = df.loc[df['year'] == years[yr]][var] - weights
            
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