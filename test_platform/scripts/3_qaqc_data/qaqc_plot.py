"""
This is a script where Stage 3: QA/QC related common plotting functions stored for ease of use
for the Historical Observations Platform.
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
    from qaqc_unusual_gaps import *
except:
    print("Error importing qaqc_unusual_gaps.py")

#============================================================================================================
# All plots helper plotting function for labeling, units, min, maxes
def _plot_format_helper(var):
    """
    Helper function for unusual large jumps plots
    
    Input:
    -----
        var [str] : variable name being plotted
    Ouput:
    ----- ylab, unit, mins[var]['North_America'], maxes[var]['North_America'])
        ylab [str] : variable name for y-axis
        unit [str] : units of variable
        miny [float] : min var value for y axis
        maxy [float] : max var value for y axis      
    """

    pr_vars = ['pr', 'pr_5min', 'pr_1h', 'pr_24h', 'pr_localmid']
    ps_vars = ['ps', 'psl', 'psl_altimeter', 'ps_derived']
    
    if var == 'tas':
        ylab = 'Air Temperature at 2m'
        unit = 'K'
        
    elif var == 'tdps' or var == 'tdps_derived':
        ylab = 'Dewpoint Temperature'
        unit = 'K'
        
    elif var == 'sfcWind':
        ylab = 'Surface Wind Speed'
        unit = '${m s^-1}$'
        
    elif var == 'sfcWind_dir':
        ylab = 'Surface Wind Direction'
        unit = 'degrees'
        
    elif var == 'rsds':
        ylab = 'Surface Radiation'
        unit = '${W m^-2}$'
        
    elif var == 'hurs':
        ylab = 'Humidity'
        unit = '%'
        
    elif var in pr_vars:
        ylab = 'Precipitation' # should be which precip var it is
        unit = 'mm'

    elif var in ps_vars:
        ylab = 'Pressure' # should eventually be what pressure var it is
        unit = 'Pa'
        
    T_X = {"North_America":329.92} #K
    T_N = {"North_America":210.15} #K
    D_X = {"North_America":329.85} #K
    D_N = {"North_America":173.15} #K
    W_X = {"North_America":113.2}  #m/s
    W_N = {"North_America":0.}     #m/s
    S_X = {"North_America":108330} #Pa
    S_N = {"North_America":87000}  #Pa

    maxes = {"tas": T_X, "tdps": D_X, "tdps_derived": D_X, "sfcWind": W_X, "psl": S_X, "ps": S_X}
    mins =  {"tas": T_N, "tdps": D_N, "tdps_derived": D_N, "sfcWind": W_N, "psl": S_N, "ps": S_N}
    miny = mins[var]['North_America']
    maxy = maxes[var]['North_America']
    
    return ylab, unit, miny, maxy

#============================================================================================================
## flagged timeseries plot
def flagged_timeseries_plot(df, vars_to_check, flag_to_viz):
    '''Produces a scatterplot timeseries figure of variables that have flags placed'''

    # can pass a list of flags
    for flag in flag_to_viz:
    
        # assess where each variable has flagged values
        for var in vars_to_check:
            flagged_data = df.loc[df[var+'_eraqc'] == flag]

            # only produce a plot if there is flagged values
            if len(flagged_data) == 0:
                continue

            # plot
            ax = df.plot.scatter(x='time', y=var, color='k', s=0.8, label='Valid')

            # plot flagged data
            flagged_data.plot.scatter(ax=ax, x='time', y=var, color='r', s=0.9, label='Flag: {}'.format(flag))

            # plot aesthetics
            plt.legend(loc='best', ncol=2)
            ylab, units, miny, maxy = _plot_format_helper(var)
            plt.ylabel('{} [{}]'.format(ylab, units));
            plt.xlabel('')
            plt.title('{0}'.format(df['station'].unique()[0]), fontsize=10);

            # save to AWS
            bucket_name = 'wecc-historical-wx'
            directory = '3_qaqc_wx'
            img_data = BytesIO()
            plt.savefig(img_data, format='png')
            img_data.seek(0)

            s3 = boto3.resource('s3')
            bucket = s3.Bucket(bucket_name)
            network = df['station'].unique()[0].split('_')[0]
            figname = 'flagged_timeseries_{0}_{1}'.format(df['station'].unique()[0], var)
            bucket.put_object(Body=img_data, ContentType='image/png',
                        Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                        directory, network, figname))

            # close figure for memory
            plt.close()

#============================================================================================================
## frequent values plotting functions
def frequent_plot_helper(df, var, bins, flag, yr, rad_scheme):
    '''Plotting helper with common plotting elements for all 3 versions of this plot'''
    
    # plot all valid data within year/season
    _plot = df.plot.hist(column=var, bins=bins, color='k', legend=False, alpha=0.5)
    
    # plot flagged values
    # first identify which values are flagged
    vals_to_flag = df.loc[df[var+'_eraqc'] == flag][var].unique()
        
    bars_to_flag = []
    for i in vals_to_flag:
        if math.isnan(i) == False:
            bars_to_flag.append(math.floor(i))
            
    # flag bars if too frequent
    for bar in _plot.patches:
        x = bar.get_x()
        if x in bars_to_flag: # right tail
            bar.set_color('r')
            
    # plot aesthetics
    xlab, units, miny, maxy = _plot_format_helper(var)
    plt.xlabel('{0} [{1}]'.format(xlab, units))
    yr_formatted = str(yr).replace('_', ' ') # simple formatting for plot aesthetic
    plt.annotate(yr_formatted, xy=(0.02, 0.95), xycoords='axes fraction', fontsize=10);
    plt.title('Frequent value check: {}'.format(df['station'].unique()[0]),
             fontsize=10);
    plt.legend(('Valid', 'Flagged'), loc='upper right')
    ax = plt.gca()
    leg = ax.get_legend()
    leg.legend_handles[0].set_color('k') # set valid to blue
    leg.legend_handles[-1].set_color('r') # set flagged bar to red
    
    if var == 'rsds':
        plt.annotate('Sfc. radiation option: \n{}'.format(rad_scheme), xy=(0.02, 0.85), xycoords='axes fraction', fontsize=10)
        
    # save figure to AWS
    network = df['station'].unique()[0].split('_')[0]
    
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    img_data = BytesIO()
    plt.savefig(img_data, format='png')
    img_data.seek(0)

    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_frequent_value_check_{0}_{1}_{2}'.format(df['station'].unique()[0], var, yr)
    bucket.put_object(Body=img_data, ContentType='image/png',
                     Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                     directory, network, figname))

    # close figure for memory
    plt.close()

#-----------------------------------------------------------------------------------------
def create_bins_frequent(df, var, bin_size=0.25):
    '''Create bins from data covering entire data range'''

    # ideally shouldn't have to have a separate function for this
    
    # grab data
    data = df[var]
    
    # radiation handling
    if var != 'rsds':
        # set up bins
        b_min = np.floor(np.nanmin(data))
        b_max = np.ceil(np.nanmax(data))
        bins = np.arange(b_min, b_max + 1, bin_size)
        
    else: 
        b_min = 0
        b_max = int(np.ceil(np.nanmax(data / 100))) * 100 # rounds up to next hundred
        bins = np.arange(b_min, b_max + 50, bin_size)
   
    return bins

#-----------------------------------------------------------------------------------------
def frequent_vals_plot(df, var, rad_scheme):
    '''
    Produces a histogram of the diagnostic histogram per variable, 
    and any bin that is indicated as "too frequent" by the qaqc_frequent_vals test 
    is visually flagged
    ''' 
    # bin sizes: using 1 degC for tas/tdps, and 1 hPa for ps vars
    ps_vars = ['ps', 'ps_altimeter', 'psl']
        
    if var in ps_vars: 
        bin_s = 100 # all of our pressure vars are in Pa, convert to 100 Pa bin size
    elif var == 'rsds':
        bin_s = 50
    else:
        bin_s = 1 
        
    bins = create_bins_frequent(df, var, bin_s)
    
    # first identify which values are flagged and "where"
    
    ## Year-by-year flag (23): plot all data for that year
    flag_df = df.loc[df[var+'_eraqc'] == 23]
    
    if len(flag_df) != 0:
        
        # identify year(s) with flagged data
        plot_yrs = flag_df['year'].unique()
        
        for y in plot_yrs:
            df_to_plot = df.loc[df['year']==y]
            _plot = frequent_plot_helper(df_to_plot, var, bins, flag=23, yr=y, rad_scheme=rad_scheme)
            
    ## Seasonal flag (24): plot all data for that year and season + specific handling for winter
    flag_df = df.loc[df[var+'_eraqc'] == 24]
    
    if len(flag_df) != 0:
        
        # identify unique years with flagged seasonal data
        plot_yrs = flag_df['year'].unique()
        
        for y in plot_yrs:
            df_year = df.loc[df['year']==y] # grab the entire year
            
            flagged_szns = df_year.loc[df_year[var+'_eraqc'] == 24]['month'].unique() # identify flagged months in that year
            
            if 3 in flagged_szns or 4 in flagged_szns or 5 in flagged_szns: # Spring - MAM
                df_to_plot = df_year.loc[(df_year['month']==3) | (df_year['month']==4) | (df_year['month']==5)]
                _plot = frequent_plot_helper(df_to_plot, var, bins, flag=24, yr=str(y)+'_spring', rad_scheme=rad_scheme)
                
            if 6 in flagged_szns or 7 in flagged_szns or 8 in flagged_szns: # Summer - JJA
                df_to_plot = df_year.loc[(df_year['month']==6) | (df_year['month']==7) | (df_year['month']==8)]
                _plot = frequent_plot_helper(df_to_plot, var, bins, flag=24, yr=str(y)+'_summer', rad_scheme=rad_scheme)

            if 9 in flagged_szns or 10 in flagged_szns or 11 in flagged_szns: # Autumn - SON
                df_to_plot = df_year.loc[(df_year['month']==9) | (df_year['month']==10) | (df_year['month']==11)]
                _plot = frequent_plot_helper(df_to_plot, var, bins, flag=24, yr=str(y)+'_autumn', rad_scheme=rad_scheme)
           
            if 12 in flagged_szns: # Winter - current year D + next year JF
                # special handling as follows
                # if the next year has flagged jan/feb, this will overwrite, but will be identical figure
                # some years will not have current year december and next year jan/feb so need this edge case                
                df_d = df_year.loc[df_year['month']==12] # current year dec
                df_jf = df.loc[(df['year']==y+1) & ((df['month']==1) | (df['month']==2))] # next year jan+feb
                df_to_plot = pd.concat([df_d, df_jf])
                _plot = frequent_plot_helper(df_to_plot, var, bins, flag=24, yr=str(y+1)+'_winter', rad_scheme=rad_scheme)
                
            if 1 in flagged_szns or 2 in flagged_szns: # Winter - previous year D + current year JF
                # special handling as follows
                # if the previous year has flagged december, this will overwrite, but will be identical figure
                # some years will not have previous year december and current jan/feb so need this edge case
                df_d = df.loc[(df['year']==y-1) & (df['month']==12)] # previous year dec
                df_jf = df_year[(df_year['month']==1) | (df_year['month']==2)] # current year jan+feb
                df_to_plot = pd.concat([df_d, df_jf])
                _plot = frequent_plot_helper(df_to_plot, var, bins, flag=24, yr=str(y)+'_winter', rad_scheme=rad_scheme)

#============================================================================================================
## distribution gap plotting functions
def dist_gap_part1_plot(df, month, var, flagval, iqr_thresh, network):
    '''
    Produces a timeseries plots of specific months and variables for part 1 of the unusual gaps function.
    Any variable that is flagged is noted
    '''

    # grab data by months
    df = df.loc[df['month'] == month]
        
    # grab flagged data
    flag_vals = df.loc[df[var + '_eraqc'] == flagval]
    
    # plot valid data
    ax = df.plot.scatter(x='time', y=var, label='Pass')
    
    # plot flagged data
    flag_vals.plot.scatter(ax=ax, x='time', y=var, color='r', label='Flagged')
    # should be consistent with other plots - I like Hector's open circles around flagged values

    # plot climatological median and threshold * IQR range
    mid, low_bnd, high_bnd = standardized_median_bounds(df, month, var, iqr_thresh=5)
    
    plt.axhline(y=mid, color='k', lw=0.5, label='Climatological monthly median')
    plt.fill_between(x=df['time'],
                    y1=low_bnd,
                    y2=high_bnd,
                    alpha=0.25, color='0.75', 
                    label='{} * IQR range'.format(iqr_thresh))
    
    # plot aesthetics
    plt.legend(loc='best')
    ylab, units, miny, maxy = _plot_format_helper(var)
    plt.ylabel('{} [{}]'.format(ylab, units));
    plt.xlabel('')
    plt.title('Distribution gap check pt 1: {0} / month: {1}'.format(
        df['station'].unique()[0],
        month), 
              fontsize=10);
    
    # save to AWS    
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_dist_gap_check_part1_{0}_{1}_{2}'.format(df['station'].unique()[0], var, month)
    bucket.put_object(Body=img_data, ContentType='image/png',
                 Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                 directory, network, figname))
    
    # close figures to save memory
    plt.close()

#-----------------------------------------------------------------------------------------
def dist_gap_part2_plot(df, month, var, network):
    '''
    Produces a histogram of the monthly standardized distribution
    with PDF overlay and threshold lines where pdf falls below y=0.1.
    Any bin that is outside of the threshold is visually flagged
    ''' 

    # select month
    df = df.loc[df['month'] == month]
    
    # standardize against IQR range
    df_month_iqr = standardized_iqr(df, var)
    
    # determine number of bins
    bins = create_bins(df_month_iqr)
    
    # plot histogram
    ax = plt.hist(df_month_iqr, bins=bins, log=False, density=True, alpha=0.3);
    xmin, xmax = plt.xlim()
    plt.ylim(ymin=0.1)

    # plot pdf
    mu = np.nanmean(df_month_iqr)
    sigma = np.nanstd(df_month_iqr)
    y = stats.norm.pdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'k--', linewidth=1)
    
    # add vertical lines to indicate thresholds where pdf y=0.1
    pdf_bounds = np.argwhere(y > 0.1)

    # find first index
    left_bnd = round(bins[pdf_bounds[0][0] -1])
    right_bnd = round(bins[pdf_bounds[-1][0] + 1])
    thresholds = (left_bnd - 1, right_bnd + 1)

    plt.axvline(thresholds[1], color='r') # right tail
    plt.axvline(thresholds[0], color='r') # left tail
    
    # flag (visually) obs that are beyond threshold
    for bar in ax[2].patches:
        x = bar.get_x() + 0.5 * bar.get_width()
        if x > thresholds[1]: # right tail
            bar.set_color('r')
        elif x < thresholds[0]: # left tail
            bar.set_color('r')

    # title and useful annotations
    plt.title('Distribution gap check, {0}: {1}'.format(df['station'].unique()[0], var), fontsize=10);
    plt.annotate('Month: {}'.format(month), xy=(0.025, 0.95), xycoords='axes fraction', fontsize=8);
    plt.annotate('Mean: {}'.format(round(mu,3)), xy=(0.025, 0.9), xycoords='axes fraction', fontsize=8);
    plt.annotate('Std.Dev: {}'.format(round(sigma,3)), xy=(0.025, 0.85), xycoords='axes fraction', fontsize=8);
    plt.ylabel('Frequency (obs)')
    
    # save figure to AWS
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    img_data = BytesIO()
    plt.savefig(img_data, format='png')
    img_data.seek(0)
    
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_dist_gap_check_part2_{0}_{1}_{2}'.format(df['station'].unique()[0], var, month)
    bucket.put_object(Body=img_data, ContentType='image/png',
                     Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                     directory, network, figname))

    # close figures to save memory
    plt.close()

#============================================================================================================
## unusual large jumps plotting functions
def unusual_jumps_plot(df, var, flagval=22, dpi=None, local=False, date=None):
    """
    Plots unusual large jumps qaqc result and uploads it to AWS (if local, also writes to local folder)
    Input:
    -----
        df [pd.Dataframe] : station pd.DataFrame from qaqc pipeline
        var [str] : variable name
        flagval [int] : flag value to plot (22 for unusual large jumps)
        dpi [int] : resolution for png plots
        local [bool] : if True, saves plot locally, else: only saves plot to AWS
    Ouput:
    ----- 
        None
    """
    
    # grab flagged data
    flag_vals = df.loc[df[var + '_eraqc'] == flagval]   
    
    # Create figure
    if date is not None:
        fig,ax = plt.subplots(figsize=(7,3))
    else:
        fig,ax = plt.subplots(figsize=(10,3))

    # Plot variable and flagged data
    df[var].plot(ax=ax, marker=".", ms=4, lw=1, color="k", alpha=0.5, label="Original data")
    
    flag_label = "{:.4f}% of data flagged".format(100*len(df.loc[df[var+"_eraqc"]==22, var])/len(df))
    df.loc[df[var+"_eraqc"]==22, var].plot(ax=ax, marker="o", ms=7, lw=0, mfc="none", color="C3", label=flag_label)    
    
    legend = ax.legend(loc=0, prop={'size': 8})    
        
    station = df['station'].unique()[0]
    network = station.split('_')[0]
    
    # Plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper(var)
    ylab = '{0} [{1}]'.format(ylab, units)
    
    ax.set_ylabel(ylab)
    ax.set_xlabel('')
    
    # We can set ylim since this function is supposed to be run after other QAQC functions (including world records)
    if date is not None:
        timestamp = str(date).split(":")[0].replace(" ","T")
    else:
        timestamp = "full_series"
        miny = max(miny, df[var].min())
        maxy = min(maxy, df[var].max())
        ax.set_ylim(miny,maxy)
    
    title = 'Unusual large jumps check: {0}'.format(station)
    ax.set_title(title, fontsize=10)
    
    # save to AWS
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    figname = 'qaqc_figs/qaqc_unusual_large_jumps_{0}_{1}_{2}'.format(station, var, timestamp)
    key = '{0}/{1}/{2}.png'.format(directory, network, figname)
    img_data = BytesIO()
    fig.savefig(img_data, format='png', dpi=dpi, bbox_inches="tight")
    img_data.seek(0)
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    bucket.put_object(Body=img_data, ContentType='image/png', Key=key)
    plt.close()
    if local:
        fig.savefig(figname+".png", format='png', dpi=dpi, bbox_inches="tight")
    
    return 

#============================================================================================================
def clim_outlier_plot(df, var, month, network):
    '''
    Produces a histogram of monthly standardized distribution
    with PDF overlay and threshold lines where pdf falls below y=0.1.
    Any bin that is outside of the threshold is visually flagged.
    
    Differs from dist_gap_part2_plot for the climatological outlier
    as IQR standardization does not occur within plotting
    '''
    
    # select month
    df_to_plot = df.loc[df.time.dt.month == month][var]
    
    # determine number of bins
    bins = create_bins(df_to_plot)
    
    # plot histogram
    ax = plt.hist(df_to_plot, bins=bins, log=False, density=True, alpha=0.3)
    xmin, xmax = plt.xlim()
    plt.ylim(ymin=0.1)
    
    # # plot pdf
    # mu = np.nanmean(df_to_plot)
    # sigma = np.nanmean(df_to_plot)
    # y = stats.norm.pdf(bins, mu, sigma)
    # l = plt.plot(bins, y, 'k--', linewidth=1)
    
    # # add vertical lines to indicate thresholds where pdf y=0.1
    # pdf_bounds = np.argwhere(y > 0.1)
    
    # # find first index
    # left_bnd = round(bins[pdf_bounds[0][0] - 1])
    # right_bnd = round(bins[pdf_bounds[-1][0] + 1])
    # thresholds = (left_bnd, right_bnd)
    # print(thresholds)
    
    # plt.axvline(thresholds[1], color='r') # right tail
    # plt.axvline(thresholds[0], color='r') # left tail
    
    # # flag visually obs that are beyond threshold
    # for bar in ax[2].patches:
    #     x = bar.get_x() + 0.5 * bar.get_width()
    #     if x > thresholds[1]: # right tail
    #         bar.set_color('r')
    #     elif x < thresholds[0]: # left tail
    #         bar.set_color('r')
            
    # title and useful annotations
    plt.title('Climatological outlier check, {0}: {1}'.format(df['station'].unique()[0], var), fontsize=10);
    plt.annotate('Month: {}'.format(month), xy=(0.025, 0.95), xycoords='axes fraction', fontsize=8);
    # plt.annotate('Mean: {}'.format(round(mu,3)), xy=(0.025, 0.9), xycoords='axes fraction', fontsize=8);
    # plt.annotate('Std.Dev: {}'.format(round(sigma,3)), xy=(0.025, 0.85), xycoords='axes fraction', fontsize=8);
    plt.ylabel('Frequency (obs)')
    
    # save figure to AWS
    bucket_name = 'wecc-historical-wx'
    directory = '3_qaqc_wx'
    img_data = BytesIO()
    plt.savefig(img_data, format='png')
    img_data.seek(0)
    
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket_name)
    figname = 'qaqc_climatological_outlier_{0}_{1}_{2}'.format(df['station'].unique()[0], var, month)
    bucket.put_object(Body=img_data, ContentType='image/png',
                     Key='{0}/{1}/qaqc_figs/{2}.png'.format(
                     directory, network, figname))
    
    # close figures to save memory
    plt.close()

    #============================================================================================================