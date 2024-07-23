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

## Import plotting functions
try:
    from qaqc_plot import *
except:
    print("Error importing qaqc_plot.py")

try:
    from qaqc_utils import *
except Exception as e:
    print("Error importing qaqc_utils: {}".format(e))
    
def open_log_file_frequent(file):
    global log_file
    log_file = file

# #####################################
# #FOR DEBUG
# #UNCOMMENT FOR NOTEBOOK DEBUGGING
# global log_file
# log_file = open("logtest.log","w")
# verbose=True
# #####################################
    
## frequent values + helper functions
#-----------------------------------------------------------------------------
def qaqc_frequent_vals(df, rad_scheme, plots=True, verbose=False, local=False):
    '''
    Test for unusually frequent values. This check is performed in two phases.
    Phase 1: Check is applied to all observations for a designated variable. If the current bin has >50% + >30 number of observations
    compared to +/- 3 surrounding bins, the current bin is highlighted for further check on the year-by-year basis. If the bin persists 
    as unusually frequent, the bin is flagged.
    Phase 2: Check is applied on a seasonal basis, for all observations within that season (mirroring phase 1). If a suspect bin is noted
    in the all observations stage, the check is performed on the year-by-year basis for that season. 

    This test is synergistically applied for air temperature and dew point temperature.         

    Input:
    -----
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        rad_scheme [str]: radiation handling for frequent occurence of valid zeros
        plots [bool]: if True, produces plots of any flagged data and saved to AWS

    Returns:
    -------
        qaqc success:
            df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
        qaqc failure:
            None
    
    Flag meaning:
    -------------
        24,qaqc_frequent_vals,Value flagged as unusually frequent in occurrence at the annual scale after assessing the entire observation record. Temperature and dew point temperature are synergistically flagged.
        25,qaqc_frequent_vals,Value flagged as unusually frequent in occurrence at the seasonal scale after assessing the entire observation record. Temperature and dew point temperature are synergistically flagged.
    '''
    # import pdb; pdb.set_trace()
    printf("Running: qaqc_frequent_vals", log_file=log_file, verbose=verbose)
    
    # this check is only done on air temp, dewpoint temp, and pressure
    vars_to_remove = ['qc', 'duration', 'method', 'flag', 'depth'] # list of var substrings to remove if present in var
    vars_to_include = ['tas', 'tdps', 'ps', 'psl', 'ps_altimeter', 'ps_derived', 'rsds'] 
    vars_to_check = [var for var in df.columns if any(True for item in vars_to_include if item in var) and not any(True for item in vars_to_remove if item in var)]

    try:
        printf("Running qaqc_frequent_vals on {}".format(vars_to_check), log_file=log_file, verbose=verbose)

        for var in vars_to_check:
            # if var=="rsds":
            #     import pdb; pdb.set_trace()
            printf('Running frequent values check on: {}'.format(var), log_file=log_file, verbose=verbose)
            df_valid = grab_valid_obs(df, var) # subset for valid obs
            
            # first scans suspect values using entire record -- this is the per variable check?
            # all years
            if df_valid[var].isna().all() == True:
                continue # bypass to next variable if all obs are nans 

            df_valid = frequent_bincheck(df_valid, var, data_group='all', rad_scheme=rad_scheme, verbose=verbose)

            # if no values are flagged as suspect, end function, no need to proceed
            if len(df_valid.loc[df_valid[var+'_eraqc'] == 100]) == 0:
                printf('No unusually frequent values detected for entire {} observation record'.format(var), log_file=log_file, verbose=verbose)
                # goes to seasonal check, no bypass

            else:
                # year by year
                # then scans for each value on a year-by-year basis to flag if they are a problem within that year
                    # DECISION: the annual check uses the unfiltered data
                    # previously flagged values are included here -- this would interfere with our entire workflow
                df_valid = frequent_bincheck(df_valid, var, data_group='annual', rad_scheme=rad_scheme, verbose=verbose)

            # seasonal scan (JF+D, MAM, JJA, SON) 
            # each season is scanned over entire record to identify problem values
            # only flags applied on annual basis using the three months on their own
            # NOTE: HadISD approach is to use the current year's december, rather than the preceeding december

            # seasonal version because seasonal shift in distribution of temps/dewpoints can reveal hidden values
            # all years
            df_valid = frequent_bincheck(df_valid, var, data_group='seasonal_all', rad_scheme=rad_scheme, verbose=verbose) ## DECISION: December is from the current year
            if len(df_valid.loc[df_valid[var+'_eraqc'] == 100]) == 0:
                printf('No unusually frequent values detected for seasonal {} observation record'.format(var), log_file=log_file, verbose=verbose)
                continue # bypasses to next variable

            else:
                printf('Unusually frequent values detected in seasonal distribution, continuing to annual check', log_file=log_file, verbose=verbose)
                # year by year --> December selection must be specific
                df_valid = frequent_bincheck(df_valid, var, data_group='seasonal_annual', rad_scheme=rad_scheme, verbose=verbose)    
                        
            # remove any lingering preliminary flags, data passed check
            df_valid.loc[df_valid[var+'_eraqc'] == 100, var+'_eraqc'] = np.nan

            # apply unique df_valid flags into full df
            isFlagged = df_valid.loc[df_valid[var+'_eraqc'].isnull() == False]
            for i in isFlagged.index:
                flag_to_place = isFlagged.loc[isFlagged.index == i][var+'_eraqc'].values[0]
                df.loc[isFlagged.index, var+'_eraqc'] = flag_to_place

        # synergistic flag on tas and tdps/tdps_derived
        # first establish at least tas and one tdps var present
        temp_vars = ['tas', 'tdps', 'tdps_derived']
        num_temp_vars = [var for var in vars_to_check if var in temp_vars]
        if len(num_temp_vars) != 1 and 'tas' in num_temp_vars:
            # proceed to synergistic check
            df = synergistic_flag(df, num_temp_vars)

        # plots item
        if plots:
            for var in vars_to_check:
                if 24 in df[var+'_eraqc'].unique() or 25 in df[var+'_eraqc'].unique(): # only plot a figure if a value is flagged
                    frequent_vals_plot(df, var, rad_scheme, local=local)

        # Drop month,year vars used for calculations
        # df = df.drop(columns=['month','year'])
        return df
    
    except Exception as e:
        printf("qaqc_frequent_vals failed with Exception: {}".format(e), log_file=log_file, verbose=verbose)
        return None

#-----------------------------------------------------------------------------
def frequent_bincheck(df, var, data_group, rad_scheme, verbose=False):
    '''
    Approach: 
        - histograms created with 0.5 or 1.0 or hpa increments (depending on accuracy of instrument)
        - each bin compared to the three on either side
        - if this bin contains more than half the total population of the seven bins combined
        - and more than 30 observations over the station record (20 for seasonal)
        - then histogram bin is highlighted for further investigation
        - minimum number limit imposted to avoid removing true tails of distribution

    Input:
    -----
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        var [str]: variable to run check on
        data_group [str]: annual vs. seasonal handling, options: all, annual, seasonal_all, seasonal_annual
        rad_scheme [str]: radiation handling for frequent occurence of valid zeros

    Returns:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
    '''    
    # if var=="rsds":
    #     import pdb; pdb.set_trace() 
    # seasons
    szns = [[3,4,5], [6,7,8], [9,10,11], [12,1,2]] 
    
    # bin sizes: using 1 degC for tas/tdps, and 1 hPa for ps vars
    ps_vars = ['ps', 'ps_altimeter', 'psl', 'ps_derived']
    
    if var in ps_vars: 
        bin_s = 100 # all of our pressure vars are in Pa, convert to 100 Pa bin size
    elif var == 'rsds':
        bin_s = 50
    else:
        bin_s = 1 
         
    # radiation schemes for assessment
    if var == 'rsds':
        if rad_scheme == 'all_hours':
            # all valid observations included -- frequent flag will set on 0/nighttime hours
            printf('Radiation frequent value check scheme: all_hours selected, will flag nighttime', log_file=log_file, verbose=verbose)
            df_to_test = df
        
        elif rad_scheme == "day_hours":
            # only day hours -- 7am-8pm as "day"
            printf('Radiation frequent value check scheme: day_hours selected, day set to 7am - 8pm', log_file=log_file, verbose=verbose)
            # 6am PST ~ 1400 UTC, 8pm PST ~ 0400 UTC
            df_to_test = df.loc[(df['hour'] >= 14) | (df['hour'] <=4)]
            
        elif rad_scheme == "remove_zeros":
            # remove all zeros -- may remove too many zeros, impact daytime cloudy conditions, regional (PNW)
            printf('Radiation frequent value check scheme: remove_zeros selected, may remove valid daytime (cloudy) conditions', log_file=log_file, verbose=verbose)
            df_to_test = df.loc[df[var] >= bin_s]
            
    else: # all other variables
        df_to_test = df
    
    # If df_to_test is empty, just skip the next part
    if len(df_to_test)==0:
        return df

    # all data/annual checks
    if data_group == 'all':
        bins = create_bins_frequent(df_to_test, var, bin_size=bin_s) 
        bar_counts, bins = np.histogram(df_to_test[var], bins=bins)
        flagged_bins = bins_to_flag(bar_counts, bins)
        
        # flag values in that bin as suspect
        if len(flagged_bins) != 0:
            for sus_bin in flagged_bins:
                # indicate as suspect bins
                    # DECISION: preliminary flag? and then remove if okay/reset to nan?
                df.loc[(df[var]>=sus_bin) & (df[var]<=sus_bin+1), 
                       var+'_eraqc'] = 100 # highlight for further review flag, either overwritten with real flag or removed in next step

    #============================================================================================================
       
    elif data_group == 'annual':
        for yr in df_to_test.year.unique():
            df_yr = df_to_test.loc[df_to_test['year'] == yr]
            if df_yr[var].isna().all() == True: # some vars will have nan years
                continue
            bins = create_bins_frequent(df_yr, var, bin_size=bin_s) # using 1 degC/hPa bin width
            bar_counts, bins = np.histogram(df_yr[var], bins=bins)
            flagged_bins = bins_to_flag(bar_counts, bins, bin_main_thresh=20, secondary_bin_main_thresh=10)
            
            if len(flagged_bins) != 0:
                printf('Flagging bin: {0}'.format(flagged_bins), log_file=log_file, verbose=verbose)

                for sus_bin in flagged_bins:
                    df.loc[(df['year']==yr) & (df[var]>=sus_bin) & (df[var]<=sus_bin+1), 
                           var+'_eraqc'] = 24 # see era_qaqc_flag_meanings.csv
    
    #============================================================================================================
    # seasonal checks require special handling
    elif data_group == 'seasonal_all':
        for szn in szns:
            df_szn = df_to_test.loc[(df_to_test['month']==szn[0]) | (df_to_test['month']==szn[1]) | (df_to_test['month']==szn[2])]
            if df_szn[var].isna().all() == True:
                continue
            bins = create_bins_frequent(df_szn, var, bin_size=bin_s) # using 1 degC/hPa bin width
            bar_counts, bins = np.histogram(df_szn[var], bins=bins)
            flagged_bins = bins_to_flag(bar_counts, bins, bin_main_thresh=20, secondary_bin_main_thresh=20)
            
            if len(flagged_bins) != 0:
                for sus_bin in flagged_bins:
                    df.loc[((df['month']==szn[0]) | (df['month']==szn[1]) | (df['month']==szn[2])) & 
                           (df[var]>=sus_bin) & (df[var]<=sus_bin+1),
                           var+'_eraqc'] = 100 # highlight for further review flag, either overwritten with real flag or removed in next step
                    
    #============================================================================================================
                
    elif data_group == 'seasonal_annual':        
        for yr in df_to_test.year.unique():
            for szn in szns:
                  
                # all seasons except winter
                if szn != [12,1,2]:
                    df_szn = df_to_test.loc[(df_to_test['year']==yr) & 
                                    ((df_to_test['month']==szn[0]) | (df_to_test['month']==szn[1]) | (df_to_test['month']==szn[2]))] 
                    
                    if df_szn[var].isna().all() == True: # some vars will have nan years
                        continue

                    if yr==df_szn.loc[df_szn.index[-1],'year']:
                        if len(df_szn)==0:
                            break # after last season in last year

                    bins = create_bins_frequent(df_szn, var, bin_size=bin_s) # using 1 degC/hPa bin width
                    bar_counts, bins = np.histogram(df_szn[var], bins=bins)
                    flagged_bins = bins_to_flag(bar_counts, bins, bin_main_thresh=15, secondary_bin_main_thresh=10)

                    if len(flagged_bins) != 0:
                        printf('Flagging bins: {0}'.format(flagged_bins), log_file=log_file, verbose=verbose)

                        for sus_bin in flagged_bins:
                            df.loc[(df['year']==yr) & 
                                  ((df['month']==szn[0]) | (df['month']==szn[1]) | (df['month']==szn[2])) &
                                   (df[var]>=sus_bin) & (df[var]<=sus_bin+1),
                                  var+'_eraqc'] = 25 # see era_qaqc_flag_meanings.csv

                # special handling for winter because of december
                else:
                    df_yr = df_to_test.loc[df_to_test['year'] == yr] # that year's jan, feb, and wrong dec            
                    df_jf = df_yr.loc[(df_yr['month']==1) | (df_yr['month']==2)] # that specific year's jan and feb

                    df_d = df_to_test.loc[(df_to_test['year'] == yr-1) & (df_to_test['month'] == 12)] # previous year's dec
                    if len(df_d) == 0: # catching very first year instance
                        df_djf = df_jf 
                        printf('Winter season: proceeding with just Jan/Feb, no previous Dec', log_file=log_file, verbose=verbose) ## DECISION

                    else:
                        printf('Winter season: concatenating previous Dec', log_file=log_file, verbose=verbose)
                        df_djf = pd.concat([df_d, df_jf])
                        
                    if df_djf[var].isna().all() == True: # some vars will have nan years
                        continue
                                        
                    bins = create_bins_frequent(df_djf, var, bin_size=bin_s) # using 1 degC/hPa bin width
                    bar_counts, bins = np.histogram(df_djf[var], bins=bins)
                    flagged_bins = bins_to_flag(bar_counts, bins, bin_main_thresh=15, secondary_bin_main_thresh=10)

                    if len(flagged_bins) != 0:
                        printf('Flagging bins: {0}'.format(flagged_bins), log_file=log_file, verbose=verbose)

                        for sus_bin in flagged_bins:
                            # flag jan feb
                            df.loc[(df['year']==yr) & 
                                   ((df['month']==szn[1]) | (df['month']==szn[2])) &
                                   ((df[var]>=sus_bin) & (df[var]<=sus_bin+1)),
                                  var+'_eraqc'] = 25 # see era_qaqc_flag_meanings.csv
                            # flag correct dec
                            df.loc[((df['year']==yr-1) & (df['month']==szn[0])) &
                                   ((df[var]>=sus_bin) & (df[var]<=sus_bin+1)),
                                   var+'_eraqc'] = 25 # see era_qaqc_flag_meanings.csv

    return df

#-----------------------------------------------------------------------------
def synergistic_flag(df, num_temp_vars):  
    '''
    In frequent values, if air temp is flagged, dew point is also flagged, and vice versa.
    Applies appropriate flag in corresponding vars

    Input:
    -----
        df [pd.DataFrame]: station dataset converted to dataframe through QAQC pipeline
        num_temp_vars [list]: list of temperature vars

    Returns:
    -------
        df [pd.DataFrame]: QAQC dataframe with flagged values (see below for flag meaning)
    '''

    # need to identify which flag is placed
    # 24 for all obs/years check
    # 25 for all seasons/years check
    flags_to_set = [24, 25]
    
    for flag_to_set in flags_to_set:

        if 'tas' in num_temp_vars and 'tdps' in num_temp_vars:
            df.loc[df['tas_eraqc'] == flag_to_set, 'tdps_eraqc'] = flag_to_set
            df.loc[df['tdps_eraqc'] == flag_to_set, 'tas_eraqc'] = flag_to_set

        if 'tas' in num_temp_vars and 'tdps_derived' in num_temp_vars:
            df.loc[df['tas_eraqc'] == flag_to_set, 'tdps_derived_eraqc'] = flag_to_set
            df.loc[df['tdps_derived_eraqc'] == flag_to_set, 'tas_eraqc'] = flag_to_set    
            
    return df

#-----------------------------------------------------------------------------
def bins_to_flag(bar_counts, bins, bin_main_thresh=30, secondary_bin_main_thresh=30):
    '''Returns the specific bins to flag as suspect'''
    bins_to_flag = [] # list of bins that will be flagged
        
    for i in range(0, len(bar_counts)):
        # identify main bin + 3 on either side
        bin_end = i+4

        # need handling for first 3 blocks as there is no front
        if i < 3:
            bin_start = 0
        else:
            bin_start = i-3

        bin_block_sum = bar_counts[bin_start:bin_end].sum() # num of obs in the 7-bin block
        bin_main_sum = bar_counts[i] # num of obs in main bin

        # determine whether main bin is more than half sum in 7-block bin
        bin_block_50 = bin_block_sum * 0.5 # primary check at 50%
        bin_block_90 = bin_block_sum * 0.9 # secondary check at 90%
        
        if (bin_main_sum > bin_block_50) == True: 
            # ensure that bin_main_sum is greater than bin_main_thresh
            if bin_main_sum > bin_main_thresh:
                bins_to_flag.append(math.floor(bins[i]))
                
                # annual/seasonal check
                if (bin_main_sum > bin_block_90) == True:
                    if bin_main_sum > secondary_bin_main_thresh:
                        bins_to_flag.append(math.floor(bins[i])) 
                
            else: # less than bin_main_thresh obs in bin_main_sum, do not indicate as suspect
                continue
                                
    return bins_to_flag # returns a list of values that are suspicious
