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
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
try:
    from qaqc_plot import *
except:
    print("Error importing qaqc_plot.py")
## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

try:
    from qaqc_utils import *
except Exception as e:
    print("Error importing qaqc_utils: {}".format(e))
    
def open_log_file_streaks(file):
    global log_file
    log_file = file

#----------------------------------------------------------------------
def infere_freq(df):
    """
    """
    # Calculate time differences in data
    time_diff = pd.Series(df['time']).diff()

    # Count resolution occurences and use only those that occure in 5% of series or more
    frequencies = time_diff.value_counts()
    frequencies = frequencies[frequencies > len(df)//5]

    # Create resolutions dictionary
    frequencies = {
                   np.round(n/len(df), decimals=5):float(f/1e9) for 
                   n,f in zip(frequencies.values,
                             frequencies.index.values)
                  }
    return frequencies

#----------------------------------------------------------------------
def infere_res_var(df, var):
    """
    """
    # If var is tdps_derived, use temp measurements
    # if var == "tdps_derived":
    #     var = "tas"
    
    # Extract var data   
    data = df.copy()[var]
    
    #TODO: use function from calc_clean.py
    # Convert from Pa to hPa -- only to determine mode and resolution of data, data not returned in hPa
    if data.mean()>10000 and (var=="ps" or var=="psl" or var=='ps_altimeter' or var=='ps_derived'):
        data = data/100

    data_diff = data.sort_values().diff().clip(lower=0.01, upper=1).dropna()
    if len(data_diff)<=10:
        return 0.5
    else:
        # Calculate modified mode from avg mode and median
        # mode0 = data.sort_values().diff().replace(0, pd.NaT).mode().values[0]
        # mode1 = data.sort_values().diff().replace(0, pd.NaT).median()
        mode0 = data.sort_values().diff().mode().values[0]
        mode1 = data.sort_values().diff().median()
        mode = (mode0+mode1)/2

        # Round to the nearest 0.5
        multiplied = round(mode * 2)
        rounded_to_whole = multiplied / 2

        # If mode is 0.25 or less, round to 0.1
        if rounded_to_whole <= 0.25:
            mode = 0.1
        else:
            mode = rounded_to_whole

        if mode<=1:
            return mode
        else:
            return 1.0

#----------------------------------------------------------------------
def infere_res(df, verbose=False):
    """
    """
    # check_vars = ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter", "sfcWind"]
    check_vars = ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter", "sfcWind", "rsds"]
    variables = [var for var in check_vars if var in df.columns]
    # variables = [var for var in df.columns if any(True for item in check_vars if item in var)]
    # printf(variables, log_file=log_file, verbose=verbose)
    
    resolutions = {}
    for var in variables:
        # if var in df.columns:
        if df[var].isnull().all() or len(np.where(df[var].isnull())[0])<100:
            resolutions[var] = 0.1
        else:
            resolutions[var] = infere_res_var(df, var)

    return resolutions

#----------------------------------------------------------------------
# Straight repeat streak criteria
straight_repeat_criteria = {"tas" : {1  : [40, 14],   # 40 values of 14 days
                                     0.5 : [30, 10],  # 30 values or 10 days
                                     0.1 : [24,  7],  # 24 values or 7 days
                           },
                            "tdps" : {1  : [80, 14],  # of
                                      0.5 : [60, 10], # or
                                      0.1 : [48,  7], # or
                                             },
                            "psl" : {1   : [120, 28],  # of
                                     0.5 : [100, 21],  # or
                                     0.1 : [ 72, 14]}, # or
         
                            "sfcWind" : {1   : [40, 14], # of
                                         0.5 : [30, 10], # or
                                         0.1 : [24,  7], # or
                                        },
                           }
straight_repeat_criteria['tdps_derived'] = straight_repeat_criteria['tdps']
straight_repeat_criteria['ps'] = straight_repeat_criteria['psl']
straight_repeat_criteria['ps_derived'] = straight_repeat_criteria['psl']
straight_repeat_criteria['ps_altimeter'] = straight_repeat_criteria['psl']
straight_repeat_criteria['rsds'] = straight_repeat_criteria['psl']

#----------------------------------------------------------------------
# Hour repeat streak criteria
hour_repeat_criteria = {"tas" : {1   : 25,  # 40 days
                                 0.5 : 20,  # 20 days
                                 0.1 : 15,  # 15 days
                                }
                       }
# All variables have the same hourly criteria
hour_repeat_criteria['tdps'] = hour_repeat_criteria['tas']
hour_repeat_criteria['tdps_derived'] = hour_repeat_criteria['tas']
hour_repeat_criteria['psl'] = hour_repeat_criteria['tas']
hour_repeat_criteria['ps'] = hour_repeat_criteria['tas']
hour_repeat_criteria['ps_altimeter'] = hour_repeat_criteria['tas']
hour_repeat_criteria['ps_derived'] = hour_repeat_criteria['tas']
hour_repeat_criteria['sfcWind'] = hour_repeat_criteria['tas']
hour_repeat_criteria['rsds'] = hour_repeat_criteria['sfcWind']

#----------------------------------------------------------------------
# Day repeat streak criteria
day_repeat_criteria = {"tas" : {1   : 10,  # 10 days
                                0.5 :  7,  #  7 days
                                0.1 :  5,  #  5 days
                               }
                       }
# All variables have the same daily criteria
day_repeat_criteria['tdps'] = day_repeat_criteria['tas']
day_repeat_criteria['tdps_derived'] = day_repeat_criteria['tas']
day_repeat_criteria['psl'] = day_repeat_criteria['tas']
day_repeat_criteria['ps'] = day_repeat_criteria['tas']
day_repeat_criteria['ps_altimeter'] = day_repeat_criteria['tas']
day_repeat_criteria['ps_derived'] = day_repeat_criteria['tas']
day_repeat_criteria['sfcWind'] = day_repeat_criteria['tas']
day_repeat_criteria['rsds'] = day_repeat_criteria['sfcWind']

#----------------------------------------------------------------------
# Min wind value for straight repeat test
# TODO: HadISD thresholds: does it make sense to change them in future versions?
# More analysis needs to be done to ensure what is a good threshold for calm wind conditions for this test
WIND_MIN_VALUE = {1:1.0, 0.5:0.5, 0.1:0.5}

#---------------------------------------------------------------------------------------------------
# Function to create a new column for consecutive months
def consecutive_months(series):
    
    indices = np.where(np.diff(series.values) > 1)[0] + 1
    clusters = np.split(series.values, indices)
    isin = [series.isin(c) for c in clusters]
    groups = np.zeros_like(series.values, dtype="int")
    
    for i,ind in enumerate(isin):
        groups[ind.values] = int(i)
    return pd.Series(groups, index=series.index)

#---------------------------------------------------------------------------------------------------
def qaqc_unusual_repeated_streaks(df, plot=False, local=False, verbose=False, min_sequence_length=10):
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

    try:
        
        # Infere resolution from data
        resolutions = infere_res(df)
        
        # Save original df multiindex and create station column
        new_df = df.copy()
        new_df['hours'] = pd.to_datetime(df['time']).dt.hour.values
        new_df['day']   = pd.to_datetime(df['time']).dt.day.values
        new_df['month'] = pd.to_datetime(df['time']).dt.month.values
        new_df['year']  = pd.to_datetime(df['time']).dt.year.values
        new_df['date']  = pd.to_datetime(df['time']).dt.date.values
        
        # Define test variables and check if they are in the dataframe
        check_vars = ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter", "sfcWind", "rsds"]
        # check_vars = ["ps"]
        variables = [var for var in check_vars if var in new_df.columns]
        printf("Running {} on {}".format("qaqc_unusual_repeated_streaks", variables), verbose=verbose, log_file=log_file)
        
        # Loop through test variables
        for var in variables:
            
            # Create a copy of the original dataframe and drop NaNs in the testing variable
            test_df = new_df.copy().dropna(subset=var)
        
            # Use only values that have not been flagged by previous QAQC tests
            valid = np.where(np.isnan(test_df[var+"_eraqc"]))[0]
            test_df = test_df.iloc[valid]
            
            # first scans suspect values using entire record
            if test_df[var].isna().all() == True:
                printf("All values for {} are flagged, bypassing qaqc_unusual_repeated_streaks".format(var), verbose=verbose, log_file=log_file)
                continue # bypass to next variable if all obs are nans
                
            # Choose resolution
            res = resolutions[var]
            
            # --------------------------------------------------------
            # Hour repeat streak criteria
            threshold = hour_repeat_criteria[var][res]
            bad_hourly = hourly_repeats(test_df, var, threshold)            
            ind = df['time'].isin(bad_hourly['time']) 
            df.loc[ind, var+"_eraqc"] = 27    # Flag _eraqc variable

            # --------------------------------------------------------
            # Straight repeat streak criteria
            threshold = straight_repeat_criteria[var][res]
            if var=="sfcWind":
                wind_min_value = WIND_MIN_VALUE[res]
            else:
                wind_min_value = None
            bad_straight = consecutive_repeats(test_df, var, threshold, 
                                               wind_min_value, 
                                               min_sequence_length=min_sequence_length)
            ind = df['time'].isin(bad_straight['time']) 
            df.loc[ind, var+"_eraqc"] = 28    # Flag _eraqc variable
        
            # --------------------------------------------------------
            # Whole day replication for a streak of days
            threshold = day_repeat_criteria[var][res]
            bad_whole = consecutive_fullDay_repeats(test_df, var, threshold)
            ind = df['time'].isin(bad_whole['time']) 
            df.loc[ind, var+"_eraqc"] = 29    # Flag _eraqc variable
           
            # --------------------------------------------------------
            # Groups for zoom plots
            bad_months = np.concatenate((bad_hourly['month'].values,
                                         bad_straight['month'].values,
                                         bad_whole['month'].values))
            bad_years = np.concatenate((bad_hourly['year'].values,
                                        bad_straight['year'].values,
                                        bad_whole['year'].values))
            bad_times = np.concatenate((bad_hourly['time'].values,
                                        bad_straight['time'].values,
                                        bad_whole['time'].values))
            
            bad = pd.DataFrame({"year":bad_years, "month":bad_months, "time":bad_times})
            
            bad['consecutive_month_group'] = bad.copy().groupby('year')['month'].transform(consecutive_months)
            bad_min = bad.groupby(by=["year","consecutive_month_group"])['time'].min()
            bad_max = bad.groupby(by=["year","consecutive_month_group"])['time'].max()
            bad = pd.DataFrame(data = {"min_date":bad_min.values,
                                       "max_date":bad_max.values})
        
            # --------------------------------------------------------
            if plot:
                printf("Full timeseries plot", verbose=verbose, log_file=log_file)
                unusual_streaks_plot(df, var, local=local)
                
                for i in bad.index:
                    printf("Subset plots", verbose=verbose, log_file=log_file)
                    da = bad.loc[i]
                    
                    min_date = da.min_date - np.timedelta64(3,'D')
                    max_date = da.min_date + np.timedelta64(3,'D')
                    subset = df.loc[(df['time'] >= min_date) and (df['time'] <= max_date)]
                    # subset = np.logical_and(df['time'] >= min_date, 
                    #                         df['time'] <= max_date)
                    unusual_streaks_plot(subset, var, date=min_date+np.timedelta64(3,'D'), local=local)
        
        return df
    except Exception as e:
        printf("qaqc_unusual_repeated_streaks failed with Exception: {}".format(e), verbose=verbose, log_file=log_file)
        return None

#---------------------------------------------------------------------------------------------------
# Find clusters of equal values
def select_streaks(x,y):
    """
    """
    if len(y)>0:
        return x[y]
    else:
        return []
    
#---------------------------------------------------------------------------------------------------
# Find clusters of equal values
def find_streaks_index(y, threshold=5):
    """
    """
    y = np.concatenate([y])
    y = y.astype("int")
    y[np.where(y==0)[0]] = 1
    split_indices = np.where(np.diff(y) != 0)[0] + 1
    clusters = np.split(y, split_indices)

    ind = [list(np.arange(split_indices[i-1], split_indices[i-1]+len(clusters[i]), 1)) 
           for i in range(1,len(clusters))]
    result = [i for i in ind if len(i)>threshold]
    
    if len(result)>0:
        return np.concatenate(result)
    else:
        return np.array([])
    
#---------------------------------------------------------------------------------------------------
def hourly_repeats(df, var, threshold):
    """
    Same hour observation replication over a number of days 
         (either using a threshold of a certain number of observations, 
          or for sparser records, a number of days during which all the
          observations have the same value):
          
    1 - Group variable by hour/value
    2 - Count consecutive cluster of repeated values at the same hour of the day
    3 - Select only clusters where count is higher than threshold in dict `hour_repeat_criteria`
    
    This test is done for ["tas", "tdps", "tdps_derived", "ps", "psl", "ps_derived", "ps_altimeter", "sfcWind"]
    
    Input:
    -----
        df [pandas dataframe] : station dataset converted to dataframe through QAQC pipeline
        var [str] : variable to test
        threshold [int] : comes from hour_repeat_criteria[var][res]
   Output:
    ------
        bad [numpy array] : dates that mark the flagged values (from df.index)
                  
    NOTES (TODO:)
    """
    
    da = df.copy()
    da['hours'] = pd.to_datetime(da['time']).dt.hour.values
    counts = pd.DataFrame(da.groupby(by=["hours",var], group_keys=True).apply(lambda x: np.array(x['time'].tolist())).rename("dates"))
    # counts = da.groupby(by=["hours",var], group_keys=False).apply(lambda x: np.array(x['time'].tolist()))#.rename("dates")
    counts['date_diff'] = counts['dates'].transform(lambda x: pd.Series(x).diff().values.astype("timedelta64[D]"))
    counts['streak_index'] = counts['date_diff'].apply(find_streaks_index, args=(7,))
    counts['streaks'] = counts.apply(lambda x: select_streaks(x.dates, x.streak_index), axis=1)
    
    groups = counts[counts['streaks'].apply(len)>0]
    if len(groups)>0:
        bad = pd.DataFrame({"time"  : np.concatenate(groups['streaks'].values)})
        bad['month'] = pd.to_datetime(bad.loc[:, 'time']).dt.month.values
        bad['year']  = pd.to_datetime(bad.loc[:, 'time']).dt.year.values
    else:
        bad = pd.DataFrame({"time"  : np.array([], dtype='datetime64[ns]'),
                            "month" : np.array([], dtype=np.int64),
                            "year"  : np.array([], dtype=np.int64)})  
    return bad

#---------------------------------------------------------------------------------------------------
def consecutive_repeats(df, var, threshold, wind_min_value = None, 
                        min_sequence_length = 10):
    """
    Consecutive observation replication 
         (either using a threshold of a certain number of observations, 
          or for sparser records, a number of days during which all the
          observations have the same value):
          
    1 - 
    2 - 
    3 - 
    
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

    da = df.copy()[[var,"time"]]
    
    # If variable is wind, only use values above min wind value
    if wind_min_value is not None:
        da = da[da[var]>wind_min_value]
    # Identify sequences of similar values
    da.loc[:, 'group'] = (da[var] != da[var].shift()).cumsum()
    # da.loc[:, 'start_date'] = da.index.values
    # da.loc[:, 'end_date'] = da.index.values
    da.loc[:, 'start_date'] = da['time'].values
    da.loc[:, 'end_date'] = da['time'].values
    
    start_date = da.copy().groupby([var, 'group']).min().sort_values(by="group").copy()
    end_date   = da.copy().groupby([var, 'group']).max().sort_values(by="group").copy()
    
    # Calculate the length of each sequence
    sequence_lengths = da.copy().groupby([var, 'group']).size().reset_index(name='sequence_length').sort_values(by="group")
    sequence_lengths.loc[:, 'start_date'] = start_date.loc[:, 'start_date'].values
    sequence_lengths.loc[:, 'end_date'] = end_date.loc[:, 'end_date'].values
    
    # Filter sequences with a minimum length
    min_sequence_length = 10
    filtered_sequences = sequence_lengths.copy()[sequence_lengths.loc[:, 'sequence_length'] >= min_sequence_length]
    filtered_sequences.loc[:, 'dt'] = (filtered_sequences.loc[:, 'end_date'] -
                                filtered_sequences.loc[:, 'start_date']).values.astype("timedelta64[D]")
    
    # Get limits
    nvalues, ndays = threshold
    
    #TODO: fix this commented line, 2x sequence_length to test, it should be dt
    condition = np.logical_or(filtered_sequences.loc[:, 'sequence_length'] > nvalues,
                              filtered_sequences['dt'].astype("int") > ndays)
    
    # Find bad groups and index in the original dataset
    bad_groups = filtered_sequences[condition].group.values
    bad = da[da['group'].isin(bad_groups)].copy()
    
    if len(bad)>0:
        bad['month'] = pd.to_datetime(bad.loc[:, 'time']).dt.month.values
        bad['year']  = pd.to_datetime(bad.loc[:, 'time']).dt.year.values
    else:
        bad = pd.DataFrame({"time"  : np.array([], dtype='datetime64[ns]'),
                            "month" : np.array([], dtype=np.int64),
                            "year"  : np.array([], dtype=np.int64)})    
    return bad[['time', 'month', 'year']]

#---------------------------------------------------------------------------------------------------
def full_day_compare(series0,series1):
    
    ind = []
    groups = []
    g = 0
    for a,b in zip(series0.values, series1.values):
        if type(a)==np.ndarray and type(b)==np.ndarray and len(a)==len(b):
            if ((a==b).all()):
                groups.append(g)
            else:
                groups.append(-1)
                g += 1
        else:
            groups.append(-1)
            g += 1
    if len(np.unique(groups))>1:
        return np.array(groups)-np.max(groups)
    else:
        return groups

#---------------------------------------------------------------------------------------------------
def consecutive_fullDay_repeats(df, var, threshold):
    """
    Consecutive observation replication 
         (either using a threshold of a certain number of observations, 
          or for sparser records, a number of days during which all the
          observations have the same value):
          
    1 - 
    2 - 
    3 - 
    
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
    da = df.copy()[[var, 'time','year','month','day', 'hours']].dropna()
    datavar = da.groupby(by=['year','month','day','hours'])[var].apply(np.nanmean)
    othervars = da.drop(columns=var).groupby(by=['year','month','day','hours']).first()
    data = {"hours":othervars.index.get_level_values(-1), 
            var:datavar.values,
            "time":othervars.time.values}
    da = pd.DataFrame(data)
    da['date'] = pd.to_datetime(da['time']).dt.date.values
    
    # Whole days to analysis
    whole_days = da.groupby(by=['date'])[var].apply(lambda x: np.round(x.values, decimals=1))
    whole_days = pd.DataFrame({var:whole_days, 'date':whole_days.index.values})
    whole_days['group'] = full_day_compare(whole_days[var], whole_days[var].shift())
    
    sequence_lengths = whole_days.groupby(['group']).size().reset_index(name='length').sort_values(by="group")
    sequence_lengths = sequence_lengths[sequence_lengths['group']>=0]
    bad_groups = np.array([g for g,l in 
                           zip(sequence_lengths['group'].values, 
                               sequence_lengths['length'].values) 
                           if l>threshold])
    
    bad_dates = whole_days[whole_days['group'].isin(bad_groups)][['date','group']]
    
    df['date'] = pd.to_datetime(df['time']).dt.date.values    
    bad = df.copy()[df['date'].isin(bad_dates['date'])]
    bad['group'] = [bad_dates['group'].loc[d] for d in bad['date']]

    if len(bad)>0:
        bad['month'] = pd.to_datetime(bad.loc[:, 'time']).dt.month.values
        bad['year']  = pd.to_datetime(bad.loc[:, 'time']).dt.year.values
    else:
        bad = pd.DataFrame({"time"  : np.array([], dtype='datetime64[ns]'),
                            "month" : np.array([], dtype=np.int64),
                            "year"  : np.array([], dtype=np.int64)})    
    return bad[['time', 'month', 'year']]
