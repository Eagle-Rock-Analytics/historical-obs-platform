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

## Set AWS credentials
s3 = boto3.resource("s3")
s3_cl = boto3.client('s3') # for lower-level processes

## Set relative paths to other folders and objects in repository.
bucket_name = "wecc-historical-wx"
wecc_terr = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_land.shp"
wecc_mar = "s3://wecc-historical-wx/0_maps/WECC_Informational_MarineCoastal_Boundary_marine.shp"

#======================================================================
# QA/QC Helper functions
def get_file_paths(network):
    rawdir = "1_raw_wx/{}/".format(network)
    cleandir = "2_clean_wx/{}/".format(network)
    qaqcdir = "3_qaqc_wx/{}/".format(network)
    mergedir = "4_merge_wx/{}/".format(network)
    return rawdir, cleandir, qaqcdir, mergedir

#----------------------------------------------------------------------
def get_wecc_poly(terrpath, marpath):
    """
    Identifies a bbox of WECC area to filter stations against
    Input vars: shapefiles for maritime and terrestrial WECC boundaries
    Returns: spatial objects for each shapefile, and bounding box for their union.
    """
    t = gp.read_file(terrpath)  ## Read in terrestrial WECC shapefile.
    m = gp.read_file(marpath)   ## Read in marine WECC shapefile.
    bbox = t.union(m).bounds    ## Combine polygons and get bounding box of union.
    return t,m, bbox

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
    if var == "tdps_derived":
        var = "tas"
    
    # Extract var data   
    data = df.copy()[var]
    
    # Convert from Pa to hPa
    if data.mean()>10000 and (var=="ps" or var=="slp"):
        data = data/100
    
    # Calculate modified mode from avg mode and median
    mode0 = data.sort_values().diff().replace(0, pd.NaT).mode().values[0]
    mode1 = data.sort_values().diff().replace(0, pd.NaT).median()
    mode = (mode0+mode1)/2
    
    # Round to the nearest 0.5
    multiplied = round(mode * 2)
    rounded_to_whole = multiplied / 2
    
    # If mode is 0.25 or less, round to 0.1
    if rounded_to_whole == 0:
        mode = 0.1
    else:
        mode = rounded_to_whole
    return mode

#----------------------------------------------------------------------
def infere_res(df):
    """
    """
    variables = ["tas", "tdps", "tdps_derived", "ps", "slp", "sfcWind"]
    resolutions = {}

    for var in variables:
        if var in df.columns:
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
                            "slp" : {1   : [120, 28],  # of
                                     0.5 : [100, 21],  # or
                                     0.1 : [ 72, 14]}, # or
         
                            "sfcWind" : {1   : [40, 14], # of
                                         0.5 : [30, 10], # or
                                         0.1 : [24,  7], # or
                                        },
                           }
straight_repeat_criteria['tdps_derived'] = straight_repeat_criteria['tdps']
straight_repeat_criteria['ps'] = straight_repeat_criteria['slp']
straight_repeat_criteria['sfcWind'] = straight_repeat_criteria['tas']

#----------------------------------------------------------------------
# Hour repeat streak criteria
hour_repeat_criteria = {"tas" : {1   : 25,  # 40 days
                                 0.5 : 20,  # 20 days
                                 0.1 : 15,  # 15 days
                                }
                       }
hour_repeat_criteria['tdps'] = hour_repeat_criteria['tas']
hour_repeat_criteria['tdps_derived'] = hour_repeat_criteria['tdps']
hour_repeat_criteria['slp'] = hour_repeat_criteria['tas']
hour_repeat_criteria['ps'] = hour_repeat_criteria['slp']
hour_repeat_criteria['sfcWind'] = hour_repeat_criteria['tas']

#----------------------------------------------------------------------
# Day repeat streak criteria
day_repeat_criteria = {"tas" : {1   : 10,  # 10 days
                                0.5 :  7,  #  7 days
                                0.1 :  5,  #  5 days
                               }
                       }
day_repeat_criteria['tdps'] = day_repeat_criteria['tas']
day_repeat_criteria['tdps_derived'] = day_repeat_criteria['tdps']
day_repeat_criteria['slp'] = day_repeat_criteria['tas']
day_repeat_criteria['ps'] = day_repeat_criteria['slp']
day_repeat_criteria['sfcWind'] = day_repeat_criteria['tas']

#----------------------------------------------------------------------
# Min wind value for straight repeat test
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
def qaqc_unusual_repeated_streaks(df, plot=True, local=False, verbose=True, min_sequence_length=10):
    """
    Test for repeated streaks/unusual spell frequenc. 
    Three test are conducted here:
       - Consecutive observation replication
       - Same hour observation replication over a number of days 
         (either using a threshold of a certain number of observations, 
          or for sparser records, a number of days during which all the
          observations have the same value)
       - Whole day replication for a streak of days
    
    This test is done for ["tas", "tdps", "tdps_derived", "ps", "slp", "sfcWind"]
    
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
        28,qaqc_unusual_repeated_streaks,
        29,qaqc_unusual_repeated_streaks,

    NOTES:
    Threshold for different variables/resolutions are noted on: 
    https://doi.org/10.5194/cp-8-1649-2012 : Table 4
    
    TODO: 
    min_sequence_length can be tweaked, althought HadISD uses 10
    """

    try:
    # if True:
        
        # Infere resolution from data
        resolutions = infere_res(df)
        
        # Save original df multiindex and create station column
        df = df.copy(deep=True)

        # Define test variables and check if they are in the dataframe
        check_vars = ["tas", "tdps", "tdps_derived", "ps", "slp", "sfcWind"]
        # check_vars = ["ps"]
        variables = [var for var in check_vars if var in df.columns]
        
        if verbose:
            print("Running {} on {}".format("qaqc_unusual_repeated_streaks", variables))
        
        # Loop through test variables
        for var in variables:
            print(var)
            # Create a copy of the original dataframe and drop NaNs in the testing variable
            new_df = df.copy(deep=True)
            new_df = new_df.dropna(subset=var)#.drop(columns=["lat","lon","elevation"])
            
            # Use only values that have not been flagged by previous QAQC tests
            valid = np.where(np.isnan(new_df[var+"_eraqc"]))[0]
            new_df = new_df.iloc[valid]
            
            # Choose resolution
            res = resolutions[var]
            
            # --------------------------------------------------------
            # Hour repeat streak criteria
            threshold = hour_repeat_criteria[var][res]
            bad_hourly = hourly_repeats(new_df, var, threshold)            
            ind = df['time'].isin(bad_hourly['time']) 
            df.loc[ind, var+"_eraqc"] = 27    # Flag _eraqc variable

            # --------------------------------------------------------
            # Straight repeat streak criteria
            threshold = straight_repeat_criteria[var][res]
            if var=="sfcWind":
                wind_min_value = WIND_MIN_VALUE[res]
            else:
                wind_min_value = None
            bad_straight = consecutive_repeats(new_df, var, threshold, 
                                               wind_min_value, 
                                               min_sequence_length=min_sequence_length)
            ind = df['time'].isin(bad_straight['time']) 
            df.loc[ind, var+"_eraqc"] = 28    # Flag _eraqc variable
        
            # --------------------------------------------------------
            # Whole day replication for a streak of days
            threshold = day_repeat_criteria[var][res]
            bad_whole = consecutive_fullDay_repeats(new_df, var, threshold)
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
            # import pdb; pdb.set_trace()
            bad['consecutive_month_group'] = bad.copy().groupby('year')['month'].transform(consecutive_months)
            bad_min = bad.groupby(by=["year","consecutive_month_group"])['time'].min()
            bad_max = bad.groupby(by=["year","consecutive_month_group"])['time'].max()
            bad = pd.DataFrame(data = {"min_date":bad_min.values,
                                       "max_date":bad_max.values})
            
            # --------------------------------------------------------
            if plot:
                unusual_streaks_plot(df, var, local=local)
                
                for i in bad.index:
                    
                    da = bad.loc[i]
                    
                    min_date = da.min_date - np.timedelta64(3,'D')
                    max_date = da.min_date + np.timedelta64(3,'D')
                    subset = np.logical_and(df['time'] >= min_date, 
                                            df['time'] <= max_date)
                    unusual_streaks_plot(df[subset], var, date=min_date+np.timedelta64(3,'D'), local=local)
        return df
    except Exception as e:
        print("qaqc_unusual_repeated_streaks failed with Exception: {}".format(e))
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
    
    This test is done for ["tas", "tdps", "tdps_derived", "ps", "slp", "sfcWind"]
    
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
    
    df = df.copy()
    df['hours'] = pd.Series(df['time']).dt.hour.values
    counts = pd.DataFrame(df.groupby(by=["hours",var]).apply(lambda x: np.array(x['time'].tolist())).rename("dates"))
    counts['date_diff'] = counts['dates'].transform(lambda x: pd.Series(x).diff().values.astype("timedelta64[D]"))
    counts['streak_index'] = counts['date_diff'].apply(find_streaks_index, args=(7,))
    counts['streaks'] = counts.apply(lambda x: select_streaks(x.dates, x.streak_index), axis=1)
    
    groups = counts[counts['streaks'].apply(len)>0]
    if len(groups)>0:
        bad = pd.DataFrame({"time"  : np.concatenate(groups['streaks'].values)})
        bad['month'] = bad.loc[:, 'time'].dt.month.values
        bad['year']  = bad.loc[:, 'time'].dt.year.values
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
    
    This test is done for ["tas", "tdps", "tdps_derived", "ps", "slp", "sfcWind"]
    
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
        bad['month'] = bad.loc[:, 'time'].dt.month.values
        bad['year']  = bad.loc[:, 'time'].dt.year.values
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
    
    This test is done for ["tas", "tdps", "tdps_derived", "ps", "slp", "sfcWind"]
    
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
    da['date'] = pd.Series(da['time']).dt.date.values
    
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
    
    df['date'] = pd.Series(df['time']).dt.date.values    
    bad = df.copy()[df['date'].isin(bad_dates['date'])]
    bad['group'] = [bad_dates['group'].loc[d] for d in bad['date']]

    if len(bad)>0:
        bad['month'] = bad.loc[:, 'time'].dt.month.values
        bad['year']  = bad.loc[:, 'time'].dt.year.values
    else:
        bad = pd.DataFrame({"time"  : np.array([], dtype='datetime64[ns]'),
                            "month" : np.array([], dtype=np.int64),
                            "year"  : np.array([], dtype=np.int64)})    
    return bad[['time', 'month', 'year']]

#---------------------------------------------------------------------------------------------------
def unusual_streaks_plot(df, var, flagvals=(27,28,29), dpi=None, local=False, date=None):
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
    flag_vals_0 = df.loc[df[var + '_eraqc'] == flagvals[0]]   
    flag_vals_1 = df.loc[df[var + '_eraqc'] == flagvals[1]]   
    flag_vals_2 = df.loc[df[var + '_eraqc'] == flagvals[2]]   
    
    # Create figure
    if date is not None:
        fig,ax = plt.subplots(figsize=(7,3))
    else:
        fig,ax = plt.subplots(figsize=(10,3))

    # Plot variable and flagged data
    df.plot("time", var, ax=ax, marker=".", ms=4, lw=1, color="k", alpha=0.5, label="Original data")
    
    # Amount of data flagged
    nflags = len(flag_vals_0) + len(flag_vals_1) + len(flag_vals_2)
    title = "{:.4f}% of data flagged".format(100*nflags/len(df))
    flag_label_0 = "Same hour replication"
    flag_label_1 = "Consecutive replication"
    flag_label_2 = "Whole-day replication"
    df.loc[df[var+"_eraqc"]==flagvals[0]].plot("time", var, ax=ax, marker="s", ms=7, lw=0, mfc="none", color="C3", label=flag_label_0)    
    df.loc[df[var+"_eraqc"]==flagvals[1]].plot("time", var, ax=ax, marker="x", ms=7, lw=0, mfc="none", color="C4", label=flag_label_1)    
    df.loc[df[var+"_eraqc"]==flagvals[2]].plot("time", var, ax=ax, marker="o", ms=7, lw=0, mfc="none", color="C4", label=flag_label_2)    
    legend = ax.legend(loc=0, prop={'size': 8})    
    title = ax.set_title(title)    
        
    station = df['station'].unique()[0]
    network = station.split('_')[0]
    
    # Plot aesthetics
    ylab, units, miny, maxy = _plot_format_helper_spikes(var)
    ylab = '{} [{}]'.format(ylab, units)
    
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
    figname = 'qaqc_figs/qaqc_unusual_repeated_streaks_{0}_{1}_{2}'.format(station, var, timestamp)
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

#-----------------------------------------------------------------------------------------
# Plot helper for unusual large jumps plots
def _plot_format_helper_spikes(var):
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
    from calc_qaqc import _plot_format_helper_spikes

    pr_vars = ['pr', 'pr_5min', 'pr_1h', 'pr_24h', 'pr_localmid']
    ps_vars = ['ps', 'psl', 'psl_altimeter']
    
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



#########################################################################################################
#========================================================================================================
# TEMPORARY FUNCTIONS FOR TESTING

#---------------------------------------------------------------------------------------------------
def qaqc_world_record(df, verbose=True):
    '''
    Checks if temperature, dewpoint, windspeed, or sea level pressure are outside North American world records
    If outside minimum or maximum records, flags values
    '''
    try:
        # world records from HadISD protocol, cross-checked with WMO database
        # https://wmo.asu.edu/content/world-meteorological-organization-global-weather-climate-extremes-archive
        T_X = {"North_America":329.92} #K
        T_N = {"North_America":210.15} #K
        D_X = {"North_America":329.85} #K
        D_N = {"North_America":173.15} #K
        W_X = {"North_America":113.2} #m/s
        W_N = {"North_America":0.} #m/s
        S_X = {"North_America":108330} #Pa
        S_N = {"North_America":87000} #Pa

        maxes = {"tas": T_X, "tdps": D_X, "tdps_derived": D_X, "sfcWind": W_X, "psl": S_X, "ps": S_X}
        mins =  {"tas": T_N, "tdps": D_N, "tdps_derived": D_N, "sfcWind": W_N, "psl": S_N, "ps": S_N}

        # variable names to check against world record limits
        wr_vars = ['tas', 'tdps_derived', 'tdps', 'sfcWind', 'psl', 'ps']

        for var in wr_vars:
            if var in list(df.columns):
                                
                isOffRecord = np.logical_or(df[var] < mins[var]['North_America'],
                                            df[var] > maxes[var]['North_America'])
                if isOffRecord.any():
                    df.loc[isOffRecord, var + '_eraqc'] = 11
        return df
    except Exception as e:
        if verbose:
            print("qaqc_world_record failed with Exception: {}".format(e))
        return None

#---------------------------------------------------------------------------------------------------
def xarray_to_pandas_qaqc(ds):
    """
    """

    ## Add qc_flag variable for all variables, including elevation; 
    ## defaulting to nan for fill value that will be replaced with qc flag
    exclude_qaqc = ["time", "station", "lat", "lon", 
                    "qaqc_process", "sfcWind_method"] # lat and lon have a different qc check

    raw_qc_vars = [] # qc_variable for each data variable, will vary station to station
    era_qc_vars = [] # our qc variable
    for var in ds.data_vars:
        if 'q_code' in var:
            raw_qc_vars.append(var) # raw qc variable, need to keep for comparison, then drop
        if '_qc' in var:
            raw_qc_vars.append(var) # raw qc variables, need to keep for comparison, then drop

    for var in ds.data_vars:
        if var not in exclude_qaqc and var not in raw_qc_vars:
            qc_var = var + "_eraqc" # variable/column label
            era_qc_vars.append(qc_var)
            # adds new variable in shape of original variable with designated nan fill value
            ds = ds.assign({qc_var: xr.ones_like(ds[var])*np.nan})

    # Save attributes to inheret them to the QAQC'ed file
    attrs = ds.attrs
    var_attrs = {var:ds[var].attrs for var in list(ds.data_vars.keys())}

    df = ds.to_dataframe()
    df['anemometer_height_m'] = np.ones(ds['time'].shape)*ds.anemometer_height_m
    df['thermometer_height_m'] = np.ones(ds['time'].shape)*ds.thermometer_height_m
    
    # Save station/time multiindex
    MultiIndex = df.index
    station = df.index.get_level_values(0)
    df['station'] = station
    
    # Station pd.Series to str
    station = station.unique().values[0]
    
    # Convert time/station index to columns and reset index
    df = df.droplevel(0).reset_index()
    
    return df