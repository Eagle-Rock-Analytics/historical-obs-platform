"""
This is a script where Stage 3: QA/QC function(s) on unusually frequent values in the data observations are flagged. 
For use within the PIR-19-006 Historical Obsevations Platform.
"""

## Import Libraries
import numpy as np

#-----------------------------------------------------------------------------
## logic check: dew point must not exceed air temperature
def qaqc_crossvar_logic_tdps_to_tas(df, verbose=True):
    """
    Checks that dewpoint temperature does not exceed air temperature.
    If fails, only dewpoint temperature is flagged.
    """ 
    try:
        # First check that tdps and/or tdps_derived are provided
        dew_vars = [col for col in df.columns if 'tdps' in col]
        all_dew_vars = [var for var in dew_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc

        # dew point is not present
        if not all_dew_vars:
            if verbose:
                print('station does not report dew point temperature - bypassing temperature cross-variable logic check')
        # dew point is present
        else:
            for var in all_dew_vars: 
                isBad = df[var] > df['tas']
                df.loc[isBad, var + '_eraqc'] = 12 # see qaqc_flag_meanings.csv
                if verbose:
                    print('{0} eraqc flags (any other value than nan is an active flag!): {1}'.
                          format(var, df[var + '_eraqc'].unique()))
        return df
    
    except Exception as e:
        print("qaqc_crossvar_logic_tdps_to_tas failed with Exception: {}".format(e))
        return None    

#----------------------------------------------------------------------
## logic check: precip does not have any negative values
def qaqc_precip_logic_nonegvals(df, verbose=True):
    """
    Ensures that precipitation values are positive. Negative values are flagged as impossible.
    Provides handling for the multiple precipitation variables presently in the cleaned data. 
    """
    # pr_24h: Precipitation accumulated from last 24 hours
    # pr_localmid: Precipitation accumulated from local midnight
    # pr: Precipitation accumulated since last record
    # pr_1h: Precipitation accumulated in last hour
    # pr_5min: Precipitation accumulated in last 5 minutes
    
    # identify which precipitation vars are reported by a station
    all_pr_vars = [var for var in df.columns if 'pr' in var] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    pr_vars = [var for var in pr_vars if 'method' not in var]
    pr_vars = [var for var in pr_vars if 'duration' not in var]

    if not pr_vars: # precipitation variable(s) is not present
        print('station does not report precipitation - bypassing precip logic nonnegvals check')
        return None
    else:
        for item in pr_vars:
            if verbose:
                print('Precip range: ', df[item].min(), '-', df[item].max()) # testing
            isNeg = df[item] < 0
            df.loc[isNeg, item+'_eraqc'] = 10 # see era_qaqc_flag_meanings.csv

            if verbose:
                print('Precipitation eraqc flags (any other value than nan is an active flag!):' + 
                      '{}'.format(df[item+'_eraqc'].unique())) # testing
    return df

#----------------------------------------------------------------------
## logic check: precip accumulation amounts balance for time period
def qaqc_precip_logic_accum_amounts(df, verbose=True):
    """
    Ensures that precipitation accumulation amounts are consistent with reporting time frame.
    Only needs to be applied when 2 or more precipitation duration specific
    variables are present (pr_5min, pr_1h, pr_24h)
    For example: pr_5min should not be larger than pr_1h
    
    # pr: Precipitation accumulated since last record
    # pr_5min: Precipitation accumulated in last 5 minutes
    # pr_1h: Precipitation accumulated in last hour
    # pr_24h: Precipitation accumulated from last 24 hours
    # pr_localmid: Precipitation accumulated from local midnight
        
    # rules
    # pr_5min < pr_1h < pr_24h
    # pr_localmid should never exceed pr_24h
    """

    # identify which precipitation vars are reported by a station
    all_pr_vars = [var for var in df.columns if 'pr' in var] # can be variable length depending if there is a raw qc var
    pr_vars = [var for var in all_pr_vars if 'qc' not in var] # remove all qc variables so they do not also run through: raw, eraqc, qaqc_process
    pr_vars = [var for var in pr_vars if 'method' not in var]
    pr_vars = [var for var in pr_vars if 'duration' not in var]

    if not pr_vars: # precipitation variable(s) is not present
        print('station does not report precipitation - bypassing precip logic accum check')
        return None
    
    # if station does not report any precipitation values, or only one, bypass
    if len(pr_vars) == 0 or len(pr_vars) == 1:
        return df

    # checks accumulated precip vars against each other
    # noting that these flags are essentially identical in operation
    # flag assignment is logically dependent on the first var to determine 
    # which flag is placed (i.e. to determine if too larges/small)

    # checks accumulated precip vars against each other
    # noting that these flags are essentially identical in operation
    # flag assignment is logically dependent on the first var to determine which flag is placed (i.e. to determine if too larges/small)
    
    #:::::
    if 'pr_5min' in pr_vars:
        if 'pr_1h' in pr_vars:
            isBad = df['pr_5min'] > df['pr_1h']
            df.loc[isBad, 'pr_5min_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
            
        if 'pr_24h' in pr_vars:
            isBad = df['pr_5min'] > df['pr_24h']
            df.loc[isBad, 'pr_5min_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
        if verbose:
            print('Precip 5min eraqc flags (any other value than nan is an active flag!):' + 
                  '{}'.format(df['pr_5min_eraqc'].unique())) # testing

    #:::::
    if 'pr_1h' in pr_vars:
        if 'pr_5min' in pr_vars:
            isBad = df['pr_1h'] < df['pr_5min']
            df.loc[isBad, 'pr_1h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv
            
        if 'pr_24h' in pr_vars:
            isBad = df['pr_1h'] > df['pr_24h']
            df.loc[isBad, 'pr_1h_eraqc'] = 15 # see era_qaqc_flag_meanings.csv
        if verbose:
            print('Precip 1h eraqc flags (any other value than nan is an active flag!):' + 
                  '{}'.format(df['pr_1h_eraqc'].unique())) # testing

    #:::::
    if 'pr_24h' in pr_vars:
        if 'pr_5min' in pr_vars:
            isBad = df['pr_24h'] < df['pr_5min']
            df.loc[isBad, 'pr_24h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv
        
        if 'pr_1h' in pr_vars:
            isBad = df['pr_24h'] < df['pr_1h']
            df.loc[isBad, 'pr_24h_eraqc'] = 16 # see era_qaqc_flag_meanings.csv        
        
        if 'pr_localmid' in pr_vars:
            isBad = df['pr_24h'] < df['pr_localmid']
            df.loc[isBad, 'pr_24h_eraqc'] = 17 # see era_qaqc_flag_meanings.csv
            
        if verbose:
            print('Precip 24h eraqc flags (any other value than nan is an active flag!):' + 
                  '{}'.format(np.unique(df['pr_24h_eraqc']))) # testing

    return df

#----------------------------------------------------------------------
## logic check: wind direction must be 0 if wind speed is 0
def qaqc_crossvar_logic_calm_wind_dir(df, verbose=True):
    """
    Checks that wind direction is zero when wind speed is also zero.
    If fails, wind direction is flagged. # only flag wind direction?
    """
    try:
        # Noting that a wind direction value of 0 is a valid value
        # Only a problem when wind speed is also 0, where 0 now means no winds for there to be a direction

        # First check that wind direction is provided
        if 'sfcWind_dir' not in df.columns:
            if verbose:
                print('station does not report wind direction - bypassing wind cross-variable logic check')
                return df
            
        # First, identify calm winds but with incorrect wind directions
        isCalm = df['sfcWind'] == 0
        isDirNotZero = df['sfcWind_dir'] != 0
        isNotNan = ~df['sfcWind_dir'].isnull()
        isBad = isCalm & isDirNotZero & isNotNan
        
        df.loc[isBad, 'sfcWind_dir_eraqc'] = 13 # see qaqc_flag_meanings.csv
        
        # Next, identify non-zero winds but with incorrect wind directions
        # Non-zero northerly winds should be coded as 360 deg, not 0 deg
        isNotCalm = df['sfcWind'] != 0
        isDirZero = df['sfcWind_dir'] == 0
        isBad = isNotCalm & isDirZero & isNotNan
        
        df.loc[isBad, 'sfcWind_dir_eraqc'] = 14 # see qaqc_flag_meanings.csv
        
        if verbose:
            print('sfcWind_dir eraqc flags (any value other than nan is an active flag!): {}'.
                  format(df['sfcWind_dir_eraqc'].unique()))
        return df
        
    except Exception as e:
        print("qaqc_crossvar_logic_calm_wind_dir failed with Exception: {}".format(e))
        return None